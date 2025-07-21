use strict;
use warnings;
use Cwd;
use Getopt::Long;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);
use Bio::SeqIO;
use utils qw(make_dir 
            construct_paired_filename
            append_assemblies);

exit main();

sub main {
    my $config = parse_arguments();
    my @commands = write_commands($config);

    my $num_parallel = $config->{params}{num_selections};
    my @pids;

    # Step 1: Run even-indexed commands in parallel
    for (my $i = 0; $i < ($num_parallel + 1) * 2; $i += 2) {
        last if $i > $#commands;
        my $pid = fork();
        if (!defined $pid) {
            die "Failed to fork: $!";
        } elsif ($pid == 0) {
            print "Running command (even): $commands[$i]\n";
            system($commands[$i]) == 0 or die "Failed to execute command: $commands[$i]";
            exit(0);
        } else {
            push @pids, $pid;
        }
    }
    # Wait for even-indexed processes to finish
    for my $pid (@pids) {
        waitpid($pid, 0);
    }
    @pids = ();  # Reset for next batch
    # Step 2: Run odd-indexed commands in parallel
    for (my $i = 1; $i < ($num_parallel + 1) * 2; $i += 2) {
        last if $i > $#commands;
        my $pid = fork();
        if (!defined $pid) {
            die "Failed to fork: $!";
        } elsif ($pid == 0) {
            print "Running command (odd): $commands[$i]\n";
            system($commands[$i]) == 0 or die "Failed to execute command: $commands[$i]";
            exit(0);
        } else {
            push @pids, $pid;
        }
    }
    # Wait for odd-indexed processes to finish
    for my $pid (@pids) {
        waitpid($pid, 0);
    }
    # Step 3: Run remaining commands sequentially
    for (my $i = ($num_parallel + 1) * 2; $i <= $#commands; $i++) {
        print "Running command (sequential): $commands[$i]\n";
        system($commands[$i]) == 0 or die "Failed to execute command: $commands[$i]";
    }
}

sub parse_arguments {
    my %config = (
        params => {
            threads => 24,
            memory => 500,
            min_length => 500,
            num_selections => 1,
        },
        dirs => {
            output => "output",
            assembly => "assembly",
        },
        input => {
            control_reads_2 => "use_regex",
            selection_reads_2 => "use_regex",
        },
        programs => {
            assembly => "metaspades.py",
            clean_assembly => "assembly/clean_assembly.pl",
        }
    );
    GetOptions(
        "control-reads-1=s" => \$config{input}{control_reads_1},
        "control-reads-2=s" => \$config{input}{control_reads_2},
        "selection-reads-1=s" => \$config{input}{selection_reads_1},
        "selection-reads-2=s" => \$config{input}{selection_reads_2},
        "output-dir=s" => \$config{dirs}{output},
        "num-selections=s" => \$config{params}{num_selections},
        "threads=i" => \$config{params}{threads},
        "memory=i" => \$config{params}{memory},
        "min-length=i" => \$config{params}{min_length},
        "assembler=s" => \$config{programs}{assembly},
        "help|h" => \$config{help}
    ) or usage();

    if ($config{input}{control_reads_1}) {
        $config{input}{control_reads_2} = construct_paired_filename($config{input}{control_reads_1})
            if $config{input}{control_reads_2} eq "use_regex";
    }
    if ($config{input}{selection_reads_1}) {
        $config{input}{selection_reads_2} = construct_paired_filename($config{input}{selection_reads_1})
            if $config{input}{selection_reads_2} eq "use_regex";
    }
    if ($config{programs}{assembly} !~ /^(metaspades\.py)$/) {
        if ($config{programs}{assembly} !~ /^(rnaspades\.py)$/) {
            print "Supported assemblers are metaspades.py and rnaspades.py.\n";
            print "To use a different assembler please install the program\n";
            print "and modify the script to include the correct command.\n";
            exit(1);
        }
        $config{programs}{assembly} = "rnaspades.py";
    }

    $config{dirs}{assembly} = $config{dirs}{output} . "/assembly";
    make_dir($config{dirs}{output});
    make_dir($config{dirs}{assembly});
    
    usage() if $config{help};

    return \%config;
}

sub usage {
    print <<EOF;
Usage: $0 [options]

Required:
    --control-reads-1 FILE    Forward reads for control sample
    --selection-reads-1 FILE  Forward reads for selection sample
    --control-reads-2 FILE    Reverse reads for control sample
    --selection-reads-2 FILE  Reverse reads for selection sample
    --output-dir DIR          Output directory path

Optional:
    --assembler PROGRAM     Assembler program (default: metaspades.py)
    --num-selections STR    Number of selection files inputted into assembly
    --threads INT           Number of threads (default: 24)
    --memory INT            Memory in GB (default: 900)
    --min-length INT        Minimum contig length (default: 500)
    --help                  Show this help message

Example:
    $0 --control-reads-1 control_R1.fq.gz \\
       --control-reads-2 control_R2.fq.gz \\
       --selection-reads-1 selection_R1.fq.gz \\
       --selection-reads-2 selection_R2.fq.gz \\
       --output-dir results
EOF
    exit(1);
}

sub write_commands {
    my ($config) = @_;
    my @commands;
    
    my $generic_control = $config->{input}{control_reads_1};
    $generic_control =~ s/.1_val_1.fq.gz//;
    $generic_control =~ s/.*\///g;

    my @selection_files = split /,/, $config->{input}{selection_reads_1};
    my @generic_selections = map {
        my $name = $_;
        $name =~ s/.1_val_1.fq.gz//;
        $name =~ s/.*\///g;
        $name;
    } @selection_files;

    my $generic = $generic_control; $generic =~ s/control/experiment/;
    my $generic_prefix = $config->{dirs}{assembly}."/".$generic;

    my $control_prefix =  $config->{dirs}{assembly}."/".$generic_control;

    my $assembly_final_info = $generic_prefix."control_and_selected.fasta";
    my @control_commands = assemble_reads($config->{programs}{assembly}, 
                                        $config->{input}{control_reads_1}, $config->{input}{control_reads_2},
                                        [$control_prefix], $config->{params}{memory}, $config->{params}{threads},
                                        $config->{programs}{clean_assembly}, $config->{params}{min_length}, 1, "control");
    push @commands, @control_commands;    
    
    my $num_selections = $config->{params}{num_selections};
    my @selection_prefixes = map {$config->{dirs}{assembly} . "/" . $generic_selections[$_] } (0 .. $num_selections-1);
    my @selection_commands =  assemble_reads($config->{programs}{assembly}, 
                                            $config->{input}{selection_reads_1}, $config->{input}{selection_reads_2},
                                            \@selection_prefixes, $config->{params}{memory}, $config->{params}{threads},
                                            $config->{programs}{clean_assembly}, $config->{params}{min_length}, $num_selections, "selection");
    push @commands, @selection_commands;

    # Build list of selection fasta files
    my @selection_fasta_files = map { $selection_prefixes[$_].".fasta" } (0 .. $num_selections-1);
    my $selection_concat = join(" ", @selection_fasta_files);

    push @commands, "cat $selection_concat $control_prefix.fasta > $assembly_final_info";
    push @commands, "cd-hit-est -i ".$assembly_final_info.
        " -o ".$generic_prefix."control_and_selected.nd.fasta".
        " -c 0.95 -n 10 -M 0 -T ".$config->{params}{threads}.
        " -d 0 -r 1 -G 1 -g 1";
    push @commands, "samtools faidx ".$assembly_final_info;
    push @commands, "bwa index ".$assembly_final_info;
    push @commands, "metaquast.py " . $assembly_final_info.
        " -o ".$config->{dirs}{assembly}."/quast";

    return @commands
}

sub assemble_reads {
    my ($assembly_program, $read1, $read2, $prefix, $memory, $threads, $clean_assembly, $min_length, $num_files, $tag) = @_;
    my @commands;
    
    my @read1_files = split /,/, $read1;
    my @read2_files = split /,/, $read2;
    my @prefixes = @{$prefix};
    for (my $i = 0; $i < $num_files; $i++) {
        my $r1 = $read1_files[$i];
        my $r2 = ($read2_files[$i] // "");
        my $sub_prefix = $prefixes[$i];
        push @commands, $assembly_program." -1 ".$r1.
            " -2 ".$r2.
            " -o ".$sub_prefix."_assembly".
            " --only-assembler".
            " --memory ".$memory.
            " --threads ".$threads;
        push @commands, "perl ".$clean_assembly.
            " -file ".$sub_prefix."_assembly/contigs.fasta".
            " -tag $tag".
            " --out ".$sub_prefix.".fasta".
            " --min_length ".$min_length;
    }
    return @commands;
}
