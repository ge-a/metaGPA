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

    for my $cmd (@commands) {
        print "Running command: $cmd\n";
        system($cmd) == 0 or die "Failed to execute command: $cmd";
    }
}

sub parse_arguments {
    my %config = (
        params => {
            threads => 24,
            memory => 500,
            min_length => 500,
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
        "selection-num=i" => \$config{selection_num},
        "output-dir=s" => \$config{dirs}{output},
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
    --selection-num INT       The Nth selection being ran
    --output-dir DIR          Output directory path

Optional:
    --control-reads-2 FILE    Reverse reads for control sample (default: constructed from forward reads)
    --selection-reads-2 FILE  Reverse reads for selection sample (default: constructed from forward reads)
    --assembler PROGRAM       Assembler program (default: metaspades.py)
    --threads INT            Number of threads (default: 24)
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

    my $generic_selection = $config->{input}{selection_reads_1};
    $generic_selection =~ s/.1_val_1.fq.gz//;
    $generic_selection =~ s/.*\///g;

    my $generic = $generic_control; $generic =~ s/control/experiment/;
    my $generic_prefix = $config->{dirs}{assembly}."/".$generic."_".$config->{selection_num};

    my $control_prefix =  $config->{dirs}{assembly}."/".$generic_control;
    my $selection_prefix = $config->{dirs}{assembly}."/".$generic_selection."_".$config->{selection_num};

    my $assembly_final_info = $generic_prefix."control_and_selected.fasta";

    my $path = $control_prefix."_"."assembly";
    my $not_ran_control = (-e $path) ? 0 : 1;

    if ($not_ran_control == 1) {
        my @control_commands = assemble_reads($config->{programs}{assembly}, 
                                            $config->{input}{control_reads_1}, $config->{input}{control_reads_2},
                                            $control_prefix, $config->{params}{memory}, $config->{params}{threads},
                                            $config->{programs}{clean_assembly}, $config->{params}{min_length});
        push @commands, @control_commands;
    }
    
    my @selection_commands =  assemble_reads($config->{programs}{assembly}, 
                                            $config->{input}{selection_reads_1}, $config->{input}{selection_reads_2},
                                            $selection_prefix, $config->{params}{memory}, $config->{params}{threads},
                                            $config->{programs}{clean_assembly}, $config->{params}{min_length});
    push @commands, @selection_commands;

    push @commands, "cat ".$selection_prefix.".fasta ".$control_prefix.".fasta > ".
        $assembly_final_info;
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
    my ($assembly_program, $read1, $read2, $prefix, $memory, $threads, $clean_assembly, $min_length) = @_;
    my @commands;
    push @commands, $assembly_program." -1 ".$read1.
        " -2 ".$read2.
        " -o ".$prefix."_assembly".
        " --only-assembler".
        " --memory ".$memory.
        " --threads ".$threads;
    push @commands, "perl ".$clean_assembly.
        " -file ".$prefix."_assembly/contigs.fasta".
        " -tag control".
        " --out ".$prefix.".fasta".
        " --min_length ".$min_length;
    return @commands;
}