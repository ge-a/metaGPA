use strict;
use warnings;
use Cwd;
use Getopt::Long;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);
use Bio::SeqIO;

exit main();

sub main {
    my $config = parse_arguments();
    my @commands = write_commands($config);

    my @pids;
    for my $i (0, 1) {
        my $pid = fork();
        if (!defined $pid) {
            die "Failed to fork: $!";
        } elsif ($pid == 0) {
            print "Running command: $commands[$i]\n";
            system($commands[$i]) == 0 or die "Failed to execute command: $commands[$i]";
            exit(0);
        } else {
            push @pids, $pid;
        }
    }
    for my $pid (@pids) {
        waitpid($pid, 0);
    }
    for my $i (2 .. $#commands) {
        print "Running command: $commands[$i]\n";
        system($commands[$i]) == 0 or die "Failed to execute command: $commands[$i]";
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
    system("mkdir -p ".$config{dirs}{output}) unless -d $config{dirs}{output};
    system("mkdir -p ".$config{dirs}{assembly}) unless -d $config{dirs}{assembly};
    
    usage() if $config{help};

    return \%config;
}

sub usage {
    print <<EOF;
Usage: $0 [options]

Required:
    --control-reads-1 FILE    Forward reads for control sample
    --selection-reads-1 FILE  Forward reads for selection sample
    --output-dir DIR         Output directory path

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

sub construct_paired_filename {
    my ($read1_file) = @_;
    my $read2_file = $read1_file;
    $read2_file =~ s/1_val_1.fq.gz/2_val_2.fq.gz/;
    return $read2_file;
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
    my $generic_prefix = $config->{dirs}{assembly}."/".$generic;

    my $control_prefix =  $config->{dirs}{assembly}."/".$generic_control;
    my $selection_prefix = $config->{dirs}{assembly}."/".$generic_selection;

    my $assembly_final_info = $generic_prefix."control_and_selected.fasta";
    my $translated_assembly_final = $config->{dirs}{assembly}."/".$generic."control_and_selected.faa";
    
    push @commands, $config->{programs}{assembly}." -1 ".$config->{input}{control_reads_1}.
        " -2 ".$config->{input}{control_reads_2}.
        " -o ".$control_prefix."_assembly".
        " --only-assembler".
        " --memory ".$config->{params}{memory}.
        " --threads ".$config->{params}{threads};
    push @commands, $config->{programs}{assembly}." -1 ".$config->{input}{selection_reads_1}.
        " -2 ".$config->{input}{selection_reads_2}.
        " -o ".$selection_prefix."_assembly".
        " --only-assembler".
        " --memory ".$config->{params}{memory}.
        " --threads ".$config->{params}{threads};
    push @commands, "perl ".$config->{programs}{clean_assembly}.
        " -file ".$control_prefix."_assembly/contigs.fasta".
        " -tag control".
        " --out ".$control_prefix.".fasta".
        " --min_length ".$config->{params}{min_length};
    push @commands, "perl ".$config->{programs}{clean_assembly}.
        " -file ".$selection_prefix."_assembly/contigs.fasta".
        " -tag selection".
        " --out ".$selection_prefix.".fasta".
        " --min_length ".$config->{params}{min_length};
    push @commands, "cat ".$selection_prefix.".fasta ".$control_prefix.".fasta > ".
        $generic_prefix."control_and_selected.fasta";
    push @commands, "cd-hit-est -i ".$generic_prefix."control_and_selected.fasta".
        " -o ".$generic_prefix."control_and_selected.nd.fasta".
        " -c 0.95 -n 10 -M 0 -T ".$config->{params}{threads}.
        " -d 0 -r 1 -G 1 -g 1";
    push @commands, "samtools faidx ".$generic_prefix."control_and_selected.fasta";
    push @commands, "bwa index ".$generic_prefix."control_and_selected.fasta";

    return @commands
}