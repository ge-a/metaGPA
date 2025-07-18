use strict;
use warnings;
use Cwd;
use Getopt::Long;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);
use Bio::SeqIO;
use utils qw(make_dir 
            construct_paired_filename);

exit main();

sub main {
    my $config = parse_arguments();
    my @commands = write_commands($config);

    my $num_selection = $config->{params}{num_selections}

    # Run commands 0, 1, 2 together (group 1)
    my @pids1;
    for my $i (0 .. 2) {
        my $pid = fork();
        if (!defined $pid) {
            die "Failed to fork: $!";
        } elsif ($pid == 0) {
            print "Running command: $commands[$i]\n";
            system($commands[$i]) == 0 or die "Failed to execute command: $commands[$i]";
            exit(0);
        } else {
            push @pids1, $pid;
        }
    }
    for my $pid (@pids1) {
        waitpid($pid, 0);
    }
    # Run commands 4 to 4+num_selection-1 concurrently (group 2)
    my @pids2;
    for my $i (4 .. 4 + $num_selection - 1) {
        my $pid = fork();
        if (!defined $pid) {
            die "Failed to fork: $!";
        } elsif ($pid == 0) {
            print "Running command: $commands[$i]\n";
            system($commands[$i]) == 0 or die "Failed to execute command: $commands[$i]";
            exit(0);
        } else {
            push @pids2, $pid;
        }
    }
    for my $pid (@pids2) {
        waitpid($pid, 0);
    }

    # Run commands 4+num_selection to 4+(3*num_selection)-1 in pairs of 2 concurrently (group 3)
    my $start = 4 + $num_selection;
    my $end = 4 + (3 * $num_selection) - 1;
    my @indices = ($start .. $end);

    while (@indices) {
        my @pair = splice(@indices, 0, 2);
        my @pids3;
        for my $i (@pair) {
            next if $i > $#commands;
            my $pid = fork();
            if (!defined $pid) {
                die "Failed to fork: $!";
            } elsif ($pid == 0) {
                print "Running command: $commands[$i]\n";
                system($commands[$i]) == 0 or die "Failed to execute command: $commands[$i]";
                exit(0);
            } else {
                push @pids3, $pid;
            }
        }
        for my $pid (@pids3) {
            waitpid($pid, 0);
        }
    }
    return 0;
}

sub parse_arguments {
    my %config = (
        params => {
            num_selections => 1,
        },
        dirs => {
            output => "output",
            mapping => "mapping"
        },
        input => {
            control_reads_2 => "use_regex",
            selection_reads_2 => "use_regex"
        },
        programs => {
            mapping => "mapping/run_dna_map.pl",
            picard => "/mnt/home/ettwiller/anaconda3/share/picard-2.27.1-0/picard.jar",
            mapping_func => "bwa"
        }
    );
    GetOptions(
        "control-reads-1=s"   => \$config{input}{control_reads_1},
        "control-reads-2=s"   => \$config{input}{control_reads_2},
        "selection-reads-1=s" => \$config{input}{selection_reads_1},
        "selection-reads-2=s" => \$config{input}{selection_reads_2},
        "prefix=s"            => \$config{params}{prefix},
        "num-selections"      => \$config{params}{num_selections}
        "output-dir=s"        => \$config{dirs}{output},
        "mapping-func=s"      => \$config{programs}{mapping_func},
        "help|h"              => \$config{help},
    ) or usage();

    if ($config{input}{control_reads_1}) {
        $config{input}{control_reads_2} = construct_paired_filename($config{input}{control_reads_1})
            if $config{input}{control_reads_2} eq "use_regex";
    }
    if ($config{input}{selection_reads_1}) {
        $config{input}{selection_reads_2} = construct_paired_filename($config{input}{selection_reads_1})
            if $config{input}{selection_reads_2} eq "use_regex";
    }

    $config{dirs}{mapping} = $config{dirs}{output}."/mapping";
    make_dir($config{dirs}{output});
    make_dir($config{dirs}{mapping});

    usage() if $config{help};

    return \%config;
}

sub usage {
    print <<EOF;
Usage: $0 [options]

Required:
    --control-reads-1 FILE      Forward reads for control sample
    --selection-reads-1 FILE    Forward reads for selection sample
    --output-dir DIR           Output directory path

Optional:
    --control-reads-2 FILE     Reverse reads for control sample (default: auto-detect)
    --selection-reads-2 FILE   Reverse reads for selection sample (default: auto-detect)
    --num-selections STR    Number of selection files inputted into assembly
    --mapping-func FUNC        Mapping function to use (default: bwa)
    --help                     Show this help message

Example:
    $0 --control-reads-1 control_1_val_1.fq.gz \\
       --selection-reads-1 selection_1_val_1.fq.gz \\
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

    my $generic = $config->{params}{prefix}

    my $assembly_final = $config->{dirs}{output}."/assembly/".$generic."control_and_selected.fasta";
    my $control_bam = $config->{dirs}{mapping}."/".$generic_control."mapped_to_".$generic.".bam";
    my $selection_bam = $config->{dirs}{mapping}."/".$generic_selection."mapped_to_".$generic.".bam";
    my $control_dedup = $config->{dirs}{mapping}."/".$generic_control."mapped_to_".$generic."duplicated_remove.bam";
    my $selection_dedup = $config->{dirs}{mapping}."/".$generic_selection."mapped_to_".$generic."duplicated_remove.bam";
    my $control_stats = $config->{dirs}{mapping}."/".$generic_control."mapped_to_".$generic."duplicated_remove.txt";
    my $selection_stats = $config->{dirs}{mapping}."/".$generic_selection."mapped_to_".$generic."duplicated_remove.txt";

    my $num_control = 1;
    my $num_selection = $config->{params}{num_selections}

    push @commands, mapping_commands(
        $config->{programs}{mapping},
        $config->{input}{control_reads_1},
        $config->{input}{control_reads_2},
        $assembly_final,
        $config->{programs}{mapping_func},
        $control_bam,
        $num_control
    );
    push @commands, "java -jar ".$config->{programs}{picard}.
        " MarkDuplicates REMOVE_DUPLICATES=true".
        " I=".$control_bam.
        " O=".$control_dedup.
        " M=".$control_stats;
    push @commands, "samtools index ".$control_dedup;

    push @commands, mapping_commands(
        $config->{programs}{mapping},
        $config->{input}{selection_reads_1},
        $config->{input}{selection_reads_2},
        $assembly_final,
        $config->{programs}{mapping_func},
        $selection_bam,
        $num_selection
    );
    for (my $i = 0; $i < $num_selection; $i++) {
        my $selection_bam_file = $selection_bam . ($num_selection > 1 ? "_$i.bam" : ".bam");
        my $selection_dedup_file = $selection_dedup . ($num_selection > 1 ? "_$i.bam" : ".bam");
        my $selection_stats_file = $selection_stats . ($num_selection > 1 ? "_$i.txt" : ".txt");
        push @commands, "java -jar ".$config->{programs}{picard}.
            " MarkDuplicates REMOVE_DUPLICATES=true".
            " I=".$selection_bam_file.
            " O=".$selection_dedup_file.
            " M=".$selection_stats_file;
        push @commands, "samtools index ".$selection_dedup_file;
    }
    return @commands;
}

sub mapping_commands {
    my ($mapping_program, $read1, $read2, $genome, $mapping_func, $out_prefix, $num_files) = @_;
    my @commands;

    my @read1_files = split /,/, $read1;
    my @read2_files = split /,/, $read2;

    for (my $i = 0; $i < $num_files; $i++) {
        my $r1 = $read1_files[$i] // "";
        my $r2 = $read2_files[$i] // "";
        my $sub_out = $out_prefix . ($num_files > 1 ? "_$i.bam" : ".bam");
        push @commands, "perl $mapping_program".
            " --fq1 $r1".
            " --fq2 $r2".
            " --genome $genome".
            " --mapping-func $mapping_func".
            " --out $sub_out";
    }
    return @commands;
}