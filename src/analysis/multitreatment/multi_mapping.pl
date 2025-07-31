use strict;
use warnings;
use Cwd;
use FindBin qw($Bin); 
use File::Spec;
use Getopt::Long;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);
use Bio::SeqIO;
use lib "$Bin/../..";
use utils qw(make_dir 
            construct_paired_filename);

exit main();


sub main {
    my $config = parse_arguments();
    my @commands = write_commands($config);

    my $total_commands = scalar @commands;
    my @pids;

    # Start every 3rd command in a separate process
    for (my $i = 0; $i < $total_commands; $i += 3) {
        my $pid = fork();
        if (!defined $pid) {
            die "Failed to fork: $!";
        } elsif ($pid == 0) {
            # Child process: run up to 3 commands sequentially
            for my $j (0 .. 2) {
                my $cmd_index = $i + $j;
                last if $cmd_index >= $total_commands;
                print "Running command [$cmd_index]: $commands[$cmd_index]\n";
                system($commands[$cmd_index]) == 0 or die "Failed to execute command [$cmd_index]: $commands[$cmd_index]";
            }
            exit(0);
        } else {
            # Parent process: track child PID
            push @pids, $pid;
        }
    }
    # Wait for all child processes to finish
    for my $pid (@pids) {
        waitpid($pid, 0);
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
            mapping => File::Spec->catfile($Bin, "..", "..", "mapping", "run_dna_map.pl"),
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
        "num-selections=s"     => \$config{params}{num_selections},
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

    my $num_selections = $config->{params}{num_selections};

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

    my $generic = $config->{params}{prefix};
    my $assembly_final = $config->{dirs}{output}."/assembly/".$generic."control_and_selected.nd.fasta";

    my $mapping_dir = $config->{dirs}{mapping};
    my $control_bam = $mapping_dir."/".$generic_control."mapped_to_".$generic.".bam";
    my $control_dedup = $mapping_dir."/".$generic_control."mapped_to_".$generic."duplicated_remove.bam";
    my $control_stats = $mapping_dir."/".$generic_control."mapped_to_".$generic."duplicated_remove.txt";

    my @selection_bams = map {
        $mapping_dir."/".$generic_selections[$_]."mapped_to_".$generic.".bam"
    } (0 .. $num_selections - 1);
    my @selection_dedups = map {
        $mapping_dir."/".$generic_selections[$_]."mapped_to_".$generic."duplicated_remove.bam"
    } (0 .. $num_selections - 1);
    my @selection_stats = map {
        $mapping_dir."/".$generic_selections[$_]."mapped_to_".$generic."duplicated_remove.txt"
    } (0 .. $num_selections - 1);

    push @commands, "perl ".$config->{programs}{mapping}.
            " --fq1 ".$config->{input}{control_reads_1}.
            " --fq2 ".$config->{input}{control_reads_2}.
            " --genome $assembly_final".
            " --mapping-func ".$config->{programs}{mapping_func}.
            " --out $control_bam";
    push @commands, "java -jar ".$config->{programs}{picard}.
        " MarkDuplicates REMOVE_DUPLICATES=true".
        " I=".$control_bam.
        " O=".$control_dedup.
        " M=".$control_stats;
    push @commands, "samtools index ".$control_dedup;

    my @read1_files = split /,/, $config->{input}{selection_reads_1};
    my @read2_files = split /,/, $config->{input}{selection_reads_2};

    for (my $i = 0; $i < $num_selections; $i++) {
        my $r1 = $read1_files[$i] // "";
        my $r2 = ($read2_files[$i] // "");
        my $selection_bam_file = $selection_bams[$i];
        my $selection_dedup_file = $selection_dedups[$i];
        my $selection_stats_file = $selection_stats[$i];
        push @commands, "perl ".$config->{programs}{mapping}.
            " --fq1 $r1".
            " --fq2 $r2".
            " --genome $assembly_final".
            " --mapping-func ".$config->{programs}{mapping_func}.
            " --out $selection_bam_file";
        push @commands, "java -jar ".$config->{programs}{picard}.
            " MarkDuplicates REMOVE_DUPLICATES=true".
            " I=".$selection_bam_file.
            " O=".$selection_dedup_file.
            " M=".$selection_stats_file;
        push @commands, "samtools index ".$selection_dedup_file;
    }
    return @commands;
}