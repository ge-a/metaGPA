use strict;
use warnings;
use Cwd;
use FindBin qw($Bin); 
use File::Spec;
use Getopt::Long;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);
use Bio::SeqIO;
use lib "$Bin/..";
use utils qw(make_dir 
            construct_paired_filename);

exit main();

sub main {
    my $config = parse_arguments();
    my @commands = write_commands($config);

    my $pid1 = fork();
    if (!defined $pid1) {
        die "Failed to fork: $!";
    } elsif ($pid1 == 0) {
        foreach my $i (0, 2, 4) {
            print "Running command: $commands[$i]\n";
            system($commands[$i]) == 0 or die "Failed to execute command: $commands[$i]";
        }
        exit(0);
    }
    my $pid2 = fork();
    if (!defined $pid2) {
        die "Failed to fork: $!";
    } elsif ($pid2 == 0) {
        foreach my $i (1, 3, 5) {
            print "Running command: $commands[$i]\n";
            system($commands[$i]) == 0 or die "Failed to execute command: $commands[$i]";
        }
        exit(0);
    }

    # Wait for both child processes to finish
    waitpid($pid1, 0);
    waitpid($pid2, 0);

    return 0;
}

sub parse_arguments {
    my %config = (
        dirs => {
            output => "output",
            mapping => "mapping"
        },
        input => {
            control_reads_2 => "use_regex",
            selection_reads_2 => "use_regex"
        },
        programs => {
            mapping => File::Spec->catfile($Bin, "run_dna_map.pl"),
            picard => "/mnt/home/ettwiller/anaconda3/share/picard-2.27.1-0/picard.jar",
            mapping_func => "bwa"
        }
    );
    GetOptions(
        "control-reads-1=s"   => \$config{input}{control_reads_1},
        "control-reads-2=s"   => \$config{input}{control_reads_2},
        "selection-reads-1=s" => \$config{input}{selection_reads_1},
        "selection-reads-2=s" => \$config{input}{selection_reads_2},
        "prefix=s"            => \$config{prefix},
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

    my $generic = $config->{prefix};

    my $assembly_final = $config->{dirs}{output}."/assembly/".$generic."control_and_selected.nd.fasta";
    my $control_bam = $config->{dirs}{mapping}."/".$generic_control."mapped_to_".$generic.".bam";
    my $selection_bam = $config->{dirs}{mapping}."/".$generic_selection."mapped_to_".$generic.".bam";
    my $control_dedup = $config->{dirs}{mapping}."/".$generic_control."mapped_to_".$generic."duplicated_remove.bam";
    my $selection_dedup = $config->{dirs}{mapping}."/".$generic_selection."mapped_to_".$generic."duplicated_remove.bam";
    my $control_stats = $config->{dirs}{mapping}."/".$generic_control."mapped_to_".$generic."duplicated_remove.txt";
    my $selection_stats = $config->{dirs}{mapping}."/".$generic_selection."mapped_to_".$generic."duplicated_remove.txt";

    push @commands, "perl ".$config->{programs}{mapping}.
        " --fq1 ".$config->{input}{control_reads_1}.
        " --fq2 ".$config->{input}{control_reads_2}.
        " --genome ".$assembly_final. 
        " --mapping-func ".$config->{programs}{mapping_func}.
        " --out ".$control_bam;
    push @commands, "perl ".$config->{programs}{mapping}.
        " --fq1 ".$config->{input}{selection_reads_1}.
        " --fq2 ".$config->{input}{selection_reads_2}.
        " --genome ".$assembly_final.
        " --mapping-func ".$config->{programs}{mapping_func}.
        " --out ".$selection_bam; 
    push @commands, "java -jar ".$config->{programs}{picard}.
        " MarkDuplicates REMOVE_DUPLICATES=true".
        " I=".$control_bam.
        " O=".$control_dedup.
        " M=".$control_stats;
    push @commands, "java -jar ".$config->{programs}{picard}.
        " MarkDuplicates REMOVE_DUPLICATES=true".
        " I=".$selection_bam.
        " O=".$selection_dedup.
        " M=".$selection_stats;
    push @commands, "samtools index ".$control_dedup;
    push @commands, "samtools index ".$selection_dedup;

    return @commands;
}