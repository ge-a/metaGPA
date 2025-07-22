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
use utils qw(make_dir);

exit main();

sub main {
    my $config = parse_arguments();
    my @commands = write_commands($config);

    foreach my $command (@commands) {
        print "Running command: $command\n";
        system($command) == 0 or die "Failed to execute command: $command";
    }
    return 0;
}

sub parse_arguments {
    my %config = (
        dirs => {
            output => "output",
            mapping => "enrichment",
        },
        programs => {
            get_enrich => File::Spec->catfile($Bin, "multi_annotate_enrich.pl"),
            add_enrich => File::Spec->catfile($Bin, "multi_get_enrich_txt.pl"),
        },
        cutoff => 3,
    );
    GetOptions(
        "generic-control=s" => \$config{input}{control},
        "generic-selection=s" => \$config{input}{selection},
        "prefix=s" => \$config{prefix},
        "num-selections=s" => \$config{num_selections},
        "output-dir=s" => \$config{dirs}{output},
        "cutoff=f" => \$config{cutoff},
        "help|h" => \$config{help}
    ) or usage();
    
    $config{dirs}{enrichment} = $config{dirs}{output}."/enrichment";
    make_dir($config{dirs}{output});
    make_dir($config{dirs}{enrichment});

    usage() if $config{help};

    return \%config;
}

sub usage {
    print <<EOF;
Usage: $0 [options]

Required:
    --generic-control FILE    Generic name for control sample (e.g., PF4.1)
    --generic-selection FILE  Generic name for selection sample (e.g., PF1.1)
    --prefix STRING          Prefix for output files
    --output-dir DIR         Output directory path

Optional:
    --cutoff FLOAT           Cutoff value for enrichment, if -1 calculated from distribution (default: 3)
    --help                   Show this help message

Example:
    $0 --generic-control PF4.1 \\
       --generic-selection PF1.1 \\
       --prefix experiment \\
       --output-dir results
EOF
    exit(1);
}

sub write_commands {
    my ($config) = @_;
    my @commands;

    my $prefix = $config->{prefix};

    my $fai = $config->{dirs}{output}."/assembly/".$prefix."control_and_selected.nd.fasta.fai";
    my $assembly_final = $config->{dirs}{output}."/assembly/".$prefix."control_and_selected.nd.fasta";
    my $control_bam = $config->{dirs}{output}."/mapping/".$config->{input}{control}."mapped_to_".$prefix.".bam";
    my $control_dedup = $config->{dirs}{output}."/mapping/".$config->{input}{control}."mapped_to_".$prefix."duplicated_remove.bam";
    my $coverage_bed = $config->{dirs}{output}."/mapping/".$prefix."control_versus_selection_coverage.bed";
    my $hmmer2 = $config->{dirs}{output}."/annotation/".$prefix."control_and_selected_hmmer_format.tab";
    my $hmmer4 = $config->{dirs}{output}."/annotation/".$prefix."control_and_selected_hmmer_format_TIGRFAM.tab";
    my $pfam = $config->{dirs}{enrichment}."/".$prefix."control_and_selected_hmmer_format_PFAM.enriched";
    my $tigrfam = $config->{dirs}{enrichment}."/".$prefix."control_and_selected_hmmer_format_TIGRFAM.enriched";
    my $enrichment_info = $config->{dirs}{enrichment}."/".$prefix."enrichment_info.txt";
    
    my $num_selection = $config->{num_selections};
    my @generic_selections = split(/,/,$config->{input}{selection});
    my @selection_bams = map {
        $config->{dirs}{output}."/mapping/".$generic_selections[$_]."mapped_to_".$prefix.".bam"
    } (0 .. $num_selection-1);

    my @selection_dedups = map {
        $config->{dirs}{output}."/mapping/".$generic_selections[$_]."mapped_to_".$prefix."duplicated_remove.bam"
    } (0 .. $num_selection-1);

    # Format: control_bam control_dedup selection_bam_0 selection_dedup_0 selection_bam_1 selection_dedup_1 ...
    my $bam_list = "$control_bam $control_dedup";
    for (my $i = 0; $i < $num_selection; $i++) {
        $bam_list .= " $selection_bams[$i] $selection_dedups[$i]";
    }
    my $selection_list = join(",", @selection_bams);
    push @commands, "awk 'BEGIN{ OFS=\"\\t\"; }{ print \$1, 1, \$2, FILENAME, 100, \"+\"; }' ".
        $fai." | bedtools multicov -bed - -bams $bam_list > $coverage_bed";
    push @commands, "perl ".$config->{programs}{add_enrich}.
        " --fasta ".$assembly_final.
        " --enrichment ".$coverage_bed.
        " --control-bam ".$control_bam.
        " --selection-bam-str ".$selection_list.
        " --out ".$enrichment_info;
    push @commands, "perl ".$config->{programs}{get_enrich}.
        " --pfam ".$hmmer2.
        " --enrichment_txt ".$enrichment_info.
        " --out ".$pfam.
        " --cutoff ".$config->{cutoff};
    push @commands, "perl ".$config->{programs}{get_enrich}.
        " --pfam ".$hmmer4.
        " --enrichment_txt ".$enrichment_info.
        " --out ".$tigrfam.
        " --cutoff ".$config->{cutoff};
    return @commands;
}