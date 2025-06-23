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
            get_enrich => "enrichment/get_enriched.pl",
            add_enrich => "enrichment/write_enrichment.pl",
        },
    );
    GetOptions(
        "generic-control=s" => \$config{input}{control},
        "generic-selection=s" => \$config{input}{selection},
        "prefix=s" => \$config{prefix},
        "output-dir=s" => \$config{dirs}{output},
        "help|h" => \$config{help}
    ) or usage();
    $config{dirs}{enrichment} = $config{dirs}{output}."/enrichment";
    system("mkdir -p ".$config{dirs}{output}) unless -d $config{dirs}{output};
    system("mkdir -p ".$config{dirs}{enrichment}) unless -d $config{dirs}{enrichment};

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

    my $generic = $config->{input}{control}; $generic =~ s/control/experiment/;
    my $prefix = $config->{prefix};

    my $fai = $config->{dirs}{output}."/assembly/".$prefix."control_and_selected.fasta.fai";
    my $assembly_final = $config->{dirs}{output}."/assembly/".$generic."control_and_selected.fasta";
    my $control_bam = $config->{dirs}{output}."/mapping/".$config->{input}{control}."mapped_to_".$prefix.".bam";
    my $selection_bam = $config->{dirs}{output}."/mapping/".$config->{input}{selection}."mapped_to_".$prefix.".bam";
    my $control_dedup = $config->{dirs}{output}."/mapping/".$config->{input}{control}."mapped_to_".$generic."duplicated_remove.bam";
    my $selection_dedup = $config->{dirs}{output}."/mapping/".$config->{input}{control}."mapped_to_".$generic."duplicated_remove.bam";
    my $coverage_bed = $config->{dirs}{output}."/mapping/".$generic."control_versus_selection_coverage.bed";
    my $hmmer2 = $config->{dirs}{output}."/annotation/".$prefix."control_and_selected_hmmer_format.tab";
    my $hmmer4 = $config->{dirs}{output}."/annotation/".$prefix."control_and_selected_hmmer_format_TIGRFAM.tab";
    my $pfam = $config->{dirs}{enrichment}."/".$generic."control_and_selected_hmmer_format_PFAM.enriched";
    my $tigrfam = $config->{dirs}{enrichment}."/".$generic."control_and_selected_hmmer_format_TIGRFAM.enriched";

    push @commands, "awk 'BEGIN{ OFS=\"\\t\"; }{ print \$1, 1, \$2, FILENAME, 100, \"+\"; }' ".
        $fai." | bedtools multicov -bed - -bams ".
        $control_bam." ".
        $selection_bam." ".
        $control_dedup." ".
        $selection_dedup." > ".
        $coverage_bed;
    # This file needs to be rewritten to not append to end of file but create new file
    push @commands, "perl ".$config->{programs}{add_enrich}.
        " --fasta ".$assembly_final.
        " --enrichment ".$coverage_bed.
        " --out ".$config->{dirs}{enrichment}."/".$generic."enrichment_info.txt";
    push @commands, "perl ".$config->{programs}{get_enrich}.
        " --pfam ".$hmmer2.
        " --out ".$pfam.
        " --cutoff 3";
    push @commands, "perl ".$config->{programs}{get_enrich}.
        " --pfam ".$hmmer4.
        " --out ".$tigrfam.
        " --cutoff 3";
    return @commands;
}