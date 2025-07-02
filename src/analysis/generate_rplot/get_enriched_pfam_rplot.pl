use strict;
use Cwd;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);
use Data::Dumper;
use Math::CDF qw(:all);
use utils qw(parse_pfam_file 
            calculate_enrichment_stats
            parse_bed);

exit main();

# Plan for this is to abstract out the code from get_annotated_enrichment into util function
# Then call it in that file and this file, and then focus on rewriting the write_enriched_pfam_results
# Package this script with the helpers in its own dir
# Then think about how we can add the gene delimeters
# Programatically see which genes of interest appear most often around our domain 

sub main {
    my $config = parse_arguments();

    my $hash_contig = $config->{unwanted} ? parse_bed($config->{unwanted}) : {};

    my ($result2, $result3, $pfam2desc) = parse_pfam_file(
        $config->{pfam_file},
        \$config{enrichment_file_path},
        $hash_contig,
        $config->{read_count_min},
        $config->{contig_min},
        $config->{pfam_cutoff},
        $config->{cutoff},
        $config->{direction},
    );

    my $result4 = calculate_enrichment_stats($result2, $result3);

    write_enriched_pfam_results(
        $result4,
        $pfam2desc,
        $config->{cutoff_evalue},
        $config->{out},
        $config->{pfam_file},
        $config->{contig_program},
        $config->{read_count_min},
        $config->{contig_min},
        $config->{cutoff},
        $config->{direction},
        $config->{pfam_cutoff},
        $config->{rscript}
    );

    return 0;
}

sub parse_arguments {
    my %config = (
            cutoff_evalue  => 0.05,
            cutoff         => 3,
            direction      => "enriched",
            read_count_min => 0,
            contig_min     => 0,
            pfam_cutoff    => 1e14,
            contig_program => "get_genomic_location.pl",
            rscript        => "plot.R",
            pfam           => "/mnt/home/ettwiller/laurence/database/pfam/Pfam-A.hmm",
        );

    GetOptions(
        "pfam=s"          => \$config{pfam_file},
        "enrich_txt=s"    => \$config{enrichment_file_path},
        "out=s"           => \$config{out},
        "cutoff=s"        => \$config{cutoff},
        "cutoff_evalue=s" => \$config{cutoff_evalue},
        "direction=s"     => \$config{direction},
        "read_count_min=s"=> \$config{read_count_min},
        "contig_min=s"    => \$config{contig_min},
        "unwanted=s"      => \$config{unwanted},
        "pfam_cutoff=s"   => \$config{pfam_cutoff},
        "help|h"          => \$config{help},
    ) or usage();

    usage() if $config{help} || !$config{out} || !$config{pfam_file};
    return \%config;
}

sub usage {
    print <<EOF;
Usage: $0 [options]

Required:
    --pfam FILE                Path to PFAM tab file
    --out FILE                 Output base name
    --enrich_txt FILE          enrichment_info txt file

Optional:
    --cutoff FLOAT             Enrichment cutoff value (default: 3)
    --cutoff_evalue FLOAT      Domain association e-value cutoff (default: 0.05)
    --direction STRING         Enrichment direction: enriched or depleted (default: enriched)
    --read_count_min INT       Minimum total reads in control+enriched (default: 0)
    --contig_min INT           Minimum contig length in bp (default: 0)
    --unwanted FILE            BED file of unwanted contigs (default: none)
    --pfam_cutoff FLOAT        PFAM E-value cutoff (default: 1e14)
    --help, -h                 Show this help message

Example:
    $0 --pfam pfam_hits.tab \\
       --out enriched_pfam \\
       --cutoff 3 \\
       --direction enriched \\
       --read_count_min 100 \\
       --contig_min 500

EOF
    exit(1);
}


sub write_enriched_pfam_results {
    my (
        $result4, $pfam2desc, $cutoff_evalue, $out, $pfam_file, $contig_program,
        $read_count_min, $contig_min, $cutoff, $direction, $pfam_cutoff, $rscript
    ) = @_;

    foreach my $e (sort { $a <=> $b } keys %$result4) {
        my $pfams = $result4->{$e};
        foreach my $pfam (keys %$pfams) {
            my $x = $pfams->{$pfam};
            my $description = $pfam2desc->{$pfam};
            if ($e < $cutoff_evalue) {
                my $file_out = $out . "_$pfam.txt";
                print STDERR "$e $pfam $x $description\n";
                my $command = "perl $contig_program --out $file_out --pfam_hit $pfam_file --pfam $pfam --read_count_min $read_count_min --contig_min $contig_min --flanks 1000 --cutoff $cutoff --direction $direction --pfam_cutoff $pfam_cutoff";
                print STDERR "$command\n";
                system($command);
                my $file_png = $out . "_$pfam.png";
                my $commandR = "Rscript $rscript $file_out $file_png";
                print STDERR "$commandR\n";
                system($commandR);
            }
        }
    }
}