use strict;
use warnings;
use Cwd;
use FindBin qw($Bin); 
use File::Spec;
use Getopt::Long qw(GetOptions);
use Statistics::Descriptive;
use File::Temp qw(tempfile);
use Data::Dumper;
use Math::CDF qw(:all);
use lib "$Bin/..";
use utils qw(parse_pfam_file 
            calculate_enrichment_stats
            parse_bed);

exit main();

sub main {
    my $config = parse_args();

    my $hash_contig = $config->{unwanted} ? parse_bed($config->{unwanted}) : undef;

    if ($config->{cutoff} == -1) {
        $config->{cutoff} = get_cutoff_from_distribution($config->{enrichment_txt}, 1);
    }

    my ($result2, $result3, $pfam2desc) = parse_pfam_file(
        $config->{pfam},
        $config->{enrichment_txt},
        $hash_contig,
        $config->{read_count_min},
        $config->{contig_min},
        $config->{pfam_cutoff},
        $config->{cutoff},
        $config->{direction}
    );
    my $result4 = calculate_enrichment_stats($result2, $result3);

    open(my $out_fh, ">", $config->{out}) or die "Can't open $config->{out}\n";
    write_enrichment_output($result4, $pfam2desc, $config->{pfam_number}, $out_fh);
    close $out_fh;

    return 0;
}

sub write_enrichment_output {
    my ($pvalue_to_pfam_stats, $pfam_to_description, $min_total_count, $out_fh) = @_;

    print $out_fh join("\t", "p_value", "pfam_id", "enriched", "total", 
        "local_ratio", "global_enrichment_prob", "description"), "\n";
    foreach my $p_value (sort { $a <=> $b } keys %$pvalue_to_pfam_stats) {
        my $pfam_stats = $pvalue_to_pfam_stats->{$p_value};
        foreach my $pfam_id (keys %$pfam_stats) {
            my $stats_ref = $pfam_stats->{$pfam_id};
            my $stats_string = $stats_ref->{stats};  # e.g., "6_10_0.6000"
            my $global_enrichment_prob = $stats_ref->{global_enrichment_probability};

            my ($enriched_count, $total_count, $local_ratio) = split /_/, $stats_string;
            if ($total_count > $min_total_count) {
                my $description = $pfam_to_description->{$pfam_id} // "NA";
                print $out_fh join("\t", $p_value, $pfam_id, $enriched_count,
                    $total_count, $local_ratio, $global_enrichment_prob, $description), "\n";
            }
        }
    }
}

sub parse_args {
    my $error_sentence = "USAGE : perl $0 --pfam pfam_tab --out output_file OPTIONAL : --cutoff 10 (default 3) --direction depleted (default enriched) --read_count_min 100 (total number of reads in control+enriched, default 0) --contig_min 500 (length of the contig in bp, default 0) --unwanted contif.bed (defaut NONE) --pfam_cutoff 0.0000001 (this is the pfam Evalue. Anything matching of below this cutoff is used, DEFAULT no cutoff) --pfam_number 10 (at least 10 instances of a particular pfam, DEFAULT=0)";
    my %config = (
        cutoff        => 3,
        direction     => "enriched",
        read_count_min=> 0,
        contig_min    => 0,
        pfam_cutoff   => 100000000000000,
        pfam          => "/mnt/home/ettwiller/laurence/database/pfam/Pfam-A.hmm",
        pfam_number   => 0,
    );

    GetOptions(
        "pfam=s"          => \$config{pfam},
        "enrichment_txt=s"=> \$config{enrichment_txt},
        "out=s"           => \$config{out},
        "cutoff=s"        => \$config{cutoff},
        "direction=s"     => \$config{direction},
        "read_count_min=s"=> \$config{read_count_min},
        "contig_min=s"    => \$config{contig_min},
        "unwanted=s"      => \$config{unwanted},
        "pfam_cutoff=s"   => \$config{pfam_cutoff},
        "pfam_number=s"   => \$config{pfam_number},
        "help|h"          => \$config{help},
    ) or die $error_sentence;

    die $error_sentence unless $config{out} && $config{pfam};
    return \%config;
}

sub get_cutoff_from_distribution {
    my ($filename, $tails) = @_;
    my @ratios;

    open(my $fh, "<", $filename) or die "Can't open $filename: $!";
    while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ /^\s*$/; # skip empty lines
        my @fields = split /\t/, $line;
        push @ratios, $fields[-1];
    }
    close $fh;

    return -1 if @ratios < 5;

    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@ratios);
    my $skew = $stat->skewness();
    my $mean = $stat->mean();
    my $stddev = $stat->standard_deviation();

    my ($cutoff, $lower_cutoff, $upper_cutoff);

    if ($tails == 2) {
        my $N = (abs($skew) < 0.5) ? 2 : 3;
        $lower_cutoff = $mean - $N * $stddev;
        $upper_cutoff = $mean + $N * $stddev;
        return ($lower_cutoff, $upper_cutoff);
    } else {
        print(abs($skew)."\n");
        print(" $mean $stddev\n");
        if (abs($skew) < 0.5) {
            $cutoff = $mean + 2 * $stddev;
        } elsif ($skew > 0) {
            $cutoff = $mean + 3 * $stddev;
        } else {
            $cutoff = $mean - 3 * $stddev;
        }
        return $cutoff;
    }
}