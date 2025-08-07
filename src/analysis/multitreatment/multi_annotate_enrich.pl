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
use lib "$Bin/../..";
use utils qw(parse_pfam_file 
            calculate_enrichment_stats
            parse_bed);

exit main();

sub main {
    my ($pfam, $enrichment_file_path, $num_selection, $out, $cutoff, $direction, $read_count_min, $contig_min, $unwanted, $pfam_cutoff, $TOTAL) = parse_args();

    my $hash_contig = $unwanted ? parse_bed($unwanted) : undef;

    if ($cutoff == -1) {
        $cutoff = get_cutoff_from_distribution($enrichment_file_path, 1);
    }
    for (my $i = 1; $i <= $num_selection; $i++) {
        my $enrich_file_i = parse_multi_enrich_txt($enrichment_file_path, $i);
        my ($result2, $result3, $pfam2desc) = parse_pfam_file($pfam, $enrich_file_i, $hash_contig, $read_count_min, $contig_min, $pfam_cutoff, $cutoff, $direction, 1);
        my $result4 = calculate_enrichment_stats($result2, $result3);

        my $out_i = $out;
        $out_i =~ s/(\.\w+)?$/_$i$1/;

        open(my $out_fh, ">", $out_i) or die "Can't open $out_i\n";
        write_enrichment_output($result4, $pfam2desc, $TOTAL, $out_fh);
        close $out_fh;
    }
    return 0;
}

sub parse_multi_enrich_txt {
    my ($enrichment_file_path, $selection_num) = @_;
    open(my $fh, '<', $enrichment_file_path) or die "Cannot open file: $enrichment_file_path\n";

    my @result;
    my $header = <$fh>;  # skip header line

    while (my $line = <$fh>) {
        chomp $line;
        my @fields = split /\t/, $line;

        # Extract fields
        my $id = $fields[0];
        my $control_count = $fields[1];

        # Determine column indices for the desired selection number
        my $sel_count_idx = ($selection_num * 2);
        my $sel_ratio_idx = $sel_count_idx + 1;

        my $selection_count = $fields[$sel_count_idx];
        my $selection_ratio = $fields[$sel_ratio_idx];

        # Store as a hash (or any structure you want)
        push @result, {
            id => $id,
            control_count => $control_count,
            selection_count => $selection_count,
            selection_ratio => $selection_ratio,
        };
    }
    close $fh;
    return \@result;  # return a reference to the array
}

sub write_enrichment_output {
    my ($pvalue_to_pfam_stats, $pfam_to_description, $min_total_count, $out_fh) = @_;

    # Print header (optional, comment out if not needed)
    print $out_fh join("\t", "p_value", "pfam_id", "enriched", "total", 
        "local_ratio", "global_enrichment_prob", "description"), "\n";
    foreach my $p_value (sort { $a <=> $b } keys %$pvalue_to_pfam_stats) {
        my $pfam_stats = $pvalue_to_pfam_stats->{$p_value};
        foreach my $pfam_id (keys %$pfam_stats) {
            my $stats_ref = $pfam_stats->{$pfam_id};
            my $stats_string = $stats_ref->{stats};  # e.g., "6_10_0.6000"
            my $global_enrichment_prob = $stats_ref->{global_enrichment_probability};

            my ($enriched_count, $total_count, $local_ratio) = split /_/, $stats_string;
            # Apply filter: only write if total count > threshold
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
    my ($pfam, $enrichment_file_path, $num_selection, $out, $cutoff, $direction, $read_count_min, $contig_min, $unwanted, $pfam_cutoff, $TOTAL);
    $cutoff = 3;
    $direction = "enriched";
    $read_count_min = 0;
    $contig_min = 0;
    $pfam_cutoff = 100000000000000;
    #from : http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.gz
    $pfam = "/mnt/home/ettwiller/laurence/database/pfam/Pfam-A.hmm";
    $TOTAL = 0;

    GetOptions(
        "pfam=s"          => \$pfam, # hmmer format file (tab delimited).
        "enrichment-txt=s"=> \$enrichment_file_path, # enrichment file path.
        "num-selection=s" => \$num_selection,
        "out=s"           => \$out, # output file.
        "cutoff=s"        => \$cutoff, # cutoff for the enrichment
        "direction=s"     => \$direction, # depeted or enriched
        "read_count_min=s"=> \$read_count_min, # minimum (total) read number
        "contig_min=s"    => \$contig_min, # minimum size contig (bp)
        "unwanted=s"      => \$unwanted,
        "pfam_cutoff=s"   => \$pfam_cutoff,
        "pfam_number=s"   => \$TOTAL
    ) or die $error_sentence;

    die $error_sentence unless $out && $pfam;
    return ($pfam, $enrichment_file_path, $num_selection, $out, $cutoff, $direction, $read_count_min, $contig_min, $unwanted, $pfam_cutoff, $TOTAL);
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