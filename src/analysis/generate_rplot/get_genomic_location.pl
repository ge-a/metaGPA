#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);
use Bio::SeqIO;
use utils qw(parse_bed);

exit main();

sub main {
    my $config = parse_arguments();

    my $hash_contig = $config->{unwanted} ? parse_bed($config->{unwanted}) : {};

    my ($result, $all) = parse_pfam_hit(
        $config->{pfam_hit}, $config->{PFAM}, $config->{pfam_cutoff},
        $config->{read_count_min}, $config->{contig_min}, $hash_contig,
        $config->{cutoff}
    );

    my ($occurence, $final) = process_occurrences(
        $result, $all, $config->{flanks}, $config->{cutoff});

    assign_colors($occurence);

    write_output($final, $occurence, $config->{out});

    return 0;
}

sub parse_arguments {
    my %config = (
        read_count_min => 0,
        contig_min     => 0,
        flanks         => 500,
        direction      => "SELECTED",
        cutoff         => 3,
        pfam_cutoff    => 1e14,
    );

    GetOptions(
        "pfam_hit=s"      => \$config{pfam_hit},
        "pfam=s"          => \$config{PFAM},
        "out=s"           => \$config{out},
        "cutoff=s"        => \$config{cutoff},
        "read_count_min=s"=> \$config{read_count_min},
        "contig_min=s"    => \$config{contig_min},
        "unwanted=s"      => \$config{unwanted},
        "flanks=s"        => \$config{flanks},
        "direction=s"     => \$config{direction},
        "pfam_cutoff=s"   => \$config{pfam_cutoff},
        "help|h"          => \$config{help},
    ) or usage();

    usage() if $config{help} || !$config{pfam_hit} || !$config{PFAM} || !$config{out};
    return \%config;
}

sub usage {
    print <<EOF;
Usage: $0 [options]

Required:
    --pfam_hit FILE            PFAM hit tab file
    --pfam STRING              PFAM ID (e.g., PF000234)
    --out FILE                 Output file

Optional:
    --cutoff FLOAT             Enrichment cutoff value (default: 3)
    --read_count_min INT       Minimum total reads in control+enriched (default: 0)
    --contig_min INT           Minimum contig length in bp (default: 0)
    --flanks INT               Flanking region size (default: 500)
    --direction STRING         Enrichment direction: SELECTED or DEPLETED (default: SELECTED)
    --unwanted FILE            BED file of unwanted contigs (default: none)
    --pfam_cutoff FLOAT        PFAM E-value cutoff (default: 1e14)
    --help, -h                 Show this help message

Example:
    $0 --pfam_hit pfam_hits.tab \\
       --pfam PF000234 \\
       --out genomic_location.txt \\
       --cutoff 3 \\
       --read_count_min 100 \\
       --contig_min 500 \\
       --flanks 1000

EOF
    exit(1);
}

sub parse_pfam_hit {
    my ($pfam_hit, $PFAM, $pfam_cutoff, $read_count_min, $contig_min, $hash_contig, $cutoff) = @_;
    my (%result, %all);

    open(my $file_fh, "<", $pfam_hit) or die "can't open $pfam_hit\n";
    while (my $line = <$file_fh>) {
        chomp $line;
        next if $line =~ /^\#/;
        my @tmp = split /\s+/, $line;
        my $pfam_Evalue = $tmp[6];
        next unless $pfam_Evalue < $pfam_cutoff;

        my $contig = $tmp[0];
        $contig =~ /ratio_(.*)\-/;
        my $enrichment = $1;
        $contig =~ /length_(\S+)_/;
        my $contig_length = $1;
        $contig =~ /(\S+)_ENRICHMENT_control_(\S+)_selected_(\S+)_ratio_(\S+)/;
        my $contig_name_to_check = $1;
        my $control_reads = $2;
        my $enriched_reads = $3;
        my $read_count = $control_reads + $enriched_reads;

        if ($$hash_contig{$contig_name_to_check}) {
            print STDERR "$contig_name_to_check is DNA\n";
            next;
        }
        next unless $read_count > $read_count_min && $contig_length > $contig_min;

        my $category = ($enrichment > $cutoff) ? "SELECTED" : "DEPLETED";
        $contig =~ /length_(\S+)_/;
        my $contig_length2 = $1;
        $contig =~ /(\S+)-(\S+)/;
        my $contig_DNA = $1;
        my $orientation = "forward";
        my $frame = $2; $frame =~ s/\d+//; $orientation = "reverse" if $frame eq "R";
        my $start = $tmp[19];
        my $end = $tmp[20];
        my $pfam = $tmp[4];
        my $name = $tmp[3];

        if ($pfam eq $PFAM) {
            $result{$contig_DNA}{$start}++;
        }
        $all{$contig_DNA}{$start} = {
            category    => $category,
            pfam        => $name,
            end         => $end,
            orientation => $orientation,
        };
    }
    close $file_fh;
    return (\%result, \%all);
}

sub process_occurrences {
    my ($result, $all, $flanks, $cutoff) = @_;
    my %occurence;
    my @final;

    foreach my $contig (keys %$all) {
        my $hash1 = $all->{$contig};
        my $old = -1000000;
        if ($result->{$contig}) {
            my $hash2 = $result->{$contig};
            my $i = 0;
            foreach my $start2 (sort { $a <=> $b } keys %$hash2) {
                if (abs($start2 - $old) > $flanks) {
                    $i++;
                    $old = $start2;
                    my $new_contig = $contig . "-" . $i;
                    foreach my $start1 (keys %$hash1) {
                        if (abs($start1 - $start2) < $flanks) {
                            my $pfam = $hash1->{$start1}{"pfam"};
                            my $end1 = $hash1->{$start1}{"end"};
                            my $category = $hash1->{$start1}{"category"};
                            my $diff = $end1 - $start1;
                            my $start = $start1 - ($start2 - $flanks);
                            my $end = $start + $diff;
                            my $orientation = $hash1->{$start1}{"orientation"};
                            my $sens = ($orientation eq "reverse") ? 1 : 0;
                            $occurence{$pfam}{"selected"}++ if $category eq "SELECTED";
                            $occurence{$pfam}{"count"}++;
                            $occurence{$pfam}{"category"}++;
                            my $line = "$new_contig\t$pfam\t$start\t$end\t$orientation\t$sens\t$category";
                            push @final, $line;
                        }
                    }
                }
            }
        }
    }
    return (\%occurence, \@final);
}

sub assign_colors {
    my ($occurence) = @_;
    my @color = (
        "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#FFFF33",
        "#8DD3C7","#FFFFCC","#1B9E77","#BEAED4","#386CB0","#CCEBC5","#A65628","#FED9A6","#F4CAE4",
        "#7570B3","#FDDAEC","#E5C494","#FFD92F","#B15928"
    );
    my $i = 0;
    my $length_color = @color;
    foreach my $pfam (keys %$occurence) {
        my $v = $occurence->{$pfam}{"selected"};
        if ($v && $v > 0) {
            $occurence->{$pfam}{"color"} = $color[$i];
            $i = 0 if $i == $length_color - 1;
            $i++;
        } else {
            $occurence->{$pfam}{"color"} = "NA";
        }
    }
}

sub write_output {
    my ($final, $occurence, $out) = @_;
    open(my $out_fh, ">", $out) or die "can't open $out\n";
    foreach my $line (@$final) {
        my @tmp = split /\t/, $line;
        my $color = $occurence->{$tmp[1]}{"color"};
        print $out_fh "$line\t\"$color\"\n";
    }
    close $out_fh;
}
