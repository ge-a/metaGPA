#!/usr/bin/perl
use strict;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);
use Bio::SeqIO;

my $error_sentence = "USAGE : perl $0 --fasta fastafile.faa --pfam_hit pfamhit.tab --pfam PF000234 --read_count_min 100 --contig_min 500 --cutoff 10 --unwanted contig.bed\n";

exit main();

sub main {
    my $config = parse_args();
    residue_analysis($config);
    return 0;
}

sub usage {
    print <<EOF;
Usage: $0 [options]

Required:
    --fasta FILE        Combined assembly FASTA file (AA) with enrichment scores
    --pfam_hit FILE     PFAM hit tab file with position of pfam
    --pfam STRING       PFAM family ID (e.g., PF000234) desired

Optional:
    --cutoff FLOAT      Enrichment cutoff (default: 3)
    --read_count_min INT    Minimum total read number (default: 0)
    --contig_min INT        Minimum contig length (bp, default: 0)
    --unwanted FILE     BED file of unwanted contigs (default: NONE)

Example:
    $0 --fasta combined.faa \\
       --pfam_hit pfamhit.tab \\
       --pfam PF000234 \\
       --cutoff 10 \\
       --read_count_min 100 \\
       --contig_min 500 \\
       --unwanted contig.bed

EOF
    exit(1);
}

sub parse_args {
    my %config = (
        cutoff => 3,
        read_count_min => 0,
        contig_min => 0,
        print_only => "ALL",
    );
    GetOptions(
        "fasta=s"         => \$config{fasta_aa},
        "pfam_hit=s"      => \$config{pfam_hit},
        "pfam=s"          => \$config{pfam},
        "cutoff=f"        => \$config{cutoff},
        "read_count_min=i"=> \$config{read_count_min},
        "contig_min=i"    => \$config{contig_min},
        "unwanted=s"      => \$config{unwanted},
        "help|h"          => \$config{help},
    ) or usage();

    usage() if $config{help} || !$config{fasta_aa} || !$config{pfam_hit} || !$config{pfam};
    return \%config;
}

sub residue_analysis {
    my ($config) = @_;

    my $hash_contig;
    if (defined $config->{unwanted}) {
        $hash_contig = parse_bed($config->{unwanted});
    }

    my %id2AA;
    my $seq_in_aa = Bio::SeqIO->new(-file => $config->{fasta_aa}, -format => "fasta");
    while (my $faa = $seq_in_aa->next_seq()) {
        my $id = $faa->id;
        $id2AA{$id} = $faa;
    }

    my ($unselected, $selected) = (0, 0);

    open(my $hit_fh, "<", $config->{pfam_hit}) or die "can't open $config->{pfam_hit}\n";
    while (my $line = <$hit_fh>) {
        chomp $line;
        $line =~ s/ +/&/g;
        my @tmp = split /\&/, $line;
        my $domain = $tmp[4];

        if ($domain eq $config->{pfam}) {
            my $contig = $tmp[0];
            $contig =~ /ratio_(.*)\-(.)(.)/;
            my $contig_DNA = $contig;
            my $enrichment = $1;
            my $frame = $2;
            my $orientation = $3;
            $contig_DNA =~ s/(\_ENRICHMENT.*)//;
            $contig =~ /(\S+)_ENRICHMENT_control_(\S+)_selected_(\S+)_ratio_(\S+)/;
            my $contig_name_to_check = $1;
            my $control_reads = $2;
            my $enriched_reads = $3;
            my $read_count = $control_reads + $enriched_reads;

            $contig =~ /length_(\S+)_/;
            my $contig_length = $1;

            my $start = $tmp[17];
            my $end = $tmp[18];
            my $faa = $id2AA{$contig};
            my $seq = $faa->subseq($start, $end);

            # Filter by unwanted, contig length, and read count
            next if ($hash_contig && $$hash_contig{$contig_name_to_check});
            next if ($config->{contig_min} && $contig_length < $config->{contig_min});
            next if ($config->{read_count_min} && $read_count < $config->{read_count_min});

            my $id;
            if ($enrichment > $config->{cutoff}) {
                $selected++;
                $id = "S_" . $selected . "_" . $contig . "|group1";
            } else {
                $unselected++;
                $id = "unselected_" . $unselected . "_" . $contig . "|group2";
            }
            my $seq2 = $faa->seq;
            my ($start_domain, $end_domain) = match_positions($seq, $seq2);
            if ($start_domain && $end_domain) {
                my $AAseq = $faa->subseq($start_domain, $end_domain);
                my ($meth1, $meth2) = get_first_M($AAseq);
                $AAseq =~ s/\*//g;
                print ">$id\n$AAseq\n";
            }
        }
    }
    close $hit_fh;
}

sub get_first_M {
    my ($string) = @_;
    return if not $string =~ /[^M]*M/si;
    my $start = $-[0];
    my $end = $+[0] - 1;
    return ($start, $end);
}

sub match_positions {
    my ($regex, $string) = @_;
    return if not $string =~ /\*[^\*]*$regex[^\*]*\*/si;
    my $start = $-[0] + 1;
    my $end = $+[0];
    return ($start, $end);
}

sub parse_bed {
    my ($bed) = @_;
    my %hash;
    open(my $bed_fh, "<", $bed) or die "can't open $bed\n";
    while (my $line = <$bed_fh>) {
        chomp $line;
        my @tmp = split /\t/, $line;
        my $id = $tmp[0];
        $hash{$id}++;
    }
    close $bed_fh;
    return \%hash;
}
