#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use Getopt::Long qw(GetOptions);
use Bio::SeqUtils;
use Bio::SeqIO;
use FindBin qw($Bin);
use lib File::Spec->catdir($Bin, "..");
use utils qw(parse_pfam);

exit main();

sub main {
    my $config = parse_arguments();

    my $id2pfam = parse_pfam($config->{pfam});
    my $seq_in = Bio::SeqIO->new(-format => 'fasta', -file => $config->{fasta});
    open(my $out_fh, ">", $config->{out}) or die "Can't save into $config->{out}\n";

    filter_and_write_sequences(
        $seq_in, $out_fh, $id2pfam, $config->{min_size}, $config->{enrichment_cutoff}, $config->{reads_cutoff}
    );

    close $out_fh;
    return 0;
}

sub parse_arguments {
    my %config = (
        min_size          => 0,
        enrichment_cutoff => 3,
        reads_cutoff      => 100,
    );

    GetOptions(
        "fasta=s"     => \$config{fasta},
        "pfam=s"      => \$config{pfam},
        "min_size=i"  => \$config{min_size},
        "out=s"       => \$config{out},
        "enrichment_cutoff=f" => \$config{enrichment_cutoff},
        "reads_cutoff=i"      => \$config{reads_cutoff},
        "help|h"      => \$config{help},
    ) or usage();

    usage() if $config{help} || !$config{fasta} || !$config{pfam} || !$config{out};
    return \%config;
}

sub usage {
    print <<EOF;
Usage: $0 [options]

Required:
    --fasta FILE               Input ORF FASTA file (e.g., individual_ORF.faa)
    --pfam FILE                PFAM result tab file (e.g., pfam_result.tab)
    --out FILE                 Output FASTA file (e.g., pfam_ORF.faa)

Optional:
    --min_size INT             Minimum ORF size to include (default: 0)
    --enrichment_cutoff FLOAT  Minimum enrichment value (default: 3)
    --reads_cutoff INT         Minimum reads value (default: 100)
    --help, -h                 Show this help message

Example:
    $0 --fasta individual_ORF.faa \\
       --pfam pfam_result.tab \\
       --out pfam_ORF.faa \\
       --min_size 50

EOF
    exit(1);
}

sub filter_and_write_sequences {
    my ($seq_in, $out_fh, $id2pfam, $min_size, $enrichment_cutoff, $reads_cutoff) = @_;

    while (my $seq = $seq_in->next_seq()) {
        my $id = $seq->id();
        $id =~ /selected_(\S+)_ratio_(\S+)-/;
        my $reads = $1;
        my $enrichment = $2;
        my $size = $seq->length();
        my $seq_final = $seq->seq();

        if ($size > $min_size && $$id2pfam{$id} && $enrichment > $enrichment_cutoff && $reads > $reads_cutoff) {
            my $pfam = $$id2pfam{$id};
            print $out_fh ">$id $pfam $enrichment\n$seq_final\n";
        }
    }
}