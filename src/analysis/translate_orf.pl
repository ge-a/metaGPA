#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);
use Bio::SeqUtils;
use Bio::SeqIO;

exit main();

sub main {
    my $config = parse_args();

    my $DIR = getcwd;

    my $seq_in = Bio::SeqIO->new(
        -format => 'fasta',
        -file   => $config->{fasta},
    );

    open(my $out_fh, ">", $config->{out}) or die "Can't save into $config->{out}\n";
    while (my $seq = $seq_in->next_seq()) {5
        my $id = $seq->id();
        my $size = $seq->length();

        if ($size > $config->{min_size}) {
            my @seqs = Bio::SeqUtils->translate_6frames($seq);

            foreach my $prot (@seqs) {
                my $prot_id = $prot->id();
                my $prot_seq = $prot->seq();
                $prot_seq =~ s/^[A-Za-z]+\*//; # remove N-terminal junk before first stop
                $prot_seq =~ s/\*[A-Za-z]+$//; # remove C-terminal junk after last stop
                my @ORFS = split /\*/, $prot_seq;
                my $i = 0;
                foreach my $ORF (@ORFS) {
                    my $length = length($ORF);
                    if ($length > 20) {
                        $i++;
                        my $final_prot_id = $prot_id . "_" . $i;
                        print $out_fh ">$final_prot_id\n$ORF\n";
                    }
                }
            }
        }
    }
    close $out_fh;
    return 0;
}

sub parse_args {
    my %config = (
        min_size => 0,
        pfam     => "/mnt/home/ettwiller/laurence/database/pfam/Pfam-A.hmm",
    );
    my $error_sentence = "USAGE : perl $0 --fasta contig_file.fasta --out output_file [--min_size 100 (default 0)]\n";
    GetOptions(
        "fasta=s"   => \$config{fasta},
        "min_size=i"=> \$config{min_size},
        "out=s"     => \$config{out},
        "help|h"    => \$config{help},
    ) or usage($error_sentence);

    usage($error_sentence) if $config{help} || !$config{fasta} || !$config{out};
    return \%config;
}

sub usage {
    my ($msg) = @_;
    print <<EOF;
$msg
Required:
    --fasta FILE      Input nucleotide FASTA file
    --out FILE        Output file for translated ORFs

Optional:
    --min_size INT    Minimum contig size to process (default: 0)
    --help, -h        Show this help message

Example:
    perl $0 --fasta contig_file.fasta --out output_orfs.faa --min_size 100

EOF
    exit(1);
}