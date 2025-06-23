#!/usr/bin/perl
use strict;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);
use Bio::SeqIO;

exit main();

sub main {
    my ($file, $tag, $out, $LENGTH) = parse_args();

    clean_assembly($file, $tag, $out, $LENGTH);

	return 0;
}

#this programs cleans the assemblies
sub parse_args {
    my $error_sentence = "USAGE : perl $0 --file assembly_file --tag control --out fileout --min_length 200 (default 500)\n";
    my ($file, $tag, $out, $LENGTH);
    $LENGTH = 500; # default

    GetOptions(
        "file=s"       => \$file,
        "tag=s"        => \$tag,
        "out=s"        => \$out,
        "min_length=s" => \$LENGTH
    ) or die $error_sentence;

    die $error_sentence unless $file && $tag && $out;
    return ($file, $tag, $out, $LENGTH);
}

sub clean_assembly {
    my ($file, $tag, $out, $LENGTH) = @_;

    my $seq_in = Bio::SeqIO->new(-format => 'fasta', -file => $file);
    open(my $out_fh, ">", $out) or die "can't save in $out\n";

    while (my $seq = $seq_in->next_seq()) {
        if ($seq->length >= $LENGTH) {
            my $id = $seq->id;
            my $seqstr = $seq->seq;
            $id =~ s/cov_.*//;
            $id = $id.$tag;
            print $out_fh ">$id\n$seqstr\n";
        }
    }
    close $out_fh;
}
