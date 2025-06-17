#!/usr/bin/perl
use strict;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);
use Bio::SeqIO;

#this programs cleans the assemblies

my $error_sentence = "USAGE : perl $0 --file assembly_file --tag control --out fileout --min_length 200 (default 500)\n";

# declare the options :

my $file; #assembly file from Metaspade
my $tag; #name_on the header of each contig
my $out; #output file
my $LENGHT = 500; #minimum size of the contig

#get options :
GetOptions (    "file=s" => \$file,
		"tag=s" => \$tag,
	        "out=s" => \$out,
	        "min_length=s" => \$LENGHT
    ) or die $error_sentence;

#=================================
#if something went wrong in getting the option, notify :
if (!$file || !$tag || !$out) {die $error_sentence}
#=================================

my $seq_in = Bio::SeqIO->new( -format => 'fasta',
                              -file   => $file,
    );

open (OUT, ">$out") or die "can't save in $file\n";
while ( my $seq = $seq_in->next_seq() ) {
    if($seq->length >= $LENGHT) {
	my $length = $seq->length;
	my $id = $seq->id;
	my $seq = $seq->seq;
	$id =~ s/cov_.*//;
	$id = $id.$tag;
	
	print OUT ">$id\n$seq\n";
    }
}



