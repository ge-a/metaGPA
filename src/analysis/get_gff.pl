#!/usr/bin/perl
use strict;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);
use Bio::SeqIO;

#this programs cleans the assemblies

my $error_sentence = "USAGE : perl $0 --pfam_hit pfamhit.tab --pfam PF000234 --fasta assembl.fasta (DNA) --cutoff 0.01 (default 3) --direction depleted (default enriched)\n";

# declare the options :
my $fasta; #fasta_file
my $pfam_hit; #position of pfam
my $pfam; #pfam family desired
my $out;
my $cutoff = 3;
my $direction = "enriched";
#get options :
GetOptions (    
                "fasta=s" => \$fasta,
		"pfam_hit=s" => \$pfam_hit,
		"pfam=s" => \$pfam,
                "cutoff=s" => \$cutoff,
                "direction=s" => \$direction
    ) or die $error_sentence;


#=================================                                                                                             
#if something went wrong in getting the option, notify :                                                                       
if (!$fasta || !$pfam_hit || !$pfam) {die $error_sentence}
#================================= 
my %id2seq;
my $seq_in = Bio::SeqIO->new(-file   => $fasta,
			     -format => "fasta", );

while ( my $fa = $seq_in->next_seq() ) {
    my $id = $fa->id;
    my $seq = $fa->seq;
    $id2seq{$id}=$seq;
}

my %contigs;
open (HIT, $pfam_hit) or die "can't open $pfam_hit\n";
foreach my $line (<HIT>) {
    chomp $line;
    $line =~ s/ +/&/g;
    
    my @tmp = split /\&/, $line;
    my $domain = $tmp[4];
    
    if ($domain eq $pfam) {
		my $contig = $tmp[0];
		$contig =~/.*_ratio_(.*)\-.*/;
		my $ratio = $1;
		my $selected =0;
		if ($direction eq "enriched") {if ($ratio > $cutoff) {$selected=1;}}
		if ($direction eq "depleted") {if ($ratio < $cutoff) {$selected=1;}}
		#print "direction $direction ratio  $ratio $selected\n";
		if ($selected >0) {
			$contig=~ s/_ratio_.*//g;
			$contig =~/(.*_selection)_ENRICHMENT.*/;
			my $contig_DNA = $1;
			if (!$contig_DNA) {
			$contig =~/(.*_control)_ENRICHMENT.*/;
			$contig_DNA = $1;
			}
			$contigs{$contig} = $contig_DNA;
			my $DNA = $id2seq{$contig_DNA};
			my $file_name = $contig_DNA.".fasta";
			
			open (OUT, ">$file_name") or die "can't open $file_name\n";
			print OUT ">$contig_DNA\n$DNA\n";
			close OUT;
		}
    }
}
close HIT;
my %result; my $i=0; my %function;
open (HIT, $pfam_hit) or die "can't open $pfam_hit\n";
foreach my $line (<HIT>)
{
    chomp $line;
    $line =~ s/ +/&/g;
    
    my @tmp = split /\&/, $line;
    my $contig = $tmp[0];
    my $generic_contig = $contig;
    $generic_contig=~ s/_ratio_.*//g;
    if ($contigs{$generic_contig}) {
		$i++;
		$contig=~ /(NODE_\d+)_.*/;
		my $contig_name = $1;
		$contig=~ /ratio_(.*)\-(.)(.)/;
		my $enrichment = $1;
		my $frame = $2;
		my $orientation = "+";
		my $sens = $3; 
		my $name = $tmp[3];
		my $start = $tmp[17] *3;
		my $end = $tmp[18] *3;
		my $domain = $tmp[4];
		my $score= $tmp[6];
		$contig=~ /.*length_(\d+)_.*/;
		my $length = $1;
		if ($sens eq "R"){
			$orientation = "-"; 
			my $old_start = $start; my $old_end= $end;
			$start = $length - $old_end; $end = $length - $old_start;
		}
		my $contig_DNA = $contigs{$generic_contig};
	
		if ($score < 0.0001) {
			$domain = $domain."_".$i;
			push @{$function{$name}}, $domain;
			push @{$result{$contig_DNA}},"$contig_DNA\tmetaGPA\tCDS\t$start\t$end\t$score\t$orientation\t$frame\tID=$domain";
		}
    }
}
close HIT;
open(FUNCTION, ">function.txt") or die "can't open file";
foreach my $name (keys %function)
{
    my @lines = @{$function{$name}};
    my $size = @lines;
    if ($size > 8) { 
		foreach my $l (@lines) {
			print FUNCTION "$l,$name\n";
		}
    }
}
close FUNCTION;
foreach my $contig (keys %result)
{
    my $gff3 = $contig.".gff3";
    open (FINAL, ">$gff3") or die "can't open $gff3\n";
    my @lines = @{$result{$contig}};
    foreach my $l (@lines)
    {
	print FINAL "$l\n";
    }
    close FINAL;
}