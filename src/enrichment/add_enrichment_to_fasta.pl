#!/usr/bin/perl
use strict;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);
use Bio::SeqIO;

#this programs cleans the assemblies

my $error_sentence = "USAGE : perl $0 --fasta fastafile --enrichment bedfile --out fileout --edgR edgeRfile\n";

# declare the options :

my $fasta; #assembly file combined (control and entriched) 
my $enrichment; #bed file containing the enrichment value
my $edgeR; #optional edgR file
my $out; #output file


#get options :
GetOptions (    "fasta=s" => \$fasta,
		"enrichment=s" => \$enrichment,
		"out=s" => \$out,
		"edgR=s" => \$edgeR
    ) or die $error_sentence;


#=================================                                                                                             
#if something went wrong in getting the option, notify :                                                                       
if (!$fasta || !$enrichment || !$out) {die $error_sentence}
#================================= 
my $id2edgeR_pvalue;
if ($edgeR)
{
    $id2edgeR_pvalue = parse_edgR($edgeR);
}

my %bed;
open (COVERAGE, $enrichment) or die;
foreach my $line (<COVERAGE>)
{
    chomp $line;
    my @tmp = split /\t/, $line;
    my $control = $tmp[-2];
    my $enriched = $tmp[-1];
    my $id = $tmp[0];
    my $ratio = (($enriched+1)/($control+1)); #add a speudocount
    my $info = "control_".$control."_selected_".$enriched."_ratio_".$ratio;

    if ($edgeR)
    {
	#$result{$id}{"logFC"}=$logFC;
        #$result{$id}{"Pvalue"}=$Pvalue;
	my $logFC = $$id2edgeR_pvalue{$id}{"logFC"};
	my $Pvalue = $$id2edgeR_pvalue{$id}{"Pvalue"};
	$info = $info."_logFC_".$logFC."_Pvalue_".$Pvalue;
    }
    
    $bed{$id}=$info;
}
close COVERAGE;


my $seq_in = Bio::SeqIO->new( -format => 'fasta',
                              -file   => $fasta,
    );

open (OUT, ">$out") or die "can't save in $out\n";
while (my $seq = $seq_in->next_seq() ) {
    my $id = $seq->id;
    my $enrichment = $bed{$id};
    my $seq = $seq->seq;
    $id =~ s/cov_.*//;
    $id = $id."_ENRICHMENT_".$enrichment;
    
    print OUT ">$id\n$seq\n";
}


sub parse_edgR {
    my ($file)=@_;
    my %result;
    open (EDGR, $file) or die "can't open $file\n";
    foreach my $line (<EDGR>)
    {
	chomp $line;
	my @tmp = split/,/, $line;
	my $id = $tmp[0]; $id=~ s/\"//g;
	my $logFC = $tmp[1];
	my $logCPM = $tmp[2];
	my $Pvalue = $tmp[3];
	my $FDR = $tmp[4];
	$result{$id}{"logFC"}=$logFC;
	$result{$id}{"Pvalue"}=$Pvalue;

    }
    close EDGR;
    return (\%result)
}

    
#"","logFC","logCPM","PValue","FDR"
#"NODE_25476_length_852_selection",-7.26656621616605,0.310787993213469,3.99460840344187e-24,1.62286558841591e-18
#"NODE_4102_length_1780_control",5.91969495037133,2.51778352883633,9.91339319074969e-24,1.82118737820817e-18
#"NODE_9199_length_1336_control",7.26433015701448,0.51970903485849,2.20018006183038e-23,1.82118737820817e-18
#"NODE_286_length_3863_control",8.15170978605567,0.450453759594057,2.24426221242949e-23,1.82118737820817e-18

