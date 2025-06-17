use strict;
use Cwd;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);
use Data::Dumper;
use Math::CDF qw(:all);


my $error_sentence = "USAGE : perl $0 --pfam pfam_tab --out output_file OPTIONAL : --cutoff 10 (default 3) --direction depleted (default enriched) --read_count_min 100 (total number of reads in control+enriched, default 0) --contig_min 500 (length of the contig in bp, default 0) --unwanted contif.bed (defaut NONE) --pfam_cutoff 0.0000001 (this is the pfam Evalue. Anything matching of below this cutoff is used, DEFAULT no cutoff) --pfam_number 10 (at least 10 instances of a particular pfam, DEFAULT=0)";

# declare the options upfront :
my $pfam;
my $out;
my $cutoff = 0;
my $direction = "enriched";
my $read_count_min =0;
my $contig_min =0; #minimum contig length
my $unwanted;
my $pfam_cutoff =100000000000000;
#from : http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.gz
my $pfam = "/mnt/home/ettwiller/laurence/database/pfam/Pfam-A.hmm";
my $TOTAL =0;
#get options :
GetOptions ( 
    "pfam=s" => \$pfam, # hmmer format file (tab delimited).
    "out=s" => \$out, #output file.
    "cutoff=s" => \$cutoff, #cutoff for the enrichment
    "direction=s" => \$direction, #depeted or enriched
    "read_count_min=s" => \$read_count_min, #minimum (total) read number
    "contig_min=s" => \$contig_min, #minimum size contig (bp)
    "unwanted=s" => \$unwanted,
    "pfam_cutoff=s" => \$pfam_cutoff,
    "pfam_number=s" => \$TOTAL
    ) or die $error_sentence;

#=================================
#if something went wrong in getting the option, notify :
if (!$out || !$pfam) {die $error_sentence}
#=================================
my %result2; my %result3;my %result4;my %pfam2desc;

my $hash_contig;
if (defined $unwanted) {
    $hash_contig = parse_bed($unwanted);
}


open (PFAM, $pfam) or die;
open(OUT, ">$out") or die;
foreach my $line (<PFAM>)
{
    chomp $line;
    my @tmp = split /\s+/, $line;
    my $contig = $tmp[0];
    my $pfam_Evalue = $tmp[6];
    #NODE_45532_length_1108_all_ENRICHMENT_control_96_selected_4_ratio_0.0515463917525773-1R
   # NODE_1_length_44099_selection_ENRICHMENT_control_7170_selected_32717_ratio_4.562543578-2F
    $contig =~ s/-\S+//; 
    $contig =~ /(\S+)_ENRICHMENT_control_(\S+)_selected_(\S+)_ratio_(\S+)/;
    my $contig_name_to_check = $1;
    my $control_reads = $2;
    my $enriched_reads = $3;
    my $read_count = $control_reads + $enriched_reads;
    my $ratio = $4;
    $contig =~ /length_(\S+)_/;
    my $contig_length = $1;
    my $status ="rejected";
    
    if ($$hash_contig{$contig_name_to_check}) #unwanted contif are not present
    {
	print STDERR "$contig_name_to_check is removed\n";
    }
    
    elsif ($read_count > $read_count_min && $contig_length > $contig_min && $pfam_Evalue < $pfam_cutoff) {
	if ($direction eq "depleted")
	{
	    if ($ratio < $cutoff) {$status = "enriched";}
	    if ($ratio >= $cutoff) {$status = "depleted";}	
	}
	else { #if enriched direction
	    
	    if ($ratio >= $cutoff) {$status = "enriched";}	
	    if ($ratio < $cutoff) {$status = "depleted";}
	}
	
	my $pfam = $tmp[4];  
	my $hit = $tmp[6];
	
#	print STDERR " $contig $ratio $status $pfam $hit\n";
	
	
	    #	my $description = my $str = join('_', @tmp[22 .. $#tmp]);
	my $description = $tmp[3];
	$result2{$pfam}{$status}{$contig}++;
	$pfam2desc{$pfam}=$description; 
	$result3{$status}{$contig}++;
	
    }

}

my $Ns = $result3{"depleted"};
my $ns = $result3{"enriched"};
my $n = keys %$ns;
my $N = keys %$Ns;


foreach my $pfam (keys %result2){
    
    my $xs = $result2{$pfam}{"enriched"};
    my $ys = $result2{$pfam}{"depleted"};
    my $x = keys %$xs; $x++; #speudocount
    my $y = keys %$ys; $y++; #speudocount
#    print STDERR " $xs $ys $x $y\n";
    my $total = $x + $y;
    my $proba = $n / ($n + $N);
    my $test = 1-(pbinom($x, $total, $proba));
#    print "$pfam proba : $proba test = $test $n and $N $x $y\n";
    my $stats = $x."_".$total."_".$proba;
    
    $result4{$test}{$pfam} = $stats;
     
}

foreach my $e (sort {$a <=> $b} keys %result4){
    
    my $pfams = $result4{$e};
    foreach my $pfam (keys %$pfams){
	my $x = $$pfams{$pfam};
	my ($n1, $n2, $freq) = split /\_/, $x;
	my $tot = $n1 +$n2;
	if ($tot > $TOTAL){
	    my $description = $pfam2desc{$pfam};
	    
	    print OUT "$e $pfam $x $description\n";
	}
    }
}


sub parse_bed {
    my ($bed) = @_;
    my %hash;
    open (BED, $bed) or die "can't open $bed\n";
    foreach my $line (<BED>)
    {
	chomp $line;
	my @tmp = split/\t/, $line;
	my $id = $tmp[0];
	$hash{$id}++;
    }

    close BED;
    return \%hash;
}



sub log10 
{
    my $n = shift;
      
    # using pre-defined log function
    return log($n) / log(10);
}


sub hypergeometric {

    my ($M,$p,$F,$n) = @_;
    return unless $n>0 && $n == int($n) && $p > 0 && $p == int($p) &&
	$M > 0 && $M <= $n+$p;
    return 0 unless $p <= $M && $p == int($p);

    return pari2num((binomial($M,$p) * binomial($F-$M, $n-$p)
		     / binomial($F,$n))); 

}
