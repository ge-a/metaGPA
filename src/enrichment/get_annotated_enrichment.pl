use strict;
use warnings;
use Cwd;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);
use Data::Dumper;
use Math::CDF qw(:all);

exit main();

sub main {
    my ($pfam, $enrichment_file_path, $out, $cutoff, $direction, $read_count_min, $contig_min, $unwanted, $pfam_cutoff, $TOTAL) = parse_args();

    my $hash_contig = $unwanted ? parse_bed($unwanted) : undef;

    my ($result2, $result3, $pfam2desc) = parse_pfam_file($pfam, $enrichment_file_path, $hash_contig, $read_count_min, $contig_min, $pfam_cutoff, $cutoff, $direction);
    my $result4 = calculate_enrichment_stats($result2, $result3);

    open(my $out_fh, ">", $out) or die "Can't open $out\n";
    write_enrichment_output($result4, $pfam2desc, $TOTAL, $out_fh);
    close $out_fh;

    return 0;
}

sub parse_pfam_file {
    my ($pfam, $enrichment_file_path, $hash_contig, $read_count_min, $contig_min, $pfam_cutoff, $cutoff, $direction) = @_;
    my (%result2, %result3, %pfam2desc);
    my $enrichment_info = parse_enrichment_info($enrichment_file_path);
    open(my $pfam_fh, "<", $pfam) or die "Can't open $pfam\n";
    while (my $line = <$pfam_fh>) {
        chomp $line;
        my @tmp = split /\s+/, $line;   
        if ($tmp[0] !~ /^NODE_\d+_length_\d+_(?:selection|control)-\w+$/) {
            next;
        }
        my $contig = $tmp[0];
        my $pfam_Evalue = $tmp[6]; 
        my $contig_name_to_check = $contig;
        $contig =~ s/((?:selection|control)).*$/$1/;
        my $control_reads = $enrichment_info->{$contig}[0];
        my $enriched_reads = $enrichment_info->{$contig}[1];
        my $read_count = $control_reads + $enriched_reads;
        my $ratio = $enrichment_info->{$contig}[2]; # Confirm if this is what is being referenced for creating an enrichment score distribution
        $contig =~ /length_(\S+)_/;
        my $contig_length = $1;
        my $status ="rejected";
        if ($hash_contig && $$hash_contig{$contig_name_to_check}) {
            print STDERR "$contig_name_to_check is removed\n";
        }
        elsif ($read_count > $read_count_min && $contig_length > $contig_min && $pfam_Evalue < $pfam_cutoff) {
            if ($direction eq "depleted") {
                $status = ($ratio < $cutoff) ? "enriched" : "depleted";
            } else {
                $status = ($ratio >= $cutoff) ? "enriched" : "depleted";
            }
            my $pfam = $tmp[4];  
	        my $hit = $tmp[6];
            my $description = $tmp[3];
            $result2{$pfam}{$status}{$contig}++;
            $result3{$status}{$contig}++;
            $pfam2desc{$pfam}=$description; 
            #print STDERR " $contig $ratio $status $pfam $hit\n";
	        #my $description = my $str = join('_', @tmp[22 .. $#tmp]);
        }
    }
    close $pfam_fh;
    return (\%result2, \%result3, \%pfam2desc);
}

sub calculate_enrichment_stats {
    my ($result2, $result3) = @_;
    my %result4;

    my $Ns = $result3->{"depleted"};
    my $ns = $result3->{"enriched"};
    my $n = keys %$ns;
    my $N = keys %$Ns;

    foreach my $pfam (keys %$result2) {
        my $xs = $result2->{$pfam}{"enriched"};
        my $ys = $result2->{$pfam}{"depleted"};
        my $x = keys %$xs; $x++; # pseudocount
        my $y = keys %$ys; $y++; # pseudocount
        #print STDERR " $xs $ys $x $y\n";
        my $total = $x + $y;
        my $proba = $n / ($n + $N);
        my $test = 1 - (pbinom($x, $total, $proba));
        #print "$pfam proba : $proba test = $test $n and $N $x $y\n";
        my $stats = $x."_".$total."_".$proba;
        $result4{$test}{$pfam} = $stats;
    }
    return \%result4;
}

sub write_enrichment_output {
    my ($result4, $pfam2desc, $TOTAL, $out_fh) = @_;
    foreach my $e (sort { $a <=> $b } keys %$result4) {
        my $pfams = $result4->{$e};
        foreach my $pfam (keys %$pfams) {
            my $x = $pfams->{$pfam};
            my ($n1, $n2, $freq) = split /\_/, $x;
            my $tot = $n1 + $n2;
            if ($tot > $TOTAL) {
                my $description = $pfam2desc->{$pfam};
                print $out_fh "$e $pfam $x $description\n";
            }
        }
    }
}

sub parse_args {
    my $error_sentence = "USAGE : perl $0 --pfam pfam_tab --out output_file OPTIONAL : --cutoff 10 (default 3) --direction depleted (default enriched) --read_count_min 100 (total number of reads in control+enriched, default 0) --contig_min 500 (length of the contig in bp, default 0) --unwanted contif.bed (defaut NONE) --pfam_cutoff 0.0000001 (this is the pfam Evalue. Anything matching of below this cutoff is used, DEFAULT no cutoff) --pfam_number 10 (at least 10 instances of a particular pfam, DEFAULT=0)";
    my ($pfam, $enrichment_file_path, $out, $cutoff, $direction, $read_count_min, $contig_min, $unwanted, $pfam_cutoff, $TOTAL);
    $cutoff = 0;
    $direction = "enriched";
    $read_count_min = 0;
    $contig_min = 0;
    $pfam_cutoff = 100000000000000;
    #from : http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.gz
    $pfam = "/mnt/home/ettwiller/laurence/database/pfam/Pfam-A.hmm";
    $TOTAL = 0;

    GetOptions(
        "pfam=s"          => \$pfam, # hmmer format file (tab delimited).
        "enrichment_txt=s"=> \$enrichment_file_path, # enrichment file path.
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
    return ($pfam, $enrichment_file_path, $out, $cutoff, $direction, $read_count_min, $contig_min, $unwanted, $pfam_cutoff, $TOTAL);
}

sub parse_enrichment_info {
    my ($filename) = @_;
    my %enrichment_hash;

    open(my $fh, "<", $filename) or die "Can't open $filename: $!";
    while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ /^\s*$/; # skip empty lines
        my ($id, @values) = split /\t/, $line;
        $enrichment_hash{$id} = \@values;
    }
    close $fh;
    return \%enrichment_hash;
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

sub log10 {
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
