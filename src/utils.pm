use strict;
use warnings;
use Cwd;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);
use Data::Dumper;
use Math::CDF qw(:all);
use Statistics::Descriptive;

package utils;
use Exporter qw(import);
our @EXPORT_OK = qw(make_dir 
                parse_pfam
                make_unique_path 
                process_fastq
                construct_paired_filename
                append_assemblies
                parse_pfam_file
                calculate_enrichment_stats
                parse_bed
                get_DNA);

sub main {
    # die "Usage: $0 dir1 dir2 ... output_dir\n" unless @ARGV >= 2;
    # my $out_dir = pop @ARGV;
    # my @dirs = @ARGV;
    my $out_dir = "appended_phage_test";
    my @dirs = ("phage_test", "phage_test_2", "phage_test_3");
    append_assemblies(@dirs, $out_dir);
}

### GENERAL USAGE
sub make_dir {
    my ($dir) = @_;
    system("mkdir -p $dir") unless -d $dir;
}

sub get_DNA {
    my ($contig, $DNA)=@_;
    my %id2seq;
    my $seq_in = Bio::SeqIO->new(-file   => $DNA,
				 -format => "fasta", );
    while ( my $faa = $seq_in->next_seq() ) {
	    my $tmp = $faa->id;
	    $tmp =~ /(.*)\_ratio\_.*/;
	    my $id = $1;
	    if ($$contig{$id}) {
	        my $seq = $faa->seq;
	        $id2seq{$id}=$seq;
	    }
    }
    return (\%id2seq)
}

### FOR PIPELINE WORKFLOW 
sub make_unique_path {
    my ($path, $lite, $assemble, $annotate, $map, $enrich) = @_;

    my $output_dir = "../output";
    make_dir($output_dir);

    my $basename = $path;
    $basename =~ s/.*\///g; 
    my $full_path = "$output_dir/$basename";

    if ($lite ne "none") {
        return $full_path;
    }
    
    my $counter = 1;
    my $unique_path = $full_path;
    while (-e $unique_path) {
        $unique_path = "../output/${path}_$counter";
        $counter++;
    }
    my $prev_path = $unique_path; $prev_path =~ s/_(\d+)$/'_' . ($1 - 1)/e;
    if ($assemble) {
        return $unique_path;
    }
    if ($map && -e $prev_path."/assembly") {
        return $prev_path;
    }
    if ($annotate && -e $prev_path."/assembly") {
        return $prev_path;
    }
    if ($enrich && -e $prev_path."/assembly" && -e $prev_path."/mapping") {
        return $prev_path;
    }
    return $unique_path;
}

sub process_fastq {
    my ($fq_1, $fq_2, $generic, $outdir, $trim) = @_;

    my $data_dir = "../data";
    system("mkdir -p $data_dir") unless -d $data_dir;

    # compress fq to fq.gz if it is not already compressed
    if ($fq_1 !~ /\.gz$/) {
        die "File $fq_1 does not exist!" unless -e $fq_1;
        my ($fq1_base) = $fq_1 =~ /([^\/]+)$/;  # get just the filename
        my $fq1_gz = "$data_dir/$fq1_base.gz";
        system("gzip -c $fq_1 > $fq1_gz") == 0 or die "Failed to compress $fq_1";
        $fq_1 = $fq1_gz;
    }

    # if fq_2 is empty string, then we generate it from fq_1
    if ($fq_2 eq "") {
        $fq_2 = $fq_1;
        $fq_2 =~ s/(?<=[._-])1(?=[._-]|$)/2/g;
    }

    if ($fq_2 !~ /\.gz$/) {
        die "File $fq_2 does not exist!" unless -e $fq_2;
        my ($fq2_base) = $fq_2 =~ /([^\/]+)$/;  # get just the filename
        my $fq2_gz = "$data_dir/$fq2_base.gz";
        system("gzip -c $fq_2 > $fq2_gz") == 0 or die "Failed to compress $fq_2";
        $fq_2 = $fq2_gz;
    }
    
    # run trim_galore if not another trimmer or none if empty string
    if ($trim eq "trim_galore") {
        my $command = "trim_galore --paired $fq_1 $fq_2 --basename $generic -o $outdir";
        system($command) == 0 or die "Failed to run trim_galore: $command";
        $fq_1 = "$generic.1_val_1.fq.gz";
        $fq_2 = "$generic.2_val_2.fq.gz";
    } else{
        # if no valid trimmer is specified, we assume the files are already trimmed we do nothing
        print "No valid trimmer specified, assuming files are already trimmed.\n";
    }
    
    return ($fq_1, $fq_2);
}

### FOR ASSEMBLY
sub construct_paired_filename {
    my ($read1_file) = @_;
    my $read2_file = $read1_file;
    $read2_file =~ s/1_val_1.fq.gz/2_val_2.fq.gz/;
    return $read2_file;
}

sub append_assemblies {
    my $out_dir = pop @_;
    my @dirs = @_;
    my $base_path = "../output/";
    my $appended_fai;
    my $appended_fasta;
    my $appended_faa;
    foreach my $dir (@dirs) {
        my $assembly_path = $base_path.$dir."/assembly/";
        opendir(my $dh, $assembly_path) or die "Cannot open directory $assembly_path: $!";
        my @files = readdir($dh);
        closedir($dh);
        foreach my $file (@files) {
            my $full_path = $assembly_path . $file;
            if ($file =~ /control_and_selected\.fasta$/) {
                $appended_fasta .= slurp_file($full_path);
            }
            if ($file =~ /control_and_selected\.faa$/) {
                $appended_faa .= slurp_file($full_path);
            }
            if ($file =~ /control_and_selected\.fasta\.fai$/) {
                $appended_fai .= slurp_file($full_path);
            }
        }
    }
    system("mkdir -p ../output/".$out_dir."/assembly") unless -d "../output/".$out_dir."/assembly";

    write_file("../output/".$out_dir."/assembly/appended_control_and_selected.fasta", $appended_fasta);
    write_file("../output/".$out_dir."/assembly/appended_control_and_selected.faa",   $appended_faa);
    write_file("../output/".$out_dir."/assembly/appended_control_and_selected.fasta.fai", $appended_fai);
}

sub parse_pfam {
    my ($pfam)=@_;
    my %id2pfam;
    open (PFAM, $pfam) or die "can['t open $pfam\n";
    foreach my $line (<PFAM>) {
        chomp $line;
        my @tmp = split/\s+/, $line;
        my $id = $tmp[0];

        my $pfam = $tmp[4];
        $id2pfam{$id}=$pfam;
    }
    close PFAM;
    return (\%id2pfam);
}

sub slurp_file {
    my ($filename) = @_;
    local $/ = undef;
    open(my $fh, "<", $filename) or die "Can't open $filename: $!";
    my $content = <$fh>;
    close $fh;
    return $content;
}

sub write_file {
    my ($filename, $content) = @_;
    open(my $fh, ">", $filename) or die "Can't write to $filename: $!";
    print $fh $content;
    close $fh;
}

sub downsample_fq {
    my ($fq1, $fq2, $downsample_val) = @_;

    my ($fq1_dir, $fq1_base) = ($fq1 =~ m|^(.*?)/?([^/]+)$|);
    my ($fq2_dir, $fq2_base) = ($fq2 =~ m|^(.*?)/?([^/]+)$|);

    $fq1_dir ||= ".";  # default to current directory if no path
    $fq2_dir ||= ".";

    # Remove file extension for output base name
    (my $fq1_stub = $fq1_base) =~ s/\.(fastq|fq|fasta|fa)(\.gz)?$//i;
    (my $fq2_stub = $fq2_base) =~ s/\.(fastq|fq|fasta|fa)(\.gz)?$//i;

    my $out1 = "$fq1_dir/${fq1_stub}_downsampled.fq";
    my $out2 = "$fq2_dir/${fq2_stub}_downsampled.fq";

    system("seqtk sample $fq1 $downsample_val > $out1") == 0 or die "Failed to downsample $fq1";
    system("seqtk sample $fq2 $downsample_val > $out2") == 0 or die "Failed to downsample $fq2";

    return ($out1, $out2);
}

### FOR ENRICHMENT

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
        my $total = $x + $y;
        my $proba = $n / ($n + $N);
        my $test = 1 - (Math::CDF::pbinom($x, $total, $proba));
        my $stats = $x."_".$total."_".$proba;
        $result4{$test}{$pfam} = $stats;
    }
    return \%result4;
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