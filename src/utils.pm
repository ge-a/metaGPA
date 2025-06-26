use strict;
use warnings;

package utils;
use Exporter qw(import);
our @EXPORT_OK = qw(make_dir 
                parse_pfam
                make_unique_path 
                process_fastq
                construct_paired_filename
                append_assemblies);

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

### FOR PIPELINE WORKFLOW 
sub make_unique_path {
    my ($path, $assemble, $annotate, $map, $enrich) = @_;

    my $output_dir = "../output";
    make_dir($output_dir);

    my $basename = $path;
    $basename =~ s/.*\///g; 
    my $full_path = "$output_dir/$basename";

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