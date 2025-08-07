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
                get_control_enriched_reads
                parse_enrichment_info
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

    my $output_dir = "output";
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
        $unique_path = "output/${path}_$counter";
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

    my $data_dir = "data";
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

sub get_control_enriched_reads {
    my ($contig, $enrichment_info) = @_;
    my $control_reads   = $enrichment_info->{$contig}[0];
    my $enriched_reads  = $enrichment_info->{$contig}[1];
    my $ratio           = $enrichment_info->{$contig}[2];
    return ($control_reads, $enriched_reads, $ratio);
}

sub parse_pfam_file {
    my (
        $pfam_file_path,         # Path to Pfam output file
        $enrichment_file_path,   # Path to file with read enrichment info
        $excluded_contigs_ref,   # Hashref of contig names to exclude
        $min_read_count,         # Minimum total read count threshold
        $min_contig_length,      # Minimum contig length threshold
        $pfam_evalue_cutoff,     # Maximum acceptable Pfam e-value
        $enrichment_cutoff,      # Enrichment ratio cutoff
        $direction,               # 'enriched' or 'depleted' based on how ratio should be interpreted
        $multiselection,         # whether or not we are running a multiselection
    ) = @_;
    $multiselection //= 0;
    my (%pfam_to_status_contigs, %status_to_contigs, %pfam_to_description);
    my $enrichment_data = $multiselection
        ? parse_multi_enrichment_info($enrichment_file_path)
        : parse_enrichment_info($enrichment_file_path);
    open(my $pfam_fh, "<", $pfam_file_path) or die "Can't open Pfam file: $pfam_file_path\n";

    while (my $line = <$pfam_fh>) {
        chomp $line;
        my @fields = split /\s+/, $line;
        next unless $fields[0] =~ /^NODE_\d+_length_\d+_(?:selection|control)-\w+$/;
        my $full_contig_name = $fields[0];
        my $pfam_id          = $fields[4];
        my $pfam_description = $fields[3];  
        my $pfam_evalue      = $fields[6];
        # Normalize contig name to group by selection/control
        my $normalized_contig_name = $full_contig_name;
        $normalized_contig_name =~ s/((?:selection|control)).*$/$1/;
        if ($excluded_contigs_ref && $$excluded_contigs_ref{$full_contig_name}) {
            print STDERR "$full_contig_name is removed (excluded)\n";
            next;
        }
        my ($control_reads, $selection_reads, $enrichment_ratio) = get_control_enriched_reads($normalized_contig_name, $enrichment_data);
        my $total_read_count = $control_reads + $selection_reads;
        my $contig_length;
        if ($normalized_contig_name =~ /length_(\d+)_/) {
            $contig_length = $1;
        } else {
            warn "Could not extract contig length from $normalized_contig_name\n";
            next;
        }
        if ($total_read_count > $min_read_count && $contig_length > $min_contig_length && $pfam_evalue < $pfam_evalue_cutoff) {
            my $status;
            if ($direction eq "depleted") {
                # If direction is "depleted", low ratio is "enriched"
                $status = ($enrichment_ratio < $enrichment_cutoff) ? "enriched" : "depleted";
            } else {
                # Default: high ratio is "enriched"
                $status = ($enrichment_ratio >= $enrichment_cutoff) ? "enriched" : "depleted";
            }
            # Store the filtered result in structured hashes
            $pfam_to_status_contigs{$pfam_id}{$status}{$normalized_contig_name}++;
            $status_to_contigs{$status}{$normalized_contig_name}++;
            $pfam_to_description{$pfam_id} = $pfam_description;
        }
    }
    close $pfam_fh;
    return (\%pfam_to_status_contigs, \%status_to_contigs, \%pfam_to_description);
}

sub calculate_enrichment_stats {
    my ($pfam_to_status_contigs_ref, $status_to_contigs_ref) = @_;
    my %pvalue_to_pfam_stats;

    # Total enriched and depleted contigs globally
    my $global_enriched_contigs = $status_to_contigs_ref->{"enriched"} || {};
    my $global_depleted_contigs = $status_to_contigs_ref->{"depleted"} || {};

    my $total_enriched_global = keys %$global_enriched_contigs;
    my $total_depleted_global = keys %$global_depleted_contigs;
    my $total_global = $total_enriched_global + $total_depleted_global;

    # Global background enrichment probability (used for binomial test)
    my $global_enrichment_probability = $total_global > 0
        ? $total_enriched_global / $total_global
        : 0.5;  # fallback to neutral probability if total is zero

    foreach my $pfam_id (keys %$pfam_to_status_contigs_ref) {
        my $enriched_contigs = $pfam_to_status_contigs_ref->{$pfam_id}{"enriched"} || {};
        my $depleted_contigs = $pfam_to_status_contigs_ref->{$pfam_id}{"depleted"} || {};

        # Count of enriched and depleted contigs for this Pfam domain
        my $enriched_count = keys %$enriched_contigs;
        my $depleted_count = keys %$depleted_contigs;

        $enriched_count++; # pseudocounts
        $depleted_count++;
        my $total_contigs_with_pfam = $enriched_count + $depleted_count;

        # Local ratio: proportion of enriched contigs for this Pfam
        my $local_enrichment_ratio = $enriched_count / $total_contigs_with_pfam;

        # Binomial test: p-value for observing ≥ enriched_count enriched contigs
        my $p_value = 1 - Math::CDF::pbinom($enriched_count - 1, $total_contigs_with_pfam, $global_enrichment_probability);
        my $stats_string = join("_", $enriched_count, $total_contigs_with_pfam, sprintf("%.4f", $local_enrichment_ratio));
        $pvalue_to_pfam_stats{$p_value}{$pfam_id} = {
            stats => $stats_string,
            global_enrichment_probability => sprintf("%.4f", $global_enrichment_probability)
        };
    }
    return \%pvalue_to_pfam_stats;
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
        next if $line =~ /^id\b/i;          # skip header line
        my ($id, @values) = split /\t/, $line;
        $enrichment_hash{$id} = \@values;
    }
    close $fh;
    return \%enrichment_hash;
}

sub parse_multi_enrichment_info {
    my ($enrichment_arrayref) = @_;
    my %enrichment_hash;

    foreach my $entry (@$enrichment_arrayref) {
        my $id = $entry->{id};
        my $control_count = $entry->{control_count};
        my $selection_count = $entry->{selection_count};
        my $selection_ratio = $entry->{selection_ratio};

        $enrichment_hash{$id} = [$control_count, $selection_count, $selection_ratio];
    }

    return \%enrichment_hash;
}