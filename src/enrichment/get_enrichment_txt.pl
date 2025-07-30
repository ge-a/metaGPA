#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);
use Bio::SeqIO;

exit main();

sub main {
    my $config = parse_args();

    my $id2edgeR_pvalue = $config->{edgR} ? parse_edgR($config->{edgR}) : undef;

    my %bed = parse_enrichment(
        $config->{enrichment},
        $config->{control_bam},
        $config->{selection_bam},
        $id2edgeR_pvalue
    );

    write_enriched_fasta($config->{fasta}, \%bed, $config->{out});

    return 0;
}

sub parse_args {
    my %config;

    GetOptions(
        "fasta=s"        => \$config{fasta},
        "enrichment=s"   => \$config{enrichment},
        "out=s"          => \$config{out},
        "control-bam=s"  => \$config{control_bam},
        "selection-bam=s"=> \$config{selection_bam},
        "edgR=s"         => \$config{edgR}, 
        "help|h"         => \$config{help},
    ) or usage();

    usage() unless $config{fasta} && $config{enrichment} && $config{out};

    return \%config;
}

sub usage {
    print <<EOF;
Usage: perl $0 --fasta FASTA_FILE --enrichment BED_FILE --out OUTPUT_FILE --control-bam BAM --selection-bam BAM [options]

Required:
    --fasta FILE             Input FASTA file
    --enrichment FILE        BED file with enrichment regions
    --out FILE               Output file name
    --control-bam FILE       BAM file for the control sample
    --selection-bam FILE     BAM file for the selection sample

Optional:
    --edgR FILE              Output from edgeR differential analysis (optional)
    --help, -h               Show this help message

Example:
    perl $0 --fasta contigs.fa --enrichment enriched_regions.bed \\
        --out results.txt --control-bam ctrl.bam --selection-bam sel.bam \\
        --edgR differential_results.txt

EOF
    exit(1);
}

sub parse_enrichment {
    my ($enrichment, $control_bam, $selection_bam, $id2edgeR_pvalue) = @_;
    my %bed;
    open(my $cov_fh, "<", $enrichment) or die "Can't open $enrichment: $!";
    my $control_all = get_total_mapped_reads($control_bam);
    my $enriched_all = get_total_mapped_reads($selection_bam);
    while (my $line = <$cov_fh>) {
        chomp $line;
        my @tmp = split /\t/, $line;
        my $control = $tmp[-2];
        my $enriched = $tmp[-1];
        my $id = $tmp[0];
        my $length = $tmp[2];

        my $rpkm_enriched = (($enriched + 1) * 1e9) / ($enriched_all * $length);
        my $rpkm_control = (($control + 1) * 1e9) / ($control_all * $length);
        my $rpkm_ratio = ($rpkm_enriched) / ($rpkm_control);
        my $ratio = (($enriched+1)/($control+1)); # pseudocount
        $ratio  = $rpkm_ratio;

        my @info = ($control, $enriched, $ratio);

        if ($id2edgeR_pvalue) {
            my $logFC = $$id2edgeR_pvalue{$id}{"logFC"};
            my $Pvalue = $$id2edgeR_pvalue{$id}{"Pvalue"};
            push @info, $logFC, $Pvalue;
        }
        $bed{$id} = \@info;
    }
    close $cov_fh;
    return %bed;
}

sub write_enriched_fasta {
    my ($fasta, $bed_ref, $out) = @_;
    my %bed = %$bed_ref;

    my $seq_in = Bio::SeqIO->new(-format => 'fasta', -file => $fasta);

    # Print column headers
    print STDERR "Writing output with headers: id, control_count, selected_count, rpkm_ratio\n";
    open(my $info_fh, ">", $out) or die "Can't write to $out: $!";

    print $info_fh join("\t", "id", "control_count", "selected_count", "rpkm_ratio"), "\n";

    while (my $seq = $seq_in->next_seq()) {
        my $id = $seq->id;
        my $enrichment_ref = $bed{$id};
        my $seqstr = $seq->seq;
        $id =~ s/cov_.*//;
        my $enrichment = $enrichment_ref ? join("\t", @$enrichment_ref) : "";
        print $info_fh "$id\t$enrichment\n";
    }
    close $info_fh;
}

sub get_total_mapped_reads {
    my ($bam) = @_;
    my $cmd = "samtools view -F 260 -c $bam";
    my $count = `$cmd`;
    chomp $count;
    die "Failed to get mapped read count from $bam" unless $count =~ /^\d+$/;
    return $count;
}

sub parse_edgR {
    my ($file) = @_;
    my %result;
    open(my $edgr_fh, "<", $file) or die "can't open $file\n";
    while (my $line = <$edgr_fh>) {
        chomp $line;
        my @tmp = split /,/, $line;
        my $id = $tmp[0]; $id =~ s/\"//g;
        my $logFC = $tmp[1];
        my $logCPM = $tmp[2];
        my $Pvalue = $tmp[3];
        my $FDR = $tmp[4];
        $result{$id}{"logFC"} = $logFC;
        $result{$id}{"Pvalue"} = $Pvalue;
    }
    close $edgr_fh;
    return \%result;
}


#"","logFC","logCPM","PValue","FDR"
#"NODE_25476_length_852_selection",-7.26656621616605,0.310787993213469,3.99460840344187e-24,1.62286558841591e-18
#"NODE_4102_length_1780_control",5.91969495037133,2.51778352883633,9.91339319074969e-24,1.82118737820817e-18
#"NODE_9199_length_1336_control",7.26433015701448,0.51970903485849,2.20018006183038e-23,1.82118737820817e-18
#"NODE_286_length_3863_control",8.15170978605567,0.450453759594057,2.24426221242949e-23,1.82118737820817e-18
