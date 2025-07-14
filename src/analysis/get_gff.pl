#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Bio::SeqIO;
use utils qw(parse_enrichment_info);

exit main();

sub main {
    my $config = parse_arguments();
    my %id2seq = parse_fasta_sequences($config->{fasta});
    my %contigs = select_contigs_by_pfam($config, \%id2seq);
    my ($result, $function) = extract_gff_features($config, \%contigs);
    write_function_file($function);
    write_gff_files($result);
    return 0;
}

sub parse_arguments {
    my %config;
    my $error_sentence = "USAGE : perl $0 --pfam_hit pfamhit.tab --pfam PF000234 --fasta assembly.fasta --enrichment_txt enrichment.txt --cutoff 3 --direction enriched\n";
    GetOptions(
        "fasta=s"         => \$config{fasta},
        "pfam_hit=s"      => \$config{pfam_hit},
        "pfam=s"          => \$config{pfam},
        "enrichment_txt=s"=> \$config{enrichment_txt},
        "cutoff=f"        => \$config{cutoff},
        "direction=s"     => \$config{direction},
        "help|h"          => \$config{help},
    ) or usage($error_sentence);

    $config{cutoff}    //= 3;
    $config{direction} //= "enriched";

    usage($error_sentence) if $config{help} || !$config{fasta} || !$config{pfam_hit} || !$config{pfam} || !$config{enrichment_txt};
    return \%config;
}

sub usage {
    my ($msg) = @_;
    print <<EOF;
$msg

Required:
    --fasta FILE            Assembly FASTA file (DNA)
    --pfam_hit FILE         PFAM hit tab file
    --pfam STRING           PFAM family desired
    --enrichment_txt FILE   Enrichment info file

Optional:
    --cutoff FLOAT          Enrichment cutoff (default: 3)
    --direction STRING      Enrichment direction: enriched or depleted (default: enriched)
    --help, -h              Show this help message

Example:
    $0 --pfam_hit pfamhit.tab --pfam PF000234 --fasta assembly.fasta --enrichment_txt enrichment.txt --cutoff 3 --direction enriched

EOF
    exit(1);
}
sub parse_fasta_sequences {
    my ($fasta) = @_;
    my %id2seq;
    my $seq_in = Bio::SeqIO->new(-file => $fasta, -format => "fasta");
    while (my $fa = $seq_in->next_seq()) {
        $id2seq{$fa->id} = $fa->seq;
    }
    return %id2seq;
}

# needs a helper function to pull info from an enrichment txt file

sub select_contigs_by_pfam {
    my ($config, $id2seq) = @_;
    my %contigs;
    my $enrichment_info = parse_enrichment_info($config->{enrichment_txt});
    open(my $hit_fh, "<", $config->{pfam_hit}) or die "can't open $config->{pfam_hit}\n";
    while (my $line = <$hit_fh>) {
        chomp $line;
        $line =~ s/ +/&/g;
        my @tmp = split /\&/, $line;
        my $domain = $tmp[4];
        next unless $domain eq $config->{pfam};
        my $contig = $tmp[0];
        $contig =~ /.*_ratio_(.*)\-.*/;
        my $ratio = $enrichment_info->{$contig}[2];
        my $selected = 0;
        $selected = 1 if ($config->{direction} eq "enriched" && $ratio > $config->{cutoff});
        $selected = 1 if ($config->{direction} eq "depleted" && $ratio < $config->{cutoff});
        next unless $selected;
        $contig =~ s/_ratio_.*//g;
        $contig =~ /(.*_selection)_ENRICHMENT.*/;
        my $contig_DNA = $1 || ($contig =~ /(.*_control)_ENRICHMENT.*/ ? $1 : undef);
        next unless $contig_DNA;
        $contigs{$contig} = $contig_DNA;
        my $DNA = $id2seq->{$contig_DNA};
        my $file_name = $contig_DNA . ".fasta";
        open(my $out_fh, ">", $file_name) or die "can't open $file_name\n";
        print $out_fh ">$contig_DNA\n$DNA\n";
        close $out_fh;
    }
    close $hit_fh;
    return %contigs;
}

sub extract_gff_features {
    my ($config, $contigs) = @_;
	my $enrichment_info = parse_enrichment_info($config->{enrichment_txt});
    my %result;
    my %function;
    my $i = 0;
    open(my $hit_fh, "<", $config->{pfam_hit}) or die "can't open $config->{pfam_hit}\n";
    while (my $line = <$hit_fh>) {
        chomp $line;
        $line =~ s/ +/&/g;
        my @tmp = split /\&/, $line;
        my $contig = $tmp[0];
        my $generic_contig = $contig; $generic_contig =~ s/_ratio_.*//g;
        if ($contigs->{$generic_contig}) {
			$i++;
			$contig =~ /(NODE_\d+)_.*/;
			my $contig_name = $1;
			$contig =~ /-(\d)([FR])$/;
			my $enrichment = $enrichment_info->{$contig}[2]
			my $frame = $1;
			my $orientation = "+";
			my $sens = $2;
			my $name = $tmp[3];
			my $start = $tmp[17] * 3;
			my $end = $tmp[18] * 3;
			my $domain = $tmp[4];
			my $score = $tmp[6];
			$contig =~ /.*length_(\d+)_.*/;
			my $length = $1;
			if ($sens eq "R") {
				$orientation = "-";
				my $old_start = $start; my $old_end = $end;
				$start = $length - $old_end; $end = $length - $old_start;
			}
			my $contig_DNA = $contigs->{$generic_contig};
			if ($score < 0.0001) {
				$domain = $domain . "_" . $i;
				push @{$function{$name}}, $domain;
				push @{$result{$contig_DNA}}, "$contig_DNA\tmetaGPA\tCDS\t$start\t$end\t$score\t$orientation\t$frame\tID=$domain";
			}
		}
    }
    close $hit_fh;
    return (\%result, \%function);
}

sub write_function_file {
    my ($function) = @_;
    open(my $function_fh, ">", "function.txt") or die "can't open function.txt\n";
    foreach my $name (keys %$function) {
        my @lines = @{$function->{$name}};
        my $size = @lines;
        if ($size > 8) {
            foreach my $l (@lines) {
                print $function_fh "$l,$name\n";
            }
        }
    }
    close $function_fh;
}

sub write_gff_files {
    my ($result) = @_;
    foreach my $contig (keys %$result) {
        my $gff3 = $contig . ".gff3";
        open(my $final_fh, ">", $gff3) or die "can't open $gff3\n";
        my @lines = @{$result->{$contig}};
        foreach my $l (@lines) {
            print $final_fh "$l\n";
        }
        close $final_fh;
    }
}