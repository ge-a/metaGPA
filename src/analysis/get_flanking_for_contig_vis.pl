#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Bio::SeqIO;
use utils qw(get_DNA);

exit main();

sub main {
    my $config = parse_arguments();

    my ($contigs, $contig2pfam) = get_contig_with_pfam($config->{pfam_hit}, $config->{pfam}, $config->{cutoff});
    my $DNA_contig = get_DNA($contigs, $config->{contig});

    print_contig_info($contigs, $contig2pfam, $DNA_contig, $config->{cutoff});

    return 0;
}

sub parse_arguments {
    my %config = (
        flanking => 2000,
        cutoff   => 3,
    );

    GetOptions(
        "pfam_hit=s" => \$config{pfam_hit},
        "pfam=s"     => \$config{pfam},
        "flanking=i" => \$config{flanking},
        "cutoff=f"   => \$config{cutoff},
        "contig=s"   => \$config{contig},
        "help|h"     => \$config{help},
    ) or usage();

    usage() if $config{help} || !$config{pfam_hit} || !$config{pfam} || !$config{contig};
    return \%config;
}

sub usage {
    print <<EOF;
Usage: $0 [options]

Required:
    --pfam_hit FILE        PFAM hit tab file (e.g., pfamhit.tab)
    --pfam STRING          PFAM family ID (e.g., PF000234)
    --contig FILE          DNA contig FASTA file

Optional:
    --flanking INT         Flanking region size (default: 2000)
    --cutoff FLOAT         Enrichment cutoff (default: 3)
    --help, -h             Show this help message

Example:
    $0 --pfam_hit pfamhit.tab \\
       --pfam PF000234 \\
       --contig DNA_contig_with_enrichment.fa \\
       --flanking 3000 \\
       --cutoff 3

EOF
    exit(1);
}

sub get_contig_with_pfam {
    my ($pfam_hit, $pfam, $cutoff) = @_;
    my %contig;
    my %contig2pfam;
    open(my $hit_fh, "<", $pfam_hit) or die "can't open $pfam_hit\n";
    while (my $line = <$hit_fh>) {
        chomp $line;
        $line =~ s/ +/&/g;
        my @tmp = split /\&/, $line;
        my $domain = $tmp[4];
        my $info = $tmp[0];
        $info =~ /(.*)_ratio_(.*)\-(.*)/;
        my $contig = $1;
        my $enrichment = $2;
        my $orientation = $3;
        my $start = $tmp[17];
        my $end = $tmp[18];
        my $pvalue = $tmp[6];
        my $name = $tmp[3];
        $contig2pfam{$contig}{$domain}{"start"} = $start;
        $contig2pfam{$contig}{$domain}{"end"} = $end;
        $contig2pfam{$contig}{$domain}{"orientation"} = $orientation;
        $contig2pfam{$contig}{$domain}{"name"} = $name;
        $contig2pfam{$contig}{$domain}{"pvalue"} = $pvalue;

        if ($domain eq $pfam) {
            $contig{$contig}{"start"} = $start;
            $contig{$contig}{"end"} = $end;
            $contig{$contig}{"orientation"} = $orientation;
            $contig{$contig}{"enrichment"} = $enrichment;
        }
    }
    close $hit_fh;
    return (\%contig, \%contig2pfam);
}

sub print_contig_info {
    my ($contigs, $contig2pfam, $DNA_contig, $cutoff) = @_;
    foreach my $contig (keys %$contigs) {
        my $enrichment = $contigs->{$contig}{"enrichment"};
        if ($enrichment > $cutoff) {
            print "$contig = $enrichment\n";
            my $fasta = $DNA_contig->{$contig};
            my $ref = $contig2pfam->{$contig};
            foreach my $pfam (sort { $ref->{$a}{"start"} <=> $ref->{$b}{"start"} } keys %$ref) {
                my $name = $ref->{$pfam}{"name"};
                my $start = $ref->{$pfam}{"start"};
                my $orientation = $ref->{$pfam}{"orientation"};
                my $pvalue = $ref->{$pfam}{"pvalue"};
                print "        $start $pfam = $name ($orientation) $pvalue\n";
            }
        }
    }
}