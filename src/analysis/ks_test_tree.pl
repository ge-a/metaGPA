#!/usr/bin/perl
use strict;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);
use Bio::SeqIO;
use Bio::TreeIO;
use Statistics::TTest;
use utils qw(get_control_enriched_reads
            parse_enrichment_info);


my $error_sentence = "USAGE : perl $0 --fasta fastafile.faa --pfam_hit pfamhit.tab --pfam PF000234\n";

exit main();

sub main {
    my $config = parse_arguments();
    run_ks_test_tree($config);
    return 0;
}

sub parse_arguments {
    my %config = (
        out   => "test",
        cutoff => 3,
    );
    GetOptions(
        "fasta=s"       => \$config{fasta}, 
        "enrich_txt=s"  => \$config{enrich},
        "pfam_hit=s"    => \$config{pfam_hit},
        "pfam=s"        => \$config{pfam},
        "out=s"         => \$config{out},
        "cutoff=f"      => \$config{cutoff},
        "help|h"        => \$config{help},
    ) or usage();

    usage() if $config{help} || !$config{fasta} || !$config{pfam_hit} || !$config{pfam};
    return \%config;
}

sub usage {
    print <<EOF;
Usage: $0 [options]

Required:
    --fasta FILE        Combined assembly FASTA file (control and enriched) with enrichment scores
    --enrich_txt FILE   Enrichment_info.txt file
    --pfam_hit FILE     PFAM hit tab file with position of pfam
    --pfam STRING       PFAM family ID (e.g., PF000234) desired

Optional:
    --out STRING        Output prefix (default: test)
    --cutoff FLOAT       Enrichment cutoff (default: 3)
    --help, -h          Show this help message

Example:
    $0 --fasta combined.faa \\
       --pfam_hit pfamhit.tab \\
       --pfam PF000234 \\
       --out test \\
       --cutoff 3

EOF
    exit(1);
}

sub run_ks_test_tree {
    my ($config) = @_;

    my %id2seq;
    my $seq_in = Bio::SeqIO->new(-file => $config->{fasta}, -format => "fasta");
    while (my $faa = $seq_in->next_seq()) {
        my $id = $faa->id;
        $id2seq{$id} = $faa;
    }

    my ($unselected, $selected) = (0, 0);
    my $out1 = $config->{out} . "_not_aligned.fasta";
    my $out2 = $config->{out} . "_aligned.fasta";
    my $out3 = $config->{out} . "_aligned.tree";
    my $out4 = $config->{out} . "_mapping.txt";
    my $out5 = $config->{out} . "_TT.txt";

    open(my $hit_fh, "<", $config->{pfam_hit}) or die "can't open $config->{pfam_hit}\n";
    open(my $out1_fh, ">", $out1) or die "can't open $out1\n";
    open(my $out4_fh, ">", $out4) or die "can't open $out4\n";
    print $out4_fh "name\tleaf_dot_color\tleaf_label_color\tbar1_height\tbar1_gradient\n";

    my @array_enrichment;
    my $enrichment_info = parse_enrichment_info($config->{enrich});
    while (my $line = <$hit_fh>) {
        #Ribonuc_red_lgC PF02867.18 527 NODE_1_length_63839_selection_ENRICHMENT_control_10_selected_18765_ratio_1706-0F - 21279  3.9e-100  336.0   0.0   1   2   1.2e-67   1.9e-64  218.2   0.0     2   343  8420  8733  8419  8737 0.93 Ribonucleotide reductase, barrel domain
        chomp $line;
        $line =~ s/ +/&/g;
        my @tmp = split /\&/, $line;
        my $domain = $tmp[4];

        if ($domain eq $config->{pfam}) {
            my $contig = $tmp[0];
            $contig =~ /ratio_(.*)\-/;
            my $v = $1;
            my $enrichment = sprintf("%.1f", $v);
            $enrichment = 0.1 if $enrichment == 0;
            $contig =~ s/((?:selection|control)).*$/$1/;
            my $contig_name_to_check = $contig;
            my ($control_reads, $enriched_reads) = get_control_enriched_reads($contig, $enrichment_info);
            my $read_count = $control_reads + $enriched_reads;
            $contig =~ /length_(\S+)_/;
            my $contig_length = $1;
            my $start = $tmp[17];
            my $end = $tmp[18];
            my $faa = $id2seq{$contig};
            my $seq = $faa->subseq($start, $end);
            my $diff = $end - $start;
            my $id;
            my $log_enrichment = log10($enrichment);
            if ($enrichment > $config->{cutoff}) {
                $selected++;
                $id = "S_" . $selected;
                my $name = $id . "_" . $enrichment;
                print $out1_fh ">$name\n$seq\n";
                print $out4_fh "$name\tbp_green\tptm_rose\t$enrichment\tPurples\n";
            } else {
                $unselected++;
                $id = "unselected_" . $unselected;
                my $name = $id . "_" . $enrichment;
                print $out1_fh ">$name\n$seq\n";
                print $out4_fh "$name\tk_grey\tptm_sand\t$enrichment\tPurples\n";
            }
            push @array_enrichment, $enrichment;
        }
    }
    close $hit_fh;
    close $out1_fh;
    close $out4_fh;

    my $command_alignment = "mafft --maxiterate 1000 --localpair $out1 > $out2";
    system($command_alignment);
    my $command_tree = "FastTree -gamma $out2 > $out3";
    system($command_tree);

    if ($selected > 1 && $unselected > 1) {
        open(my $out5_fh, ">", $out5) or die "can't open $out5\n";
        parse_tree($out3, \@array_enrichment, $selected, $unselected, $out5_fh, $config->{cutoff});
        close $out5_fh;
    }
}

sub parse_tree {
    my ($file, $array1, $total_selected, $total_unselected)=@_;
    my $treeio = new Bio::TreeIO(-file => $file, -format => "newick");
   
    while(my $tree = $treeio->next_tree) {
	    for my $node ($tree->get_nodes) {
            my $node_id = $node->id;
            my ($count_above_cutoff, $count_below_cutoff) = 0;
            if (!$node->is_Leaf) {
                my @array2;
                for my $child ( $node->get_all_Descendents ) {
                    my $child_id = $child->id;
                    my $child_id_name = $child_id;
                    $child_id =~ s/(.*)\_//g;
                    if ($child->is_Leaf) {
                        push @array2, $child_id;
                        if ($child_id >= $CUTOFF) {$count_above_cutoff++;}
                        if ($child_id < $CUTOFF) {$count_below_cutoff++;}
                        # print "$child_id\n";
                    }
                }
                my $fraction_above_in_leaf = ($count_above_cutoff / $total_selected)*100;
                my $fraction_below_in_leaf = ($count_below_cutoff / $total_unselected)*100;

                # if ($fraction_above_in_leaf > 50 || $fraction_below_in_leaf > 50) { 
                my $ttest = new Statistics::TTest;  
                # $ttest->set_significance(90);
                $ttest->load_data($array1,\@array2);  
                my $s1=$ttest->{s1};
                my $s2=$ttest->{s2}; # sample 2  a Statistics::PointEstimation object
                my $mean1 = $s1->mean();
                my $mean2 = $s2->mean();
                my $nmean_diff = $ttest->mean_difference();
                my $value_ttest = $ttest->{t_prob};
                if (!$node_id){$node_id= "ROOT";} 
                print OUT5 "$node_id\t$value_ttest\t$fraction_above_in_leaf\t$count_above_cutoff\t$fraction_below_in_leaf\t$count_below_cutoff\n";
            }
        }
    }
}

sub log10 {
    my $n = shift;
    return log($n)/log(10);
}