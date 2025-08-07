use strict;
use warnings;
use Cwd;
use FindBin qw($Bin);
use File::Spec;
use Bio::SeqIO;
use Getopt::Long qw(GetOptions);
use lib "$Bin/..";
use utils qw(parse_enrichment_info
			get_control_enriched_reads
            make_dir);

exit main();

sub main {
    my $config = parse_args();

    my $DIR = getcwd;
    my $domain_dir = "domain_alignments";
    if ( ! -d $domain_dir ) { system("mkdir $domain_dir"); }

	my %domain_to_output;
    # Parse pfam file and enrichment info
    my $pfam_list = parse_pfam(
        $config->{pfam_hit},
        $config->{enrichment_txt},
        $config->{read_count_min},
        $config->{contig_min}
    );
    
    foreach my $domain (keys %$pfam_list) {
        my $instances = $pfam_list->{$domain}{instance};
        my $pfamname  = $pfam_list->{$domain}{name};
        my $tree_file = run_ks_test_tree($config, $domain);
        if ($instances > 0 && defined $tree_file && -e $tree_file) {
            my $generic = $domain . "_" . $pfamname;
			my $command1 = "python ".$config->{parse_tree}." ".$tree_file;
			my $output = `$command1`;
			if ($? != 0) {
				warn "Command failed: $command1\n";
				next;
			}
			$domain_to_output{$domain} = {
                pfamname => $pfamname,
                output   => $output
            };
        }
    }
	open my $fh, '>', 'PF4_tree_eval.txt' or die "Could not open file: $!";
    foreach my $domain (sort keys %domain_to_output) {
        my $pfamname = $domain_to_output{$domain}{pfamname};
        my $output   = $domain_to_output{$domain}{output};

        print $fh "Domain: $domain\n";
        print $fh "Pfam Name: $pfamname\n";
        print $fh "Output:\n$output\n";
        print $fh "-" x 40 . "\n";
    }
    return 0;
}

sub usage {
    print <<EOF;
Usage: $0 [options]

Required:
    --fasta FILE          Combined assembly FASTA file (AA) 
    --pfam_hit FILE       PFAM hit tab file with position of pfam
    --enrichment_txt FILE Enrichment info txt file
    --out                 output dir
Optional:
    --read_count_min INT Minimum total read number (default: 100)
    --contig_min INT     Minimum contig length (bp, default: 100)
    --cutoff FLOAT       Enrichment cutoff (default: 10)

Example:
    $0 --fasta combined.faa \\
       --pfam_hit pfamhit.tab \\
       --enrichment_txt enrichment_info.txt \\
       --read_count_min 100 \\
       --contig_min 100 \\
       --cutoff 10

EOF
    exit(1);
}

sub parse_args {
    my %config = (
        read_count_min => 100,
        contig_min     => 100,
        cutoff         => 10,
		parse_tree     => "$Bin/parse_tree.py",
        file_prefix    => ""
    );
    GetOptions(
        "fasta=s"         => \$config{fasta},
        "pfam_hit=s"      => \$config{pfam_hit},
        "enrichment_txt=s"=> \$config{enrichment_txt},
        "parse_tree=s"    => \$config{parse_tree},
        "read_count_min=i"=> \$config{read_count_min},
        "contig_min=i"    => \$config{contig_min},
		"out=s"			  => \$config{out},
        "file_prefix=s"   => \$config{file_prefix},
	    "cutoff=f"        => \$config{cutoff},
        "help|h"          => \$config{help},
    ) or usage();

    usage() if $config{help} || !$config{fasta} || !$config{pfam_hit} || !$config{enrichment_txt} || !$config{out};
    return \%config;
}

sub parse_pfam {
    my ($file, $enrichment_txt, $read_min, $length_min) = @_;
    my %result;
    my $enrichment_info = parse_enrichment_info($enrichment_txt);
    open(my $fh, '<', $file) or die "Can't open $file: $!\n";
    while (my $line = <$fh>) {
        chomp $line;
        my @tmp = split /\s+/, $line;
        next if @tmp < 5 || $tmp[0] !~ /^NODE/;

        my $contig = $tmp[0];
        my ($contig_length) = $contig =~ /_length_(\d+)_/;
        next unless defined $contig_length;
        $contig =~ s/[-_]\d+[FR]$//;
        my ($control_reads, $enriched_reads, $ratio) = get_control_enriched_reads($contig, $enrichment_info);
        my $read_count = $control_reads + $enriched_reads;

        my $pfam = $tmp[4];
        my $namepfam = $tmp[3];
        $namepfam =~ s/\s+/_/g;

        if ($read_count > $read_min && $contig_length > $length_min) {
            $result{$pfam}{"instance"}++;
            $result{$pfam}{"name"} = $namepfam;
        }
    }
    close $fh;
    return \%result;
}

sub run_ks_test_tree {
    my ($config, $pfamname) = @_;

    my $outdir = $config->{out};
    make_dir($outdir) or die "Failed to create directory $outdir: $!";
    my %id2seq;
    my $seq_in = Bio::SeqIO->new(-file => $config->{fasta}, -format => "fasta");
    while (my $faa = $seq_in->next_seq()) {
        my $id = $faa->id;
        $id2seq{$id} = $faa;
    }
    my %pfam2seqs;
    my $enrichment_info = parse_enrichment_info($config->{enrichment_txt});
    open(my $hit_fh, "<", $config->{pfam_hit}) or die "can't open $config->{pfam_hit}\n";
    while (my $line = <$hit_fh>) {
        chomp $line;
        $line =~ s/ +/&/g;
        my @tmp = split /\&/, $line;
        next if @tmp < 5 || $tmp[0] !~ /^NODE/;
        my $pfam = $tmp[4];
        next unless $pfam eq $pfamname;
        my $contig = $tmp[0];
        my $faa_id = $contig;
        $contig =~ s/((?:selection|control)).*$/$1/;
        my ($control_reads, $enriched_reads, $ratio) = get_control_enriched_reads($contig, $enrichment_info);
        my $enrichment = sprintf("%.1f", $ratio);
        $enrichment = 0.1 if $enrichment == 0;
        my $start = $tmp[17];
        my $end = $tmp[18];
        my $faa = $id2seq{$faa_id};
        next unless defined $faa;
        my $seq = $faa->subseq($start, $end);
        push @{$pfam2seqs{$pfam}}, {
            seq => $seq,
            enrichment => $enrichment
        };
    }
    close $hit_fh;
    return undef unless exists $pfam2seqs{$pfamname};

    my $pfam = $pfamname;
    my $out_prefix;
    if (defined $config->{file_prefix} && $config->{file_prefix} ne '') {
        $out_prefix = $config->{file_prefix} . "_$pfam";
    } else {
        $out_prefix = $pfam;
    }
    $out_prefix = File::Spec->catfile($config->{out}, $out_prefix);
    my $out_fasta   = "$out_prefix\_not_aligned.fasta";
    my $out_aligned = "$out_prefix\_aligned.fasta";
    my $out_tree    = "$out_prefix\_aligned.tree";

    open(my $fh, ">", $out_fasta) or die "can't open $out_fasta\n";
    my ($selected, $unselected) = (0, 0);
    foreach my $record (@{$pfam2seqs{$pfam}}) {
        my $seq = $record->{seq};
        my $enrichment = $record->{enrichment};
        my $id;
        if ($enrichment > $config->{cutoff}) {
            $selected++;
            $id = "S_" . $selected;
        } else {
            $unselected++;
            $id = "unselected_" . $unselected;
        }
        my $name = $id . "_" . $enrichment;
        print $fh ">$name\n$seq\n";
    }
    close $fh;

    my $cmd_align = "mafft --maxiterate 1000 --localpair $out_fasta > $out_aligned";
    my $cmd_tree  = "FastTree -gamma $out_aligned > $out_tree";
    system($cmd_align);
    system($cmd_tree);

    return $out_tree;
}