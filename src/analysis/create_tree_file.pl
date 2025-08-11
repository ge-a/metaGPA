use strict;
use warnings;
use Cwd;
use FindBin qw($Bin);
use File::Spec;
use Bio::SeqIO;
use Parallel::ForkManager;
use File::Path qw(make_path);
use File::Basename;
use Time::HiRes qw(usleep);
use Getopt::Long qw(GetOptions);
use lib "$Bin/..";
use utils qw(parse_enrichment_info
			get_control_enriched_reads
            make_dir);

exit main();

sub main {
    my $config = parse_args();
    my $pfam_list = parse_pfam(
        $config->{pfam_hit},
        $config->{enrichment_txt},
        $config->{read_count_min},
        $config->{contig_min}
    );

    my $checkpoint_file = 'completed_domains.txt';
    my %completed = load_checkpoint($checkpoint_file);

    if ($config->{mode} eq 'fork') {
        run_fork_mode($config, $pfam_list, \%completed, $checkpoint_file);
    } elsif ($config->{mode} eq 'qsub') {
        run_qsub_mode($config, $pfam_list, \%completed);
    } elsif ($config->{mode} eq 'single') {
        run_single_mode($config, $pfam_list);
    } else {
        die "Invalid mode: $config->{mode}";
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
        mode           => "fork",
        parse_tree     => "$Bin/parse_tree.py",
        file_prefix    => "",
        domain         => undef,   # <--- add this
    );

    GetOptions(
        "fasta=s"          => \$config{fasta},
        "pfam_hit=s"       => \$config{pfam_hit},
        "enrichment_txt=s" => \$config{enrichment_txt},
        "parse_tree=s"     => \$config{parse_tree},
        "mode=s"           => \$config{mode},
        "read_count_min=i" => \$config{read_count_min},
        "contig_min=i"     => \$config{contig_min},
        "out=s"            => \$config{out},
        "file_prefix=s"    => \$config{file_prefix},
        "cutoff=f"         => \$config{cutoff},
        "domain=s"         => \$config{domain},   # <--- add this
        "help|h"           => \$config{help},
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

sub load_checkpoint {
    my ($file) = @_;
    my %completed;
    if (-e $file) {
        open my $fh, '<', $file;
        while (<$fh>) {
            chomp;
            $completed{$_} = 1;
        }
        close $fh;
    }
    return %completed;
}

sub run_fork_mode {
    my ($config, $pfam_list, $completed, $checkpoint_file) = @_;

    my $max_procs = 4;
    my %domain_to_output;
    open my $chk_out, '>>', $checkpoint_file or die "Can't write $checkpoint_file: $!";
    my $pm = Parallel::ForkManager->new($max_procs);

    $pm->run_on_finish(
        sub {
            my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data) = @_;
            return unless defined $data;

            # Merge returned hashref from child into parent hash
            my ($domain, $pfamname, $output) = @$data;
            $domain_to_output{$domain} = {
                pfamname => $pfamname,
                output   => $output
            };
        }
    );
    foreach my $domain (keys %$pfam_list) {
        next if $completed->{$domain};
        $pm->start and next;
        eval {
            local $SIG{ALRM} = sub { die "timeout\n" };
            alarm(600);  # 10 minute timeout

            my $instances = $pfam_list->{$domain}{instance};
            my $pfamname  = $pfam_list->{$domain}{name};
            my $tree_file = run_ks_test_tree($config, $domain);
            my $result;
            if ($instances > 20 && defined $tree_file && -e $tree_file) {
                my $command1 = "python $config->{parse_tree} $tree_file";
                my $output = `$command1`;
                if ($? == 0) {
                    print $chk_out "$domain\n";
                    $result = [$domain, $pfamname, $output];
                    $pm->finish(0, $result);
                } else {
                    warn "Command failed: $command1\n";
                    $pm->finish(0);
                }
            }
            alarm 0;
            exit;
        };
        warn "Domain $domain failed or timed out: $@" if $@;
        $pm->finish(0);
    }

    $pm->wait_all_children;
    close $chk_out;
    write_tree_summary('PF4_tree_eval.txt', \%domain_to_output);
}

sub run_qsub_mode {
    my ($config, $pfam_list, $completed) = @_;

    my $list_file = "pfam_domains.txt";
    open my $fh, '>', $list_file or die "Can't write $list_file: $!";
    my @domains;

    foreach my $domain (sort keys %$pfam_list) {
        next if $completed->{$domain};
        push @domains, $domain;
        print $fh "$domain\n";
    }
    close $fh;

    my $total = scalar @domains;
    my $wrapper_script = "run_tree_array.sh";

    open my $wrap, '>', $wrapper_script or die "Can't write $wrapper_script: $!";
    print $wrap generate_qsub_script($config, $total);
    close $wrap;

    print "\nJob array script written to $wrapper_script\n";
    print "Submit with: qsub $wrapper_script\n";
}

sub run_single_mode {
    my ($config, $pfam_list) = @_;
    my $domain = $config->{domain};
    die "Must provide --domain in single mode\n" unless defined $domain;
    die "Domain $domain not found or filtered out.\n" unless exists $pfam_list->{$domain};

    my $instances = $pfam_list->{$domain}{instance};
    my $pfamname  = $pfam_list->{$domain}{name};
    my $tree_file = run_ks_test_tree($config, $domain);

    if ($instances > 20 && defined $tree_file && -e $tree_file) {
        my $output = `python $config->{parse_tree} $tree_file`;
        if ($? == 0) {
            write_tree_summary('PF4_tree_eval.txt', {
                $domain => {
                    pfamname => $pfamname,
                    output   => $output
                }
            });
        } else {
            warn "Failed to parse tree for $domain\n";
        }
    } else {
        warn "No valid tree produced for $domain\n";
    }
}

sub generate_qsub_script {
    my ($config, $total) = @_;

    return qq{
#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 24
#$ -l ram=250G

DOMAIN=\$(sed -n "\${SGE_TASK_ID}p" pfam_domains.txt)

perl run_tree_builder.pl \\
    --domain \$DOMAIN \\
    --mode single \\
    --fasta "$config->{fasta}" \\
    --pfam_hit "$config->{pfam_hit}" \\
    --enrichment_txt "$config->{enrichment_txt}" \\
    --out "$config->{out}" \\
    --parse_tree "$config->{parse_tree}" \\
    --cutoff "$config->{cutoff}" \\
    --read_count_min "$config->{read_count_min}" \\
    --contig_min "$config->{contig_min}"
};
}

sub write_tree_summary {
    my ($file, $data) = @_;
    open my $fh, '>', $file or die "Could not open $file: $!";
    foreach my $domain (sort keys %$data) {
        my $pfamname = $data->{$domain}{pfamname};
        my $output   = $data->{$domain}{output};
        print $fh "Domain: $domain\nPfam Name: $pfamname\nOutput:\n$output\n";
        print $fh "-" x 40 . "\n";
    }
    close $fh;
}