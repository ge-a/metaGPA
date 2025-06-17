use strict;
use Cwd;
use Getopt::Long;

exit main();

sub main {
    my $config = parse_arguments();

    if (!$config->{commands}{do_assembly} && 
        !$config->{commands}{do_annotation} && 
        !$config->{commands}{do_mapping} && 
        !$config->{commands}{do_enrichment}) {
        
        $config->{commands}{do_assembly} = "run";
        $config->{commands}{do_annotation} = "run";
        $config->{commands}{do_mapping} = "run";
        $config->{commands}{do_enrichment} = "run";
    }

    my @fastq_control = @{$config->{input}{control_fq}};
    my @fastq_selection = @{$config->{input}{selection_fq}};
    my $size = @fastq_control;
 
    for (my $i=0; $i<$size; $i++) {
        my @commands = ();
        my $control_1 = $fastq_control[$i];
        my $selection_1 = $fastq_selection[$i];
        
        my $generic_control = $control_1;
        my $generic_selection = $selection_1;
        $generic_control =~ s/.1_val_1.fq.gz//; $generic_control =~ s/.*\///g;
        $generic_selection =~ s/.1_val_1.fq.gz//; $generic_selection =~ s/.*\///g;  

        my $prefix = $control_1;
        $prefix =~ s/.1_val_1.fq.gz//; $prefix =~ s/.*\///g;
        $prefix =~ s/control/experiment/;

        push @commands, "perl assembly/assembly.pl".
            " --control-reads-1 ".$control_1.
            " --selection-reads-1 ".$selection_1.
            " --output-dir ".$config->{dirs}{output} if $config->{commands}{do_assembly};
        push @commands, "perl annotation/annotation.pl".
            " --prefix ".$prefix.
            " --output-dir ".$config->{dirs}{output} if $config->{commands}{do_annotation};
        push @commands, "perl mapping/mapping.pl".
            " --control-reads-1 ".$control_1.
            " --selection-reads-1 ".$selection_1.
            " --output-dir ".$config->{dirs}{output} if $config->{commands}{do_mapping};
        push @commands, "perl enrichment/enrichment.pl".
            " --generic-control ".$generic_control.
            " --generic-selection ".$generic_selection.
            " --prefix ".$prefix.
            " --output-dir ".$config->{dirs}{output} if $config->{commands}{do_enrichment};
        foreach my $command (@commands) {
            print "Running command: $command\n";
            system($command) == 0 or die "Failed to execute command: $command";
        }
    }

    return 0;
}

sub parse_arguments {
    my %config = (
        commands => {
            do_assembly => "",
            do_annotation => "",
            do_mapping => "",
            do_enrichment => ""
        }
    );

    GetOptions(
        "assembly|A" => \$config{commands}{do_assembly},
        "annotation|AN" => \$config{commands}{do_annotation},
        "mapping|M" => \$config{commands}{do_mapping},
        "enrichment|E" => \$config{commands}{do_enrichment},
        "control|c=s@" => \$config{input}{control_fq},
        "selection|s=s@" => \$config{input}{selection_fq},
        "outdir|o=s" => \$config{dirs}{output},
        "help|h" => \$config{help}
    ) or usage();

    make_unique_path($config{dirs}{output}) unless -d $config{dirs}{output};
    usage() if $config{help};
    usage() if (!defined($config{input}{control_fq}) || 
                !defined($config{input}{selection_fq}));
    return \%config;
}

sub usage {
    print <<EOF;
Usage: $0 [options]

Required:
    --control, -c FILE         Path to control FASTQ files (e.g., PF4.1_val_1.fq.gz)
    --selection, -s FILE       Path to selection FASTQ files (e.g., PF1.1_val_1.fq.gz)
    --outdir, -o DIR          Output directory

Optional:
    --assembly, -A            Run assembly steps
    --annotation, -AN         Run annotation steps
    --mapping, -M             Run mapping steps
    --enrichment, -E          Run enrichment analysis steps
    --help, -h               Show this help message

Example:
    $0 --control PF4.1_val_1.fq.gz \\
       --selection PF1.1_val_1.fq.gz \\
       --outdir results \\
       --assembly --mapping --enrichment
EOF
    exit(1);
}

sub make_unique_path {
    my ($path) = @_;
    my $counter = 1;
    my $unique_path = $path;
    while (-e $unique_path) {
        $unique_path = "${path}_$counter";
        $counter++;
    }
    system("mkdir -p $unique_path");
    return $unique_path;
}
