use strict;
use Cwd;
use Getopt::Long;
use POSIX ":sys_wait_h";

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

        # Create commands
        my $assembly_cmd = "perl assembly/assembly.pl".
            " --control-reads-1 ".$control_1.
            " --selection-reads-1 ".$selection_1.
            " --output-dir ".$config->{dirs}{output};
        my $mapping_cmd = "perl mapping/mapping.pl".
            " --control-reads-1 ".$control_1.
            " --selection-reads-1 ".$selection_1.
            " --output-dir ".$config->{dirs}{output};
        my $annotation_cmd = "perl annotation/annotation.pl".
            " --prefix ".$prefix.
            " --output-dir ".$config->{dirs}{output};
        my $enrichment_cmd = "perl enrichment/enrichment.pl".
            " --generic-control ".$generic_control.
            " --generic-selection ".$generic_selection.
            " --prefix ".$prefix.
            " --output-dir ".$config->{dirs}{output};

        # Run assembly and mapping concurrently
        my $assembly_pid = fork();
        if ($assembly_pid == 0) {
            # Child process for assembly
            exec($assembly_cmd) if $config->{commands}{do_assembly};
            exit(0);
        }

        my $mapping_pid = fork();
        if ($mapping_pid == 0) {
            # Child process for mapping
            exec($mapping_cmd) if $config->{commands}{do_mapping};
            exit(0);
        }

        # Wait for assembly to complete before starting annotation
        waitpid($assembly_pid, 0) if $config->{commands}{do_assembly};
        
        # Start annotation after assembly completes
        my $annotation_pid = fork();
        if ($annotation_pid == 0) {
            # Child process for annotation
            exec($annotation_cmd) if $config->{commands}{do_annotation};
            exit(0);
        }

        # Wait for both mapping and annotation to complete
        waitpid($mapping_pid, 0) if $config->{commands}{do_mapping};
        waitpid($annotation_pid, 0) if $config->{commands}{do_annotation};

        # Run enrichment after both mapping and annotation complete
        if ($config->{commands}{do_enrichment}) {
            system($enrichment_cmd) == 0 
                or die "Failed to execute command: $enrichment_cmd";
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

    $config{dirs}{output} = make_unique_path($config{dirs}{output});
    system("mkdir -p $config{dirs}{output}") unless -d $config{dirs}{output};
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

    my $output_dir = "../output";
    system("mkdir -p $output_dir") unless -d $output_dir;

    my $basename = $path;
    $basename =~ s/.*\///g; 
    my $full_path = "$output_dir/$basename";

    my $counter = 1;
    my $unique_path = $full_path;
    while (-e $unique_path) {
        $unique_path = "${path}_$counter";
        $counter++;
    }

    return $unique_path;
}
