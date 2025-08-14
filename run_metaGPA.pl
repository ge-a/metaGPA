use strict;
use Cwd;
use Getopt::Long;
use POSIX ":sys_wait_h";
use FindBin qw($Bin);
use File::Spec;
use File::Basename;
use lib File::Spec->catdir($Bin, "src");
use utils qw(make_dir 
            process_fastq 
            make_unique_path);

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

    my $control_1 = $config->{input}{control_1};
    my $selection_1 = $config->{input}{selection_1};
    my $control2 = (defined $config->{input}{control_2}) ? $config->{input}{control_2} : "";
    my $selection_2 = (defined $config->{input}{selection_2}) ? $config->{input}{selection_2} : "";
    (my $prefix = basename($control_1)) =~ s/\.[^.]+(?:\.[^.]+)*$//;
        
    # pass control1/2 and selection1/2 to process_fastq if 2 is empty str then we generate the file ourselves
    $control_1, $control2 = process_fastq($control_1, $control2, $prefix, $config->{dirs}{output}, $config->{commands}{do_trim});
    $selection_1, $selection_2 = process_fastq($selection_1, $selection_2, $prefix, $config->{dirs}{output}, $config->{commands}{do_trim});
    # enforce in readme that if 2 file unspecified then must have same basename as 1 file with 1 replaced by 2
    # pass both files into assembly and mapping

    my $generic_control = $prefix;
    my $generic_selection = $selection_1;
    ($generic_selection = basename($generic_selection)) =~ s/\.[^.]+(?:\.[^.]+)*$//;

    my $assembly_path = File::Spec->catfile($Bin, "src", "assembly", "assembly.pl");
    my $mapping_path = File::Spec->catfile($Bin, "src", "mapping", "mapping.pl");
    my $annotation_path = File::Spec->catfile($Bin, "src", "annotation", "annotation.pl");
    my $enrichment_path = File::Spec->catfile($Bin, "src", "enrichment", "enrichment.pl");

    # Create commands
    my $assembly_cmd = "perl $assembly_path".
        " --control-reads-1 ".$control_1.
        " --control-reads-2 ".$control2.
        " --selection-reads-1 ".$selection_1.
        " --selection-reads-2 ".$selection_2.
        " --threads ".$config->{threads}.
        " --memory ".$config->{memory}.
        " --output-dir ".$config->{dirs}{output};
    my $mapping_cmd = "perl $mapping_path".
        " --control-reads-1 ".$control_1.
        " --control-reads-2 ".$control2.
        " --selection-reads-1 ".$selection_1.
        " --selection-reads-2 ".$selection_2.
        " --prefix ".$prefix.
        " --mapping-func ".$config->{mapping_func}.
        " --output-dir ".$config->{dirs}{output};
    my $annotation_cmd = "perl $annotation_path".
        " --prefix ".$prefix.
        " --output-dir ".$config->{dirs}{output};
    my $enrichment_cmd = "perl $enrichment_path".
        " --generic-control ".$generic_control.
        " --generic-selection ".$generic_selection.
        " --prefix ".$prefix.
        " --cutoff ".$config->{cutoff}.
        " --output-dir ".$config->{dirs}{output};

    # Run assembly
    my $assembly_pid = fork();
    if ($assembly_pid == 0) {
        # Child process for assembly
        exec($assembly_cmd) if $config->{commands}{do_assembly};
        exit(0);
    }

    # Wait for assembly to complete before starting annotation and mapping
    waitpid($assembly_pid, 0) if $config->{commands}{do_assembly};
    
    # Start mapping after assembly completes
    my $mapping_pid = fork();
    if ($mapping_pid == 0) {
        # Child process for mapping
        exec($mapping_cmd) if $config->{commands}{do_mapping};
        exit(0);
    }

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

    return 0;
}

sub parse_arguments {
    my %config = (
        commands => {
            do_trim => "none",
            do_assembly => "",
            do_annotation => "",
            do_mapping => "",
            do_enrichment => "",
        },
        lite => "none",
        mapping_func => "bwa", # add command line option to specify mapping function
        cutoff => 3,
        memory => 500,
        threads => 24
    );
    my $set_trim = sub { $config{commands}{do_trim} = "trim_galore" };

    # add option to map to run from a specific directory --> will create a new directory if one already exists
    GetOptions(
        "trim|t" => $set_trim,
        "assembly|A" => \$config{commands}{do_assembly},
        "annotation|AN" => \$config{commands}{do_annotation},
        "mapping|M" => \$config{commands}{do_mapping},
        "enrichment|E" => \$config{commands}{do_enrichment},
        "lite|l" => \$config{lite},
        "control_1|c1=s" => \$config{input}{control_1},
        "selection_1|s1=s" => \$config{input}{selection_1},
        "control_2|c2=s" => \$config{input}{control_2},
        "selection_2|s2=s" => \$config{input}{selection_2},
        "threads=i" => \$config{threads},
        "memory=i" => \$config{memory},
        "outdir|o=s" => \$config{dirs}{output},
        "cutoff|c=f" => \$config{cutoff},
        "help|h" => \$config{help}
    ) or usage();
    usage() if $config{help};

    my @steps = (
        $config{commands}{do_assembly},
        $config{commands}{do_annotation},
        $config{commands}{do_mapping},
        $config{commands}{do_enrichment}
    );
    my $num_enabled = grep { $_ && $_ ne "none" && $_ ne "" } @steps;
    if ($num_enabled > 0 && $num_enabled < 4) {
        $config{lite} = "lite";
    }

    $config{dirs}{output} = make_unique_path($config{dirs}{output}, $config{lite},
                                                $config{commands}{do_assembly}, 
                                                $config{commands}{do_annotation}, 
                                                $config{commands}{do_mapping}, 
                                                $config{commands}{do_enrichment});
    make_dir($config{dirs}{output});
    usage() if (!defined($config{input}{control_1}) || 
                !defined($config{input}{selection_1}));
    return \%config;
}

sub usage {
    print <<EOF;
Usage: $0 [options]

Required:
    --control_1, -c1 FILE     Path to control FASTQ files (e.g., PF4.1_val_1.fq.gz)
    --selection_1, -s1 FILE   Path to selection FASTQ files (e.g., PF1.1_val_1.fq.gz)
    --outdir, -o DIR          Output directory

Optional:
    --trim, -t <trimmer>      Specify trimmer to use (default: trim_galore)
    --assembly, -A            Run assembly steps
    --annotation, -AN         Run annotation steps
    --mapping, -M             Run mapping steps
    --enrichment, -E          Run enrichment analysis steps
    --control_2, -c2 FILE     Path to control FASTQ files opposite direction (e.g., PF4.2_val_2.fq.gz) (filename generated programatically if not specified)
    --selection_2, -s2 FILE   Path to selection FASTQ files opposite direction (e.g., PF1.2_val_2.fq.gz) (filename generated programatically if not specified)
    --threads  INTEGER        Number of threads to be ran on
    --memory   INTEGER        Amount of memory to be allocated
    --lite, l                 Overwrite an existing one if exists dir with same path instead of creating a new output dir with _{num} appended (default to create new dir)
    --cutoff -c FLOAT         Cutoff value for enrichment analysis (default: 3), if -1 then will be calculated from distribution
    --help, -h                Show this help message

Example:
    perl run_metaGPA.pl
    $0 --c1 PF4.1_val_1.fq.gz \\
       --s1 PF1.1_val_1.fq.gz \\
       --outdir results \\
       --assembly --mapping --enrichment
EOF
    exit(1);
}