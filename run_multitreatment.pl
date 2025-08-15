use strict;
use Cwd;
use POSIX ":sys_wait_h";
use FindBin qw($Bin);
use File::Basename;
use Getopt::Long;
use POSIX ":sys_wait_h";
use lib File::Spec->catdir($Bin, "src");
use utils qw(make_dir 
            process_fastq 
            make_unique_path);

exit main();

sub main {
    my $config = parse_arguments();

    # If no steps are specified then we run all of them
    if (!$config->{commands}{do_assembly} && 
        !$config->{commands}{do_annotation} && 
        !$config->{commands}{do_mapping} && 
        !$config->{commands}{do_enrichment}) {
        
        $config->{commands}{do_assembly} = "run";
        $config->{commands}{do_annotation} = "run";
        $config->{commands}{do_mapping} = "run";
        $config->{commands}{do_enrichment} = "run";
    }

    my $fastq_control1 = $config->{input}{control_1};
    my $fastq_control2 = (defined $config->{input}{control_2}) ? {$config->{input}{control_2}} : "";
    my @fastq_selection1 = @{$config->{input}{selection_1}};
    my @fastq_selection2 = (defined $config->{input}{selection_2}) ? @{$config->{input}{selection_2}} : ();
    (my $prefix = basename($fastq_control1)) =~ s/\.[^.]+(?:\.[^.]+)*$//; # generate general prefix 

    # Join lists into comma-separated strings for passing to assembly
    my ($control_1, $control_2);
    ($control_1, $control_2) = process_fastq($fastq_control1, $fastq_control2, $prefix, $config->{dirs}{output}, $config->{commands}{do_trim});

    my ($selection_1_str, $selection_2_str);
    if (@fastq_selection2 && $fastq_selection2[0] ne "") {
        # If selection2 files are provided, process and trim each pair
        my @trimmed_selection1;
        my @trimmed_selection2;
        for (my $i = 0; $i < @fastq_selection1; $i++) {
            my ($sel1, $sel2) = process_fastq($fastq_selection1[$i], $fastq_selection2[$i], $prefix."_sel".$i, $config->{dirs}{output}, $config->{commands}{do_trim});
            push @trimmed_selection1, $sel1;
            push @trimmed_selection2, $sel2;
        }
        $selection_1_str = join(",", @trimmed_selection1);
        $selection_2_str = join(",", @trimmed_selection2);
    } else {
        # Only selection1 files, process and trim each
        my @trimmed_selection1;
        my @trimmed_selection2;
        for (my $i = 0; $i < @fastq_selection1; $i++) {
            my ($sel1, $sel2) = process_fastq($fastq_selection1[$i], "", $prefix."_sel".$i, $config->{dirs}{output}, $config->{commands}{do_trim});
            push @trimmed_selection1, $sel1;
            push @trimmed_selection2, $sel2;
        }
        $selection_1_str = join(",", @trimmed_selection1);
        $selection_2_str = join(",", @trimmed_selection2);
    }

    # Generic and Selection prefix generation
    my $generic_control = $prefix;
    my @selection_files = split /,/, $selection_1_str;
    my @generic_selections = map {
        my $name = basename($_); # Strip path
        $name =~ s/\.[^.]+(?:\.[^.]+)*$//; # Remove all dot extensions
        $name;
    } @selection_files;
    my $selection_list = join(",", @generic_selections);

    my $assembly_path = File::Spec->catfile($Bin, "src", "analysis", "multitreatment", "multi_assembly.pl");
    my $mapping_path = File::Spec->catfile($Bin, "src", "analysis", "multitreatment", "multi_mapping.pl");
    my $annotation_path = File::Spec->catfile($Bin, "src", "annotation", "annotation.pl");
    my $enrichment_path = File::Spec->catfile($Bin, "src", "analysis", "multitreatment", "multi_enrichment.pl");

    my $num_selections = get_num_selection_reads($selection_1_str);
    # Create commands with full lists
    my $assembly_cmd = "perl $assembly_path".
        " --control-reads-1 ".$control_1.
        " --control-reads-2 ".$control_2.
        " --selection-reads-1 ".$selection_1_str.
        " --selection-reads-2 ".$selection_2_str.
        " --threads ".$config->{threads}.
        " --memory ".$config->{memory}.
        " --num-selections ".$num_selections.
        " --output-dir ".$config->{dirs}{output};

    my $mapping_cmd = "perl $mapping_path".
        " --control-reads-1 ".$control_1.
        " --control-reads-2 ".$control_2.
        " --selection-reads-1 ".$selection_1_str.
        " --selection-reads-2 ".$selection_2_str.
        " --num-selections ".$num_selections.
        " --prefix ".$prefix.
        " --mapping-func ".$config->{mapping_func}.
        " --output-dir ".$config->{dirs}{output};

    my $annotation_cmd = "perl $annotation_path".
        " --prefix ".$prefix.
        " --output-dir ".$config->{dirs}{output};

    my $enrichment_cmd = "perl $enrichment_path".
        " --generic-control ".$generic_control.
        " --generic-selection ".$selection_list.
        " --prefix ".$prefix.
        " --num-selections ".$num_selections.
        " --cutoff ".$config->{cutoff}.
        " --output-dir ".$config->{dirs}{output};

    # Run assembly
    my $assembly_pid = fork();
    if ($assembly_pid == 0) {
        exec($assembly_cmd) if $config->{commands}{do_assembly};
        exit(0);
    }
    waitpid($assembly_pid, 0) if $config->{commands}{do_assembly};

    # Start mapping after assembly completes
    my $mapping_pid = fork();
    if ($mapping_pid == 0) {
        exec($mapping_cmd) if $config->{commands}{do_mapping};
        exit(0);
    }

    # Start annotation after assembly completes
    my $annotation_pid = fork();
    if ($annotation_pid == 0) {
        exec($annotation_cmd) if $config->{commands}{do_annotation};
        exit(0);
    }

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
        "control_1|c1=s" => \$config{input}{control_1},
        "selection_1|s1=s@" => \$config{input}{selection_1},
        "control_2|c2=s" => \$config{input}{control_2},
        "selection_2|s2=s@" => \$config{input}{selection_2},
        "threads=i" => \$config{threads},
        "memory=i" => \$config{memory},
        "lite|l"  => \$config{lite},
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
    --control_1, -c1 FILE      Path to control FASTQ files (e.g., PF4.1_val_1.fq.gz)
    --selection_1, -s1 FILE    Path to selection FASTQ files (e.g., PF1.1_val_1.fq.gz)
    --outdir, -o DIR           Output directory

Optional:
    --trim, -t <trimmer>      Specify trimmer to use (default: trim_galore)
    --assembly, -A            Run assembly steps
    --annotation, -AN         Run annotation steps
    --mapping, -M             Run mapping steps
    --enrichment, -E          Run enrichment analysis steps
    --control_2, -c2 FILE     Path to control FASTQ files opposite direction (e.g., PF4.2_val_2.fq.gz) (filename generated programatically if not specified)
    --selection_2, -s2 FILES  Path to selection FASTQ files opposite direction (e.g., PF1.2_val_2.fq.gz) (filename generated programatically if not specified)
    --threads  INTEGER        Number of threads to be ran on
    --memory   INTEGER        Amount of memory to be allocated
    --lite, l                 Overwrite an existing one if exists dir with same path instead of creating a new output dir with _{num} appended (default to create new dir)
    --cutoff -c FLOAT         Cutoff value for enrichment analysis (default: 3), if -1 then will be calculated from distribution
    --help, -h                Show this help message

Example:
    $0 --control_1 PF4.1_val_1.fq.gz \\
       --selection_1 PF1.1_val_1.fq.gz \\
       --selection_1 PF2.1_val_1.fq.gz \\
       --selection_1 PF3.1_val_1.fq.gz \\
       --outdir results \\
EOF
    exit(1);
}

sub get_num_selection_reads {
    my ($selection_reads_1_str) = @_;
    my @reads = split /,/, $selection_reads_1_str;
    my %unique;
    $unique{$_}++ for @reads;
    return scalar(keys %unique);
}