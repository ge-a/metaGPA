use strict;
use Cwd;
use Getopt::Long;
use POSIX ":sys_wait_h";
use lib 'src';
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

    my @fastq_control1 = @{$config->{input}{control_1}};
    my @fastq_selection1 = @{$config->{input}{selection_1}};
    my @fastq_control2 = (defined $config->{input}{control_2}) ? @{$config->{input}{control_2}} : ();
    my @fastq_selection2 = (defined $config->{input}{selection_2}) ? @{$config->{input}{selection_2}} : ();

    # Join lists into comma-separated strings for passing to assembly
    my ($control_1, $control_2);
    if (@fastq_control2 && $fastq_control2[0] ne "") {
        ($control_1, $control_2) = process_fastq($fastq_control1[0], $fastq_control2[0], $prefix, $config->{dirs}{output}, $config->{commands}{do_trim});
    } else {
        ($control_1, $control_2) = process_fastq($fastq_control1[0], "", $prefix, $config->{dirs}{output}, $config->{commands}{do_trim});
    }

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
        $selection_2_str = join(",", @trimmed_selection1);
    }

    my $prefix = $fastq_control1[0];
    $prefix =~ s/.1_val_1.fq.gz//; $prefix =~ s/.*\///g;

    my $generic_control = $fastq_control1[0];
    my $generic_selection = $fastq_selection1[0];
    $generic_control =~ s/.1_val_1.fq.gz//; $generic_control =~ s/.*\///g;
    $generic_selection =~ s/.1_val_1.fq.gz//; $generic_selection =~ s/.*\///g;

    # Create commands with full lists
    my $assembly_cmd = "perl analysis/multitreatment/multi_assembly.pl".
        " --control-reads-1 ".$control_1.
        " --control-reads-2 ".$control_2.
        " --selection-reads-1 ".$selection_1_str.
        " --selection-reads-2 ".$selection_2_str.
        " --output-dir ".$config->{dirs}{output};

    my $mapping_cmd = "perl mapping/mapping.pl".
        " --control-reads-1 ".$control_1.
        " --control-reads-2 ".$control_2.
        " --selection-reads-1 ".$selection_1_str.
        " --selection-reads-2 ".$selection_2_str.
        " --prefix ".$prefix.
        " --mapping-func ".$config->{mapping_func}.
        " --output-dir ".$config->{dirs}{output};

    my $annotation_cmd = "perl annotation/annotation.pl".
        " --prefix ".$prefix.
        " --output-dir ".$config->{dirs}{output};

    my $enrichment_cmd = "perl enrichment/enrichment.pl".
        " --generic-control ".$generic_control.
        " --generic-selection ".$generic_selection.
        " --prefix ".$prefix.
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
        "control_1|c1=s@" => \$config{input}{control_1},
        "selection_1|s1=s@" => \$config{input}{selection_1},
        "control_2|c2=s@" => \$config{input}{control_2},
        "selection_2|s2=s@" => \$config{input}{selection_2},
        "outdir|o=s" => \$config{dirs}{output},
        "cutoff|c=f" => \$config{cutoff},
        "help|h" => \$config{help}
    ) or usage();

    $config{dirs}{output} = make_unique_path($config{dirs}{output}, $config{lite},
                                                $config{commands}{do_assembly}, 
                                                $config{commands}{do_annotation}, 
                                                $config{commands}{do_mapping}, 
                                                $config{commands}{do_enrichment});
    make_dir($config{dirs}{output});
    usage() if $config{help};
    usage() if (!defined($config{input}{control_1}) || 
                !defined($config{input}{selection_1}));
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
    --trim, -t <trimmer>      Specify trimmer to use (default: trim_galore)
    --assembly, -A            Run assembly steps
    --annotation, -AN         Run annotation steps
    --mapping, -M             Run mapping steps
    --enrichment, -E          Run enrichment analysis steps
    --cutoff -c FLOAT         Cutoff value for enrichment analysis (default: 3), if -1 then will be calculated from distribution
    --help, -h                Show this help message

Example:
    $0 --control PF4.1_val_1.fq.gz \\
       --selection PF1.1_val_1.fq.gz \\
       --outdir results \\
       --assembly --mapping --enrichment
EOF
    exit(1);
}