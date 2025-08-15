use strict;
use warnings;
use Cwd;
use Getopt::Long;
use FindBin qw($Bin);
use lib File::Spec->catdir($Bin, "src");
exit main();

sub main {
    my $config = parse_arguments();

    my $command = $config->{command};
    my $subhelp = $config->{subhelp};
    my @extra_args = @{$config->{extra_args}};

    my %script_map = (
        create_tree_file => File::Spec->catfile($Bin, "src", "analysis", "create_tree_file.pl"),
        get_enriched_fasta_orf => File::Spec->catfile($Bin, "src", "analysis", "get_enriched_fasta_orf.pl"),
        get_flanking_for_contig_vis => File::Spec->catfile($Bin, "src", "analysis", "get_flanking_for_contig_vis.pl"),
        get_gff => File::Spec->catfile($Bin, "src", "analysis", "get_gff.pl"),
        parse_tree => File::Spec->catfile($Bin, "src", "analysis", "parse_tree.pl"),
        residue_analysis => File::Spec->catfile($Bin, "src", "analysis", "residue_analysis.pl"),
        translate_orf => File::Spec->catfile($Bin, "src", "analysis", "translate_orf.pl"),
    );
    unless ($command && exists $script_map{$command}) {
        print STDERR "Unknown or missing command: $command\n";
        usage();
    }
    my $script_path = $script_map{$command};
    if ($subhelp) {
        my $cmd = "perl $script_path --help";
        print STDERR "Showing help for $command:\n";
        system($cmd) == 0 or die "Failed to run $cmd: $?";
        return 0;
    }
    my $cmd = "perl $script_path @extra_args";
    print STDERR "Running: $cmd\n";
    system($cmd) == 0 or die "Failed to run $cmd: $?";
    return 0;
}

sub parse_arguments {
    my $command;
    my @extra_args;
    Getopt::Long::Configure("pass_through");
    GetOptions(
        "command=s" => \$command,
        "subhelp|sh"   => \my $subhelp,
        "help|h"    => \my $help,
    ) or usage();

    @extra_args = @ARGV; # Pass all remaining arguments to the sub-script

    usage() if $help;

    return {
        command    => $command,
        subhelp    => $subhelp,
        extra_args => \@extra_args,
    };
}

sub usage {
    print <<EOF;
Usage: $0 --command <analysis_type> [options for analysis script]

Required:
    --command STRING      Analysis type (e.g., create_tree_file, get_flanking_for_contig_vis, get_gff)

Optional:
    --subhelp, -sh       Show help message for your specified command
    --help, -h           Show this help message

Example:
    $0 --command get_gff --control-reads-1 control_R1.fq.gz --output-dir results

Available commands:
    create_tree_file
    get_enriched_fasta_orf
    get_flanking_for_contig_vis
    get_gff
    parse_tree
    residue_analysis
EOF
    exit(1);
}