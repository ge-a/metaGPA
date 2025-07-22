use strict;
use warnings;
use Cwd;
use Getopt::Long;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);
use Bio::SeqIO;
use utils qw(make_dir);

exit main();

sub main {
    my $config = parse_arguments();
    my @commands = write_commands($config);
    
    foreach my $command (@commands) {
        print "Running command: $command\n";
        system($command) == 0 or die "Failed to execute command: $command";
    }
    return 0;
}

sub parse_arguments {
    my %config = (
        params => {
            min_length => 500,
        },
        dirs => {
            output => "output",
            annotation => "annotation",
        },
        programs => {
            translate => "annotation/translate.pl",
            hmmer => "annotation/run_hmmer.pl",
        }
    );
    GetOptions(
        "prefix=s" => \$config{input}{prefix},
        "output-dir=s" => \$config{dirs}{output},
        "min-length=i" => \$config{params}{min_length},
        "help|h" => \$config{help}
    ) or usage();

    $config{dirs}{annotation} = $config{dirs}{output}."/annotation";
    make_dir($config{dirs}{output});
    make_dir($config{dirs}{annotation});

    $config{dirs}{hmmer1} = $config{dirs}{annotation}."/".$config{input}{prefix}."control_and_selected_hmmer_format.hmmer"; 
    $config{dirs}{hmmer2} = $config{dirs}{annotation}."/".$config{input}{prefix}."control_and_selected_hmmer_format.tab";
    $config{dirs}{hmmer3} = $config{dirs}{annotation}."/".$config{input}{prefix}."control_and_selected_hmmer_format_TIGRFAM.hmmer";
    $config{dirs}{hmmer4} = $config{dirs}{annotation}."/".$config{input}{prefix}."control_and_selected_hmmer_format_TIGRFAM.tab";

    usage() if $config{help};

    return \%config;
}

sub usage {
    print <<EOF;
Usage: $0 [options]

Required:
    --prefix PREFIX                  Prefix for output files
    --output-dir DIR                 Output directory path

Optional:
    --min-length INT             Minimum sequence length (default: 500)
    --help                       Show this help message

Example:
    $0 --assembly_final_info assembly/experiment_control_and_selected.fasta \\
       --translated_assembly_final assembly/experiment_control_and_selected.faa \\
       --output-dir results
EOF
    exit(1);
}


sub write_commands {
    my ($config) = @_;
    my @commands;

    my $assembly_file = $config->{dirs}{output}."/assembly/".$config->{input}{prefix}."control_and_selected.nd.fasta";
    my $translated_assembly_final = $config->{dirs}{output}."/assembly/".$config->{input}{prefix}."control_and_selected.faa";

    push @commands, "perl ".$config->{programs}{translate}.
        " --fasta ".$assembly_file.
        " --out ".$translated_assembly_final. 
        " --min_size ".$config->{params}{min_length};
    push @commands, "perl ".$config->{programs}{hmmer}.
        " --faa ".$translated_assembly_final. 
        " --out1 ".$config->{dirs}{hmmer1}.
        " --out2 ".$config->{dirs}{hmmer2};
    push @commands, "perl ".$config->{programs}{hmmer}. 
        " --faa ".$translated_assembly_final.
        " --out1 ".$config->{dirs}{hmmer3}.
        " --out2 ".$config->{dirs}{hmmer4}. 
        " --pfam /mnt/home/ettwiller/laurence/database/TIGRFAM/TIGRFAM_2021_11_18.hmm";

    return @commands;
}