use strict;
use warnings;
use Cwd;
use Getopt::Long;

exit main();

sub main {
    my $config = parse_arguments();


    return 0;
}

sub parse_arguments {
    my %config = (
        commands => {
        },
    );
    
    GetOptions(

    ) or usage();

    usage() if $config{help};
    
    return \%config
}

sub usage {
    print <<EOF;
Usage: $0 [options]

Required:
    --control, -c FILE         Path to control FASTQ files (e.g., PF4.1_val_1.fq.gz)

Optional:
    --help, -h                Show this help message

Example:
    $0 --control PF4.1_val_1.fq.gz \\

EOF
    exit(1);
}