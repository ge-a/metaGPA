use strict;
use warnings;
use Cwd;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);
use Bio::SeqUtils;
use Bio::SeqIO;
use Bio::Tools::Run::Hmmer;
use Bio::SearchIO;

exit main();

sub main {
    my ($faa, $out1, $out2, $out3, $pfam) = parse_args();

    run_hmmer($faa, $out1, $out2, $out3, $pfam);

    return 0;
}

sub parse_args {
    my $error_sentence = "USAGE : perl $0 --faa contig_file.faa --out1 output_file --out2 output_file_tab OPTIONAL : --pfam pfam_database (default : /mnt/home/ettwiller/laurence/database/pfam/Pfam-A.hmm) \n";
    my ($faa, $out1, $out2, $out3);
	#from : http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.gz
    my $pfam = "/mnt/home/ettwiller/laurence/database/pfam/Pfam-A.hmm"; 

    GetOptions(
        "faa=s"  => \$faa,    # translated assembly in all 6 frames
        "out1=s" => \$out1,   # output hmmer format.
        "out2=s" => \$out2,   # output tabulated format.
        "out3=s" => \$out3,   # output in stockholm format.
        "pfam=s" => \$pfam
    ) or die $error_sentence;

    die $error_sentence unless $faa && $out1 && $out2;
    return ($faa, $out1, $out2, $out3, $pfam);
}

sub run_hmmer {
    my ($faa, $out1, $out2, $out3, $pfam) = @_;

    my $command = "hmmsearch -o $out1 --domtblout $out2 $pfam $faa";
    system($command);

    if ($out3) {
        my $command1 = "hmmsearch -A $out3 $pfam $faa";
        system($command1);
    }
}