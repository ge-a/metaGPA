use strict;
use Cwd;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);
use Bio::SeqUtils;
use Bio::SeqIO;
use Bio::Tools::Run::Hmmer;
use Bio::SearchIO;


my $error_sentence = "USAGE : perl $0 --faa contig_file.faa --out1 output_file --out2 output_file_tab OPTIONAL : --pfam pfam_database (default : /mnt/home/ettwiller/laurence/database/pfam/Pfam-A.hmm) \n";

# declare the options upfront :
my $faa;
my $out1;
my $out2;
my $out3;
#from : http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.gz
my $pfam = "/mnt/home/ettwiller/laurence/database/pfam/Pfam-A.hmm";

#get options :
GetOptions (    "faa=s" => \$faa,    # translated assembly in all 6 frames
		"out1=s" => \$out1, #output hmmer format.
		"out2=s" => \$out2, #output tabulated format.  
		"out3=s" => \$out3, #output in stockolm format.  
		"pfam=s" => \$pfam
    ) or die $error_sentence;

#=================================
#if something went wrong in getting the option, notify :
if (!$faa || !$out1 || !$out2) {die $error_sentence}
#=================================

my $DIR = getcwd;

my $command = "hmmsearch -o $out1 --domtblout $out2 $pfam $faa";
system($command);
if ($out3) {
    my $command1 = "hmmsearch -A $out3 $pfam $faa";
    system($command1);
}



