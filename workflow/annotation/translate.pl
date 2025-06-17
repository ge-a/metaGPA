use strict;
use Cwd;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);
use Bio::SeqUtils;
use Bio::SeqIO;
use Bio::Tools::Run::Hmmer;
use Bio::SearchIO;


my $error_sentence = "USAGE : perl $0 --fasta contig_file.fasta --out contig_file_translated.faa OPTIONAL --min_size 100 (default 0)\n";

# declare the options upfront :
my $fasta;
my $out;
my $min_size =0; # minimum of size of the contig
#from : http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.gz
my $pfam = "/mnt/home/ettwiller/laurence/database/pfam/Pfam-A.hmm";

#get options :
GetOptions (    "fasta=s" => \$fasta,    # the Read1 file
		"min_size=s" => \$min_size,
		"out=s" => \$out #output file name.
    ) or die $error_sentence;

#=================================
#if something went wrong in getting the option, notify :
if (!$fasta || !$out) {die $error_sentence}
#=================================

my $DIR = getcwd;

my $seq_in = Bio::SeqIO->new( -format => 'fasta',
                              -file   => $fasta,
    );

open (OUT, ">$out") or die "can't save into $out\n";
while ( my $seq = $seq_in->next_seq() ) {
    my $id = $seq->id();
    my $size = $seq->length();
    # translate a sequence in all six frames
    #hmmscan --domtblout output_scan3 ../../../database/pfam/Pfam-A.hmm test.faa
    if ($size > $min_size) {
	    my @seqs = Bio::SeqUtils->translate_6frames($seq);

	    foreach my $prot (@seqs) {
            my $prot_id = $prot->id();
            my $prot_seq = $prot->seq();
            print OUT ">$prot_id\n$prot_seq\n";
	    }
    }
    #print "\n";
}
close OUT;

# run hmmsearch (similar for hmmpfam)
#my $factory = Bio::Tools::Run::Hmmer->new(-hmm => 'model.hmm');
# Pass the factory a Bio::Seq object or a file name, returns a Bio::SearchIO
#my $searchio = $factory->hmmsearch($seq);