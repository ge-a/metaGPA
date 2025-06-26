use strict;
use Cwd;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);

exit main();

sub main {
	my ($fq1, $fq2, $genome, $bam, $mapping_func) = parse_args();
	
	map_dna_to_fq($fq1, $fq2, $genome, $bam, $mapping_func);

	return 0;
}

sub parse_args {
	my $error_sentence = "USAGE : perl $0 --fq1 fastq1 --fq2 fastq2 --genome fasta_file --out bam\n";
	my ($fq1, $fq2, $genome, $bam, $mapping_func);
	$mapping_func = "bwa";

	GetOptions(
		"fq1=s" 		 => \$fq1,    # the Read1 file 
		"fq2=s" 		 => \$fq2,# the read 2 file
		"genome=s" 		 => \$genome, #the genome_file
		"out=s" 		 => \$bam, #output file name.
		"mapping-func=s" => \$mapping_func # mapping function, default is "bwa"
    ) or die $error_sentence;

	die $error_sentence unless $fq1 && $fq2 && $bam && $genome;

	return ($fq1, $fq2, $genome, $bam, $mapping_func);
}

sub map_dna_to_fq() {
	my ($fq1, $fq2, $genome, $bam, $mapping_func) = @_;

	my $generic = $bam; $generic =~ s/\.bam//;
	my $sam = $generic.".sam";
	my $bam = $generic.".bam";
	my $index = $genome.".bwt";

	my $command1;
	if ($mapping_func eq "bwa") {
		# This if statement never runs if I am not mistaken
		if (!$index) { 
			my $command0 = "bwa index $genome"; system($command0);
		}
		$command1 = "bwa mem $genome $fq1 $fq2  > $sam";
	}
	elsif ($mapping_func eq "bowtie") {
		my $command0 = "bowtie2-build -f $genome $genome > bowtie2-build.log";
		system($command0);
		$command1 = "bowtie2 -x $genome -1 $fq1 -2 $fq2 -S $sam";
	}
	elsif ($mapping_func eq "bbmap") {
		$command1 = "bbmap.sh ref=$genome in1=$fq1 in2=$fq2 minid=0.90 maxindel=3 outm=$sam";
	} else {
		die "Unknown mapping function: $mapping_func. Supported functions are bwa, bowtie, and bbmap.";
	}
	my $command2 = "samtools view -bS $sam | samtools sort -o $bam";
	my $command3 = "samtools index $bam";

	system($command1);
	system($command2);
	system($command3);
}