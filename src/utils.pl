package src::utils;

use strict;
use warnings;
use Exporter qw(import);

our @EXPORT_OK = qw(parse_pfam);

exit main();

sub main {
    # die "Usage: $0 dir1 dir2 ... output_dir\n" unless @ARGV >= 2;
    # my $out_dir = pop @ARGV;
    # my @dirs = @ARGV;
    my $out_dir = "appended_phage_test";
    my @dirs = ("phage_test", "phage_test_2", "phage_test_3");
    append_assemblies(@dirs, $out_dir);
}

sub append_assemblies {
    my $out_dir = pop @_;
    my @dirs = @_;
    my $base_path = "../output/";
    my $appended_fai;
    my $appended_fasta;
    my $appended_faa;
    foreach my $dir (@dirs) {
        my $assembly_path = $base_path.$dir."/assembly/";
        opendir(my $dh, $assembly_path) or die "Cannot open directory $assembly_path: $!";
        my @files = readdir($dh);
        closedir($dh);
        foreach my $file (@files) {
            my $full_path = $assembly_path . $file;
            if ($file =~ /control_and_selected\.fasta$/) {
                $appended_fasta .= slurp_file($full_path);
            }
            if ($file =~ /control_and_selected\.faa$/) {
                $appended_faa .= slurp_file($full_path);
            }
            if ($file =~ /control_and_selected\.fasta\.fai$/) {
                $appended_fai .= slurp_file($full_path);
            }
        }
    }
    system("mkdir -p ../output/".$out_dir."/assembly") unless -d "../output/".$out_dir."/assembly";

    write_file("../output/".$out_dir."/assembly/appended_control_and_selected.fasta", $appended_fasta);
    write_file("../output/".$out_dir."/assembly/appended_control_and_selected.faa",   $appended_faa);
    write_file("../output/".$out_dir."/assembly/appended_control_and_selected.fasta.fai", $appended_fai);
}

sub parse_pfam {
    my ($pfam)=@_;
    my %id2pfam;
    open (PFAM, $pfam) or die "can['t open $pfam\n";
    foreach my $line (<PFAM>) {
        chomp $line;
        my @tmp = split/\s+/, $line;
        my $id = $tmp[0];

        my $pfam = $tmp[4];
        $id2pfam{$id}=$pfam;
    }
    close PFAM;
    return (\%id2pfam);
}

sub slurp_file {
    my ($filename) = @_;
    local $/ = undef;
    open(my $fh, "<", $filename) or die "Can't open $filename: $!";
    my $content = <$fh>;
    close $fh;
    return $content;
}

sub write_file {
    my ($filename, $content) = @_;
    open(my $fh, ">", $filename) or die "Can't write to $filename: $!";
    print $fh $content;
    close $fh;
}