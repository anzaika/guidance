#!/usr/local/bin/perl
use Bio::SeqIO;
use FileHandle;
use strict;

###############################################################################################################
# Takes a fasta MSA file created based on the script names2codeFasta.pl and reverts the codes to the names. 
###############################################################################################################

die "Usage: $0 <inFile> <codeFile> <outFile> " if (scalar(@ARGV) < 3);
my ($in_fileName,$code_fileName,$out_fileName) = @ARGV;


open CODE, "<$code_fileName" or die "Can't open $code_fileName";

my %codes2names;
while (my $line=<CODE>){
    $line =~ /(.+)\t(.+)/;
    my $code = $2;
    my $name = $1;
#    print "name: $name \t code: $code \n";
    $codes2names{$code}=$name;
    }

close CODE;
open IN, "<$in_fileName" or die "Can't open $in_fileName";
open OUT, ">$out_fileName" or die "Can't open $out_fileName";
my %name2seq;

foreach my $line (<IN>) {
    if ($line =~ m/^>(.*)/) {
	my $code = $1;
	my $new_name = $codes2names{$code};
#	print "$code is replaced with $new_name\n";
	print OUT ">$new_name\n";
    } else {
	print OUT "$line\n";
    }
}


close IN;
close OUT;


