#!/usr/local/bin/perl -w
use strict;
use FileHandle;
use Bio::SeqIO;

####################################################################################################################
# Convert the names in a fasta file to numbers, and creates a code file with the names and the codes (running number)
###################################################################################################################


die "Usage: name2codefasta.pl <fastaFile> <codeFile> <outputCodedFile> " if (scalar(@ARGV) < 3);

my $in_fileName = $ARGV[0];
my $code_fileName = $ARGV[1];
my $out_fileName = $ARGV[2];

my $in_file = Bio::SeqIO->new(-file => $in_fileName , '-format' => 'Fasta');
my $code_file = new FileHandle(">$code_fileName") or die "Can't write to $code_fileName";
my $out_file = new FileHandle(">$out_fileName") or die "Can't write to $out_fileName";

my $counter = 1;
my $i;

while ( my $seqObj = $in_file->next_seq() ) {
    my $name = $seqObj->display_id();
    $name.= " ".$seqObj->desc()   if ($seqObj->desc());
    print $code_file "$name\t$counter\n";
    my $seq = $seqObj->seq();
    print $out_file ">$counter\n";
	for($i=0;$i<length($seq);$i+=60){
		print $out_file substr($seq,$i,60) . "\n";
	}
	if($i<length($seq)){
		print $out_file substr($seq,$i,length($seq)-$i);
	}
	print $out_file "\n";
    $counter++;
}
