#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::Align::AlignI;
use Bio::AlignIO;
use warnings;

# no checks of the input files are done.
die "USAGE: $0 MSA_FILE SP_FILE OUT_FILE CUTOFF" if (@ARGV < 4);

my ($msaFile,$spFile,$outFile,$cutoff) = @ARGV;

my $in_fasta = Bio::SeqIO->new(-file => $msaFile, '-format' => 'fasta');
my @seqs;

# read the file into an array
while (my $seqObj = $in_fasta->next_seq()) {
    push(@seqs,$seqObj);
}
open IN, "<$spFile" or die "can't open file $spFile";
open OUT, ">$outFile";

my $numRemovedPos = 0;
foreach my $line(<IN>) {
    if ($line =~ m/^\s*(\d+)\s+(\d+(\.\d+)?)/) {
	if ($2 < $cutoff) {
	    removePos($1);
	}
    }
}
foreach my $seqObj (@seqs) { 
    my $id = $seqObj->id();
    my $seq = $seqObj->seq();
    print OUT ">$id\n$seq\n";
}

sub removePos {
    print "removed pos $_[0]\n";
    my $pos2remove = $_[0] - $numRemovedPos;
    foreach my $seqObj (@seqs) { 
	my $new_seq = "";
	if ($pos2remove>1) {
	    $new_seq = $seqObj->subseq(1,$pos2remove-1);
	}
	if ($pos2remove< $seqObj->length()) {
	    $new_seq .= $seqObj->subseq($pos2remove+1,$seqObj->length());
	}
	$seqObj->seq($new_seq);
    } 
    $numRemovedPos++;
}
