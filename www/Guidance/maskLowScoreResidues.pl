#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::Align::AlignI;
use Bio::AlignIO;
use warnings;

die "USAGE: $0 CODON_MSA_FILE GUIDANCE_RESIDUE_SCORES_FILE OUT_FILE CUTOFF ALPHABET
ALPHABET can be either aa or nuc\n" if (@ARGV < 5);

my ($msaFile,$scoreFile,$outFile,$cutoff,$alphabet) = @ARGV;

my $missingDataChar;
if ($alphabet eq "aa") {
	$missingDataChar = "X";
} elsif ($alphabet eq "nuc") {
	$missingDataChar = "N";
} else { die "ALPHABET must be either 'aa' or 'nuc'\n" }

my $in_fasta = Bio::SeqIO->new(-file => $msaFile, '-format' => 'fasta');
my @seqs;
my @ids;
while (my $seqObj = $in_fasta->next_seq()) {
	my @seq_chars = split(//,$seqObj->seq());
    push(@seqs,\@seq_chars);
    push(@ids,$seqObj->id());
}

#my $in_fasta  = Bio::AlignIO->new(-file => $msaFile , '-format' => 'fasta');
#my $aln = $in_fasta->next_aln();


open IN, "<$scoreFile" or die "can't open file $scoreFile";
#print "cutoff: $cutoff\n";
while (my $line = <IN>) {#COL_NUMBER	#ROW_NUMBER	#RES_PAIR_RESIDUE_SCORE
	chomp $line;
	next if ($line =~ m/^#/);
    if ($line =~ m/^\s*(\d+)\s+(\d+)\s+(\S+)$/) {
		if ($3 ne 'nan' and $3 < $cutoff) {
			my $col=$1-1;
			my $row=$2-1;
			$seqs[$row][$col] = $missingDataChar;
			#warn "DEBUG: masking $row,$col\n";
		}
	} else { warn "WARNING: failed to parse line: '$line'\n" }
}
close IN;

open OUT, ">$outFile";
for (my $i=0; $i<@seqs; ++$i) {
	my $id = $ids[$i];
    my $seqRef=$seqs[$i];
	my @seq_arr = @$seqRef;
	my $seq = join('',@seq_arr);
	print OUT ">$id\n";
	print OUT "$seq\n";
}
close OUT;
