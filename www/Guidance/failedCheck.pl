#!/usr/local/bin/perl -w
use strict;

@ARGV >= 1 or die "USAGE: $0 guidanceDir1 [guidanceDir2...]\n";

my @guidanceDirs = @ARGV;

my ($c,@gr);
foreach my $dir (@guidanceDirs) {
	print "*** $dir\n";
	if (-e "$dir/Col_Scores_Graph.png") {
		print "Completed successfully\n";
	} elsif (-e "$dir/MSA.PRANK.aln.std") {
		$c = "cat $dir/MSA.PRANK.aln.std";
		@gr = `$c`;
		print "> $c\n@gr";
	}# elsif {
	
}
