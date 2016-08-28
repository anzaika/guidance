#!/usr/bin/perl -w
use strict;
use warnings;

# pull out all the BS NJ trees into the BS directory
# pull out the original NJ tree (that was done on the complete MSA file)
die "Usage: $0 <noBSdir> <dataset> <num_BP_repeats> <alignment program>" if (scalar(@ARGV) < 4);

my $noBSdir = $ARGV[0];
my $dataset = $ARGV[1];
my $bp_repeats = $ARGV[2];

my $alnProg = $ARGV[3];

unless ($noBSdir =~ m/\/$/) {
    $noBSdir.="/";
}

my $BSdir = $noBSdir."BS/";
unless (-e  $BSdir) {
    system ("mkdir $BSdir");
}
# open semphy file
my $semphyLogFile = $BSdir.$dataset.".".$alnProg.".semphy.log";
print "semphy log file: $semphyLogFile\n";
open IN, "<$semphyLogFile" or die "can't open file $semphyLogFile";
my $countTrees=0;
my $readTree=0; # this is a flag. 2 == read the tree.
my $treeFile;
foreach my $line (<IN>) {
    if ($line =~ m/^\# Finished tree reconstruction\./) {
	$readTree=1;
    } elsif (($readTree==1) && ($line =~ m/^\# The tree/)) {
	$treeFile = $noBSdir.$dataset.".".$alnProg.".semphy.tree";
	$readTree=2;
    } 
    elsif (($readTree==1) && ($line =~ m/The reconsructed tree/)) {
	$readTree=2;
	my $treeDir = $BSdir."/tree_".$countTrees."/";
	unless (-e  $treeDir) {
	    system ("mkdir $treeDir");
	}
	$treeFile = $treeDir.$dataset.".".$alnProg.".semphy.tree_".$countTrees;
	$countTrees++;
    } elsif ($readTree==2) {
	# open treefile
	open OUT, ">$treeFile" or die "can't open file $treeFile";
	print "treeFile is: \n $treeFile \n";
	# write the tree into the treefile
	print OUT $line;
	# close treefile
	close OUT;
	$readTree=1;
    }
}
close IN;
if ($countTrees != $bp_repeats) {
    die "ERROR: dataset: $dataset \t countTrees: $countTrees while it should be $bp_repeats \n";
}





