#!/usr/bin/perl -w
use strict;

use Bio::SeqIO;

die " $0 <alignment program (prank/mafft/clustal)> <dataset> <noBSdir> <original alignment program - default = same as current program> <num_BP_repeats - default 20> <extra suffix - default ''> required as input\n" if (scalar(@ARGV) < 3);

my $aln_program = "";
my $prog = $ARGV[0];

if ($prog eq "prank")  {
    $aln_program = "/groups/pupko/privmane/bin/prank";
} elsif ($prog eq "mafft"){
    $aln_program = "/groups/pupko/privmane/bin/mafft";
} elsif ($prog eq "clustal"){
    $aln_program = "/groups/pupko/privmane/alignment/run/clustalw";
}
else {
    die "please specify either mafft or prank as an alignment program\n";
}

my $dataset = $ARGV[1];
my $noBSdir = $ARGV[2]."/";

my $orig_prog = $prog;
if (defined $ARGV[3]) {
    $orig_prog = $ARGV[3];
}
my $bp_repeats = 20;
if (defined $ARGV[4]) {
    $bp_repeats = $ARGV[4];
}
my $suffix = "";
if (defined $ARGV[5]) {
    $suffix = $ARGV[5];
}

my $dir = $noBSdir."BS/";
my $aminoFile = $noBSdir.$dataset.".fas".$suffix;

for (my $countTrees=0;$countTrees<$bp_repeats;++$countTrees) {
    my $treeFile = $dir."tree_".$countTrees."/".$dataset.".".$orig_prog.".semphy.tree_".$countTrees.$suffix;
    if (-e $treeFile) {
	my $alnFile = $dir."tree_".$countTrees."/".$dataset.".".$orig_prog.".tree_".$countTrees.".".$prog.".aln";
	my $cmdFile = $dir."tree_".$countTrees."/".$dataset.".".$orig_prog.".tree_".$countTrees.".".$prog.".cmd";
	my $stdFile = $dir."tree_".$countTrees."/".$dataset.".".$orig_prog.".tree_".$countTrees.".".$prog.".std";
#	if (-e $alnFile) {
#	    print "skipping tree $countTrees because $alnFile already exists.\n";
#	    next;
#	}
	my $cmd = "";
	if ($prog eq "prank") {
	    $cmd = "$aln_program -F -d=$aminoFile -o=$alnFile -t=$treeFile >& $stdFile";
	} elsif ($prog eq "mafft") {
	    my $mafftFormatTreeFile = $treeFile.".mafftFormat";
	    unless (-e $mafftFormatTreeFile) {
		die "ERROR: file does not exist: $mafftFormatTreeFile\n";
	    }
	    
	    $cmd = "($aln_program --treein $mafftFormatTreeFile $aminoFile > $alnFile) >& $stdFile";
	} elsif ($prog eq "clustal") {
	    $cmd = "$aln_program -infile=$aminoFile -usetree=$treeFile -outfile=$alnFile >& $stdFile";
	}
	print "$cmdFile\n";
	open IN, ">$cmdFile" or die "cannot open file $cmdFile\n";
	print IN "#!/bin/sh\n\ncd ".$dir."tree_".$countTrees."\n\n";
	print IN "$cmd\n";
	close IN;
	system ($cmd);
	if ($prog eq "prank") {
	    system ("cp $alnFile.1.fas $alnFile");
	}
	#my $alias = $dataset.".".$prog;
	#my $qsub = "qsub -q heavy -N $alias $cmdFile";
	# print $qsub."\n";
	#system ($qsub);
    }
}
