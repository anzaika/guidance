#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

#convert the tree to mafft format

die " $0 <BSdir> <dataset> <original alignment program> <num_BP_repeats> <extra suffix - default ''> required as input\n" if (scalar(@ARGV) < 4);

my $bsDir = $ARGV[0];
my $dataset = $ARGV[1];
my $orig_prog = $ARGV[2];
my $bp_repeats = $ARGV[3];
my $suffix = "";
if (defined $ARGV[4]) {
    $suffix = $ARGV[4];
}

$bsDir .= "/";

for (my $countTrees=0;$countTrees<$bp_repeats;++$countTrees) {
    print "**** $dataset tree_$countTrees\n";
    my $treeFile = $bsDir."tree_".$countTrees."/".$dataset.".".$orig_prog.".semphy.tree_".$countTrees.$suffix;
    if (-e $treeFile) {
	my $rootedTreeFile = $treeFile.".rooted";
	my $mafftFormatTreeFile = $treeFile.".mafftFormat";
	print "treeFile: $treeFile\nrootedTreeFile: $rootedTreeFile\nmafftFormatTreeFile: $mafftFormatTreeFile\n";

	system("/groups/pupko/privmane/pupkoSVN/trunk/scripts/rootTree.pl $treeFile $rootedTreeFile");
	open IN, "<$rootedTreeFile", or die "can't open $rootedTreeFile";
	my $newick = "";
	foreach my $line (<IN>) { 
	    $newick.=$line;
	}
	close IN;

#	if ($newick =~ m/:0[^\.]/ && $newick !~ m/:0\./) {
	if ($newick =~ m/:-/) {
	    my $rootedTreeFile_withMinusLengths = $rootedTreeFile.".withMinusLengths";
	    my $command = "mv $rootedTreeFile $rootedTreeFile_withMinusLengths";
	    system ($command);
	    $newick =~ s/:-/:/g;
	    $command = "echo '$newick' > $rootedTreeFile";
	    system ($command);
	}
	if ($newick =~ m/:\d+[^\.]/) {
	    my $rootedTreeFile_badBranchLength = $rootedTreeFile.".badBranchLength";
	    my $command = "mv $rootedTreeFile $rootedTreeFile_badBranchLength";
	    system ($command);
	    $newick =~ s/:(\d+)([^\.])/:$1\.0$2/g;
#	    print "$newick\n";
	    $command = "echo '$newick' > $rootedTreeFile";
	    system ($command);
	}

	system("/groups/pupko/privmane/programs/mafft/newick2mafft.rb $rootedTreeFile > $mafftFormatTreeFile\n"); 
    } else {
	print "File does not exist: $treeFile\n";
    }
}
