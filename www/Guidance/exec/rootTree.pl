#!/usr/bin/perl -w

use strict;
use Bio::TreeIO;

die "USAGE: $0 inTree outTree
inTree must be an unrooted tree - i.e. root node has at least 3 sons.
In outTree the root will have 2 sons
(all direct sons of root will be made biforcating - the rest of tree is left untouched)\n"
	if (@ARGV < 2);
my ($inTree, $outTree) = @ARGV;

my $in  = new Bio::TreeIO(-file   => "$inTree",
						  -format => "newick");
my $out = new Bio::TreeIO(-file   => ">$outTree",
						  -format => "newick");

while (my $tree = $in->next_tree) {
	my $root = $tree->get_root_node;
	my @sons = $root->each_Descendent;

	# Remove edges between root-sons
	foreach my $son (@sons) {
		$root->remove_Descendent($son);
	}

	# Iteratively add 
	my $currFather = $root;
	foreach my $son (@sons) {
		$currFather->add_Descendent($son);
		my $midNode = new Bio::Tree::Node();
		$currFather->add_Descendent($midNode);
$midNode->branch_length(0);
		$currFather = $midNode;
	}
$currFather->ancestor->remove_Descendent($currFather);

	$out->write_tree($tree);
}
