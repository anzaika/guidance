#!/usr/local/bin/perl
use strict;

###############################################################################################################
# Works together (or rather after) the script names2codeFasta.pl. Takes a tree created based on 
# a fasta file with codes, and reverts the codes to the names. Required input is a code file which is created by
# names2codeFasta.pl
# ** very useful for working with all phyml and such, since these programs chop the name to 10 chars
###############################################################################################################


die "Usage: code2name.pl CODE_FILE TREE_FILE NEW_FILE NAME_LENGTH" if (scalar(@ARGV) < 3);
my $code2nameFile = $ARGV[0];
my $treeFile = $ARGV[1];
my $newFile = $ARGV[2];
my $nameLength = 30;
if (defined $ARGV[3]) {
	$nameLength = $ARGV[3];
}
  
	

my %names2code;
my @fields;


open FH, "<$code2nameFile";
#my $line;
while (my $line=<FH>){
	$line =~ /(.+)\t(\d+)/;
	my $code = $2;
	my $name = $1;
	$name =~ s/[\[\]\,\:\;\(\)]/_/g; #remove characters that are newick format associated
	if ($name =~ m/(.*\|.{$nameLength})/) {
		$name = $1;
	}
	$names2code{$code}=$name;
	print "$code $name\n";
}

close FH;

open TREE, "<$treeFile";
open NEWTREE, ">$newFile";

my $full_tree = "";
my $line2;
while ($line2 = <TREE>){ # this assumes there are bootstrap values on the input tree
	chomp $line2;
	$full_tree.=$line2;
	
}

@fields = split(/:/, $full_tree);

foreach my $field (@fields) {
	if ($field =~ /[\,\(](\d+)$/) { # a leaf comes either after a "(" or a ","
		$field =~ s/(\d+)$/$names2code{$1}/;
	}
	
	print NEWTREE "$field:";	

}

print NEWTREE "\n";	
