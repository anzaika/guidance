#!/usr/bin/perl -w
use strict;


###############################################################
# convertNacessFormat.pl	
# Usage: convertNacessFormat.pl <in file> <outFile in surfracer format>
# Description: change the output format of naccess to match that of surfracer
###############################################################


if (@ARGV < 2){
    die "Usage: cconvertNacessFormat.pl <in file> <outFile in surfracer format> "
}

my $inFileName = $ARGV[0];
my $outFileName = $ARGV[1];
my $line;

open INFILE, "< $inFileName" or die "unable to open input file $!";
open OUTFILE, ">  $outFileName" or die "unable to open file $outFileName!";

while ($line = <INFILE>){
    chomp($line);
    if($line =~ m/^ATOM/){
	my @lineArray = split(/\s+/,$line);
	my $line_begin = substr($line,0,54);
	my $arraySize = @lineArray;
	my $asa = $lineArray[$arraySize-2];
	$asa = substr($asa, 0, length($asa)-1); #nacess produce 3 digits after dot while surfracer only 2
	my $radius = $lineArray[$arraySize-1]; #radius is 4 chars
	my $charsToAdd = 6 - length($asa);
	my $whites = "";
	for (my $i=0; $i < $charsToAdd; ++$i) {
	    $whites .= " ";
	}
	my $newLine = $line_begin ."  " . $radius . "   " . $asa . $whites ."0.00  " . "0.00      ". "\n";
	print OUTFILE $newLine;
    }
    else{
	print OUTFILE "$line\n";
    }
}
close(INFILE);			
close(OUTFILE)
