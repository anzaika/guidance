#!/usr/bin/perl -w

use strict;
use Bio::AlignIO;

@ARGV == 3 or die "USAGE: $0 IN_MSA_FILE OUT_HTML_FILE SCORES_FILE
SCORES_FILE - Each line should contain three values: column number, seq number, and a score between 0 and 1 (separated by white spaces)\n";

my($inMsaFile, $outHtmlFile, $scoresFile) = @ARGV;

# Read scores
open SCORES, $scoresFile or die "Can't open $scoresFile: $!";
my %scores;
foreach (<SCORES>) {
	next if (/^#/);
	s/^\s+//;
	my ($col, $seq, $score) = split;
	$scores{$seq}[$col] = $score;
}

# Read MSA
my $in  = Bio::AlignIO->new( '-format' => 'fasta' , -file => $inMsaFile) or die "Can't open $inMsaFile: $!";
my $aln = $in->next_aln;
$aln->verbose(1);
# Otherwise, bioperl adds sequence start/stop values
$aln->set_displayname_flat();


# Print HTML start
# Code from Conseq colored MSA: ~/pupkoSVN/trunk/www/conseq/runCalc_Conseq.pl line 985
my %msaColors = ();
my %msaPrintColors = ();
my $lineCounter;
my @line;
my $key;

my $fontSize=2;
my $sequenceLengthForDisplay=400000;
my @msaRightOrder=0;
my $msaRightOrderCounter=0;
my $tdWidth = 5;

my @colorstep = (); #color steps
$colorstep[0] = "#10C8D1"; #Not confident
$colorstep[1] = "#8CFFFF";
$colorstep[2] = "#D7FFFF";
$colorstep[3] = "#EAFFFF";
$colorstep[4] = "#FFFFFF"; #average
$colorstep[5] = "#FCEDF4";
$colorstep[6] = "#FAC9DE";
$colorstep[7] = "#F07DAB";
$colorstep[8] = "#A02560"; #Most confident

open MSACOLOREDHTML, ">$outHtmlFile" or die "Can't open $outHtmlFile: $!";
print MSACOLOREDHTML "<html>\n<head>\n</head>\n<body>\n\n";
print MSACOLOREDHTML "<H1 align=center><u>MSA color-coded by ABP scores</u></H1>\n\n";
print MSACOLOREDHTML "<table border=0  CELLSPACING=1  CELLPADDING=0 >\n";


# Print colored HTML

# counts how many times we print the whole section (relevants to sequences longer than the sequenceLengthForDisplay) 
for(my $blockStart=1; $blockStart<$aln->length; $blockStart+=$sequenceLengthForDisplay) {
	my $blockEnd = $blockStart+$sequenceLengthForDisplay;
	$blockEnd = $aln->length if ($blockEnd > $aln->length);

	# Iterate over sequences and print up to sequenceLengthForDisplay residues
	foreach my $seq ($aln->each_seq) {

#		next if (   $seq->id < 42
#				 || $seq->id > 73);

		# Print seq id
		print MSACOLOREDHTML "<tr>\n";
		print MSACOLOREDHTML "<td><b><font face='Courier New' color='black' size=$fontSize>", $seq->id, "</font></b></td>\n";

		# Print seq
		my @seq = split //, $seq->subseq($blockStart, $blockEnd);

		for(my $pos=0; $pos<@seq; $pos++) {
#		for(my $pos=757; $pos<875; $pos++) {
			my $res = $seq[$pos];
			if ($res eq '-') {
				print MSACOLOREDHTML "<td width=$tdWidth><b><font face='Courier New' color='black' size=$fontSize>$res</font></b></td>\n";
			} else {
				my $color = $colorstep[ int(9 * $scores{$seq->id}[$pos+1]) ];
				if($color eq "#A02560"){
					print MSACOLOREDHTML "<td width=$tdWidth><b><font face='Courier New' color='white' size=$fontSize><span style='background: $color;'>$res</span></font></b></td>\n";
				}
				else {                
					print MSACOLOREDHTML "<td width=$tdWidth><b><font face='Courier New' color='black' size=$fontSize><span style='background: $color;'>$res</span></font></b></td>\n";
				}
			}
		}
		print MSACOLOREDHTML "</tr>\n\n";
	}

	print MSACOLOREDHTML "<tr>&nbsp</tr>\n\n";
}
print MSACOLOREDHTML "</table>";

# print the color scale
print  MSACOLOREDHTML "\n<br><b><u>Legend:</u><br><br>\nThe alignment confidence scale:</b><br>\n<table border=0 cols=1 width=310>\n<tr><td align=center>\n<font face='Courier New' color='black' size=+1><center>\n";
for (my $i=8 ; $i>=0 ; $i--){
    
	if ($i == 0){
		
		print  MSACOLOREDHTML "<font face='Courier New' color='white' size=$fontSize><span style='background: $colorstep[$i];'>&nbsp;", $i+1, "&nbsp;</span></font>";
	}
	else {
		
		print  MSACOLOREDHTML "<font face='Courier New' color='black' size=$fontSize><span style='background: $colorstep[$i];'>&nbsp;", $i+1, "&nbsp;</span></font>";
	}
}
print  MSACOLOREDHTML "</font></center>\n<center><table border=0 cols=3 width=310>\n<tr>\n<td align=left><b>Confident</b></td>\n<td align=center><b><---></b></td>\n<td align=right><b>Uncertain</b></td>\n</tr>\n</table></center>\n</td>\n</tr>\n</table>\n";

print MSACOLOREDHTML "</body>\n<html>\n";
close MSACOLOREDHTML;
