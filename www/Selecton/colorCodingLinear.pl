#!/usr/local/bin/perl

use strict;
use lib "/bioseq/bioSequence_scripts_and_constants";  
use GENERAL_CONSTANTS;

#****************************************

my $run_name = shift;
my $WorkingDir = shift;
my $final_out = $WorkingDir. shift;
my $selectionFile = shift;# the file selection4site, including color range for each residue'
my $areSitesPositive = shift;

my ($html_error, $log_error);

#******************************
#variables

my @Output = ();  # To hold all the information that should be printed in the output file

my @colorstep; #color steps
$colorstep[6] = "#FFBD00"; #"[255,190,0]";      # --> Ka/Ks>1 significant
$colorstep[5] = "#FFFF78"; #"[255,255,120]";    # --> Ka/Ks>1
$colorstep[4] = "#FFFFFF"; #"[255,255,255]";    # --> Ka/Ks<1
$colorstep[3] = "#FCEDF4"; #"[252,237,244]";    # --> Ka/Ks<1
$colorstep[2] = "#FAC9DE"; #"[250,201,222]";    # --> Ka/Ks<1
$colorstep[1] = "#F07DAB"; #"[240,125,171]";    # --> Ka/Ks<1
$colorstep[0] = "#A02560"; #"[160,37,96]";      # --> Ka/Ks<1 significant


my %ColorScale_reverse = (7 => 0,
		  6 => 1,
		  5 => 2,
		  4 => 3,
		  3 => 4,
		  2 => 5,
		  1 => 6,
				);

### read the file with the amino acids and color scores
&readSelection4Site;

### create an html file with the amino acids with background coloring
&createHtml;

#--------------------------------------------------------------
sub readSelection4Site{  
    unless (open GRADES, "<$selectionFile"){
		$html_error = "sys";
		$log_error = "colorCodingLinear.pl : readSelectionFile : Cannot open the file $selectionFile for reading $!";
		&print_error_and_exit;
    } 
    my $CountSeq = 0;
	my $isASitePositive = "no";
    ### the file is in the format : position a.a. colorCode
      while (<GRADES>) {
        my $line = $_;
        chomp $line;
        my $AA;
	    my $color;
		$line =~ /^\s*(\d+)\s+(\S+)\s+(\S+)/;
		$AA = $2;
    	$color = $3;
	if (($color == 1) || ($color == 2)) { # there is a site with positive selection
	    $isASitePositive = "yes";
	}
	
	# ignore lines with '*' - meaning gaps between the seqres and its homologues
	if ($AA eq "*" || $AA eq "_" || $AA eq "-" ||$AA eq ""){
	    next;
	}
	$Output[$CountSeq]{SEQ} = $AA;
	$Output[$CountSeq]{COLOR} = $color;
	$CountSeq++;
   }
   close GRADES;
	unless (open POS, ">$areSitesPositive"){
		$html_error = "sys";
		$log_error = "colorCodingLinear.pl : readSelectionFile : Cannot open the file $areSitesPositive for writing $!";
		&print_error_and_exit;
	}
    print POS $isASitePositive."\n";
    close POS;
}
#--------------------------------------------------------------
sub createHtml{
	system 'echo "(touch '.$final_out.'; chmod oug+rxw '.$final_out.')" | /bin/tcsh';
	unless (open HTML, ">$final_out"){
		$html_error = "sys";
		$log_error = "colorCodingLinear.pl : createHtml: Cannot open the file $final_out for writing $!";
		&print_error_and_exit;
    } 
    print HTML "<html>\n<title>Selecton Results: $run_name</title>\n<body bgcolor='white'><H1 align=center><u>Selecton Results</u></H1><br><br>\n";
    print HTML "<table border=0 width=780>\n";
    my $CountSeq=1;
    foreach my $element (@Output){
	    my $proteinSize=scalar(@Output);
		###prints the amino acid numbers every in jumps of 10
	    my $x=$CountSeq % 50;
	    my $y=$CountSeq % 10;
	    if ($x==1){
		    if ($CountSeq>1) {
		    	print HTML "</tr>\n";
	    	}
		    print HTML "<tr>\n"; 
	    }
		if ($x==1){
			if ($CountSeq>1) {
				print HTML "</tr>";
			}
			print HTML "<tr>";
			for (my $count=0; $count<5; $count++){
				if (($CountSeq+10*$count)<$proteinSize){
					print HTML "<td>\n<font face='Courier New' color='black' size=+1>".($CountSeq+10*$count)."<br></td>\n";
				}
			}
			print HTML "<tr>";
		}	
		my $color = $$element{COLOR};
		my $amino= $$element{SEQ};
		if ($y==1){
			if ($CountSeq>1) {
				print HTML "</td>";
			}
			print HTML "<td>";
		}
		print HTML "<b><font face='Courier New' color='black' size=+1><span style='background: $colorstep[$ColorScale_reverse{$color}];'>$amino</span></font></b>\n";
	    $CountSeq=$CountSeq+1;
	    
    }
    print HTML "</td></tr></table>";
    print HTML "<br><br><br><b><u>Legend:</u><br><br>\n";
	print HTML "The selection scale:</b><br>\n";
	print HTML "<table border=0 cols=1 width=310>\n";
	print HTML "<tr><td align=left>\n";
	print HTML "<font face='Courier New' color='black' size=+1>\n";
	print HTML "<span style='background: $colorstep[$ColorScale_reverse{1}];'>&nbsp;1&nbsp;</span><span style='background: $colorstep[$ColorScale_reverse{2}];'>&nbsp;2&nbsp;</span><span style='background: $colorstep[$ColorScale_reverse{3}];'>&nbsp;3&nbsp;</span><span style='background: $colorstep[$ColorScale_reverse{4}];'>&nbsp;4&nbsp;</span><span style='background: $colorstep[$ColorScale_reverse{5}];'>&nbsp;5&nbsp;</span><span style='background: $colorstep[$ColorScale_reverse{6}];'>&nbsp;6&nbsp;</span><span style='background: $colorstep[$ColorScale_reverse{7}];'>&nbsp;7&nbsp;</span></font></font></center>\n";
	print HTML "</td></tr><tr><td>Positive selection&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Purifying selection";
	print HTML "</td></tr></table>";    
    print HTML "</body></html>";
    close HTML;
}
#--------------------------------------------------------------
sub print_error_and_exit{
    open ERROR, ">".$WorkingDir."error";
    print ERROR "HTML: $html_error\n";
    print ERROR "LOG: $log_error\n";
    close ERROR;
    exit;
}
#--------------------------------------------------------------
