#!/usr/bin/perl -w

package BlastResultsForm;

use strict;

#############
#subroutines#
#############

#************************************************************************************************
# read the HTML Format of Blast and creates Form Allowing Seq Selection.
# In: (1) Blast Results in HTML Format (-T T)
#     (2) Num of Blast Rounds
#     (3) Output File
# Out: return OK if Succeed otherwise the failure message

sub MakeForm {
	my $Blast_HTML_In=shift;
	my $Round=shift;
	my $Blast_Form_Out=shift;
	my $Selected_Seq_File=shift;
	my $DB=shift;
	my $Multi_Select=shift; #yes - Multi (default), no - only single
	my $Redirect_Output=shift; # The reditected output
	if (!defined $Multi_Select) {$Multi_Select="yes";}
	if (!defined $DB) {$DB="";} 
	unless (open (BLAST_FORM, ">$Blast_Form_Out")) {return "BlastResultsForm::MakeForm : Could Not Open The Blast Form Output '$Blast_Form_Out' For writing\n"};
	unless (open (BLAST_HTML_RESULTS,$Blast_HTML_In)) {return "BlastResultsForm::MakeForm : Could Not Open The Blast HTML Output '$Blast_HTML_In}' For reading\n"};
	print BLAST_FORM <<EndOfOuput;
<\?php
\$selected_seq = \$_POST[\"MarkedSeq\"];
if (!isset(\$_POST['submit'])) { // if page is not submitted to itself echo the form
\?>

<div id="Select_MSA_Seq">
<form name="Select_Blast_Hit" id="Select_Blast_Hit" method="post">
<H1 align=center><font color='red'> Please choose which sequences you want to use for ConSurf calculation</font></h1>
<H3 align=center><font color='black'> Please noticed that the query sequence is automatically included in the analysis without choosing it</font></h3><BR>
<FONT FACE= "Courier New">
EndOfOuput
	my $Counter_E_Value_Section=1;
	my $Counter_Align_Section=1;
	my $Desired_Round=0;
	my $Header=1; # Flag to indicate blast header
	if ($Round>1) # More than one Round - Go Till the results of the right run
	{
		while (my $line=<BLAST_HTML_RESULTS> and ($Header==1))
		{
			if ($line!~/Results from round/)
			{
				print BLAST_FORM $line;
			}
			elsif ($line=~/Results from round ([0-9]+)/)
			{
				if ($1==$Round)
				{
					$Desired_Round=1;
					print BLAST_FORM $line;	
				}
				$Header=0;
			}
		}
		while (my $line=<BLAST_HTML_RESULTS> and ($Desired_Round==0)) #Go till the desired round
		{
			if ($line =~/Results from round (\d+)/)
			{
				if ($1==$Round)
				{
					$Desired_Round=1;	
					print BLAST_FORM $line;
				}
			}
		}
	}
	while (my $line=<BLAST_HTML_RESULTS>)
	{
		my $printed=0; #Was this line processed
		my $SeqID;		
		if (($line=~/^([A-Za-z0-9_]+)\|(.*)/) or ($line=~/(^UniRef.*)/)){
			if ($line=~/^([A-Za-z0-9]+)\|(.*)\|/){
				$SeqID=$2;
			}
			if (($DB eq "UniProt") or ($DB eq "SWISS-PROT")){
				if ($line =~ m/([A-Za-z0-9]+)\|([A-Za-z0-9]+_[A-Za-z0-9]+)/) {$SeqID = $2;}
			}
			elsif ($DB eq "CLEAN_UNIPROT"){
				if ($line =~ m/(\S+_\S+)\|(\S+)/) {$SeqID = $1;}
			}
			elsif ($DB eq "UNIREF90"){
				if ($line =~ m/(\S+_\S+)/) {$SeqID = $1;}
			}
			
#			print "$SeqID\t";#<STDIN>;
			if ($Counter_E_Value_Section<10) {print BLAST_FORM "$Counter_E_Value_Section&nbsp;&nbsp;&nbsp;";}
			elsif ($Counter_E_Value_Section>9 and $Counter_E_Value_Section<100) {print BLAST_FORM "$Counter_E_Value_Section&nbsp;&nbsp;";}
			elsif ($Counter_E_Value_Section>=100) {print BLAST_FORM "$Counter_E_Value_Section ";}
			if ($Multi_Select eq "no") {print BLAST_FORM "<input type=\"radio\" ";}
			else {print BLAST_FORM "<input type=\"checkbox\" ";}
			print BLAST_FORM "value=\"$SeqID".".BlastHit_$Counter_E_Value_Section\" name=\"MarkedSeq[]\" id=\"$SeqID.evalue\" onchange=\"document.getElementById('$SeqID.align').checked=document.getElementById('$SeqID.evalue').checked\">$line";
			$Counter_E_Value_Section++;
			$printed=1;
		}
		elsif(($DB eq "CULLED_PDB") or ($DB eq "PDB")){
			if ($line =~ m/^([A-Z0-9]{5}) [0-9]+/) {
				$SeqID = $1;
#				print "$SeqID\t";#<STDIN>;
				if ($Counter_E_Value_Section<10) {print BLAST_FORM "$Counter_E_Value_Section&nbsp;&nbsp;&nbsp;";}
				elsif ($Counter_E_Value_Section>9 and $Counter_E_Value_Section<100) {print BLAST_FORM "$Counter_E_Value_Section&nbsp;&nbsp;";}
				elsif ($Counter_E_Value_Section>=100) {print BLAST_FORM "$Counter_E_Value_Section ";}
				if ($Multi_Select eq "no")
				{
					print BLAST_FORM "<input type=\"radio\"  value=\"$SeqID"."\" name=\"MarkedSeq\" id=\"$SeqID.evalue\" >$line";
				}
				else {
					print BLAST_FORM "<input type=\"checkbox\" value=\"$SeqID".".BlastHit_$Counter_E_Value_Section\" name=\"MarkedSeq[]\" id=\"$SeqID.evalue\" onchange=\"document.getElementById('$SeqID.align').checked=document.getElementById('$SeqID.evalue').checked\">$line";
				}
				#print BLAST_FORM "value=\"$SeqID".".BlastHit_$Counter_E_Value_Section\" name=\"MarkedSeq[]\" id=\"$SeqID.evalue\" onchange=\"document.getElementById('$SeqID.align').checked=document.getElementById('$SeqID.evalue').checked\">$line";
				$Counter_E_Value_Section++;
				$printed=1;
			}
		}
		elsif (($line=~/^<a href=(.*)>([A-Za-z0-9_]+)\|(.*)\|(.*)?(\|)?(.*)?<\/a>/) and ($DB eq "NT")){ # for blast agains nt or nr
			$SeqID=$3;
			if ($2 eq "pdb"){$SeqID=$3.$4;}
			if ($Counter_E_Value_Section<10) {print BLAST_FORM "$Counter_E_Value_Section&nbsp;&nbsp;&nbsp;";}
			elsif ($Counter_E_Value_Section>9 and $Counter_E_Value_Section<100) {print BLAST_FORM "$Counter_E_Value_Section&nbsp;&nbsp;";}
			elsif ($Counter_E_Value_Section>=100) {print BLAST_FORM "$Counter_E_Value_Section ";}
			if ($Multi_Select eq "no") {print BLAST_FORM "<input type=\"radio\" ";}
			else {print BLAST_FORM "<input type=\"checkbox\" ";}
			print BLAST_FORM "value=\"$SeqID".".BlastHit_$Counter_E_Value_Section\" name=\"MarkedSeq[]\" id=\"$SeqID.evalue\" onchange=\"document.getElementById('$SeqID.align').checked=document.getElementById('$SeqID.evalue').checked\">$line";
			$Counter_E_Value_Section++;
			$printed=1;
		}
		elsif (($line=~/^<a href=(.*)>([A-Za-z0-9_]+)\|(.*)\|(.*)?(\|)?(.*)?<\/a>/) and ($DB eq "NR_PROT_DB")){ # for blast agains or nr
            $SeqID=$3;
            if ($2 eq "pdb"){$SeqID=$3.$4;}
            if ($Counter_E_Value_Section<10) {print BLAST_FORM "$Counter_E_Value_Section&nbsp;&nbsp;&nbsp;";}
            elsif ($Counter_E_Value_Section>9 and $Counter_E_Value_Section<100) {print BLAST_FORM "$Counter_E_Value_Section&nbsp;&nbsp;";}
            elsif ($Counter_E_Value_Section>=100) {print BLAST_FORM "$Counter_E_Value_Section ";}
            if ($Multi_Select eq "no") {print BLAST_FORM "<input type=\"radio\" ";}
            else {print BLAST_FORM "<input type=\"checkbox\" ";}
            print BLAST_FORM "value=\"$SeqID".".BlastHit_$Counter_E_Value_Section\" name=\"MarkedSeq[]\" id=\"$SeqID.evalue\" onchange=\"document.getElementById('$SeqID.align').checked=document.getElementById('$SeqID.evalue').checked\">$line";
            $Counter_E_Value_Section++;
            $printed=1;
        }

		if ($line=~/^>(<a name = [0-9]+><\/a>)(.*)/){
			$SeqID=$2;
			if (($DB eq "NT")or ($DB eq "NR_PROT_DB"))
			{
				if ($2=~/<a href=(.*)>(.*)\|([A-Za-z0-9._]+)\|(.*)?(\|)?(.*)?/){
					$SeqID=$3;
					if ($2 eq "pdb"){$SeqID=$3.$4;}
				}
			}
			if (($DB eq "UniProt") or ($DB eq "SWISS-PROT")){
				if ($2 =~ m/([A-Za-z0-9]+)\|([A-Za-z0-9]+_[A-Za-z0-9]+)/) {$SeqID = $2;}
			}
			elsif ($DB eq "CLEAN_UNIPROT"){
				if ($2 =~ m/(\S+_\S+)\|(\S+)/) {$SeqID = $1;}
			}
			elsif ($DB eq "UNIREF90"){
				if ($2 =~ m/(\S+_\S+)/) {$SeqID = $1;}
			}
			elsif(($DB eq "CULLED_PDB") or ($DB eq "PDB")){
				if ($2 =~ m/^([A-Z0-9]{5}) [0-9]+/) {$SeqID = $1;}
			}
			#elsif ($line=~/^>(<a name = [0-9]+><\/a>)([A-Za-z0-9]+)\|(.*)\|/){
			# $SeqID=$3;
#			print "$SeqID\t$Counter_Align_Section\n";#<STDIN>;
			if ($Multi_Select eq "no") 
			{
				print BLAST_FORM "$Counter_Align_Section<input type=\"radio\" value=\"$SeqID\" name=\"MarkedSeq\" id=\"$SeqID.align\"\">$line";
			}
			else 
			{
				print BLAST_FORM "$Counter_Align_Section<input type=\"checkbox\" value=\"\" name=\"MarkedSeq[]\" id=\"$SeqID.align\" onchange=\"document.getElementById('$SeqID.evalue').checked=document.getElementById('$SeqID.align').checked\">$line";
			}
			#print BLAST_FORM "value=\"\" name=\"MarkedSeq[]\" id=\"$SeqID.align\" onchange=\"document.getElementById('$SeqID.evalue').checked=document.getElementById('$SeqID.align').checked\">$line"; # this line not enter value to the array because it is synchronized with the e-value section that enter the seqid to the array;
			$Counter_Align_Section++;
			$printed=1;
		}
		if ($printed==0)
		{
			print BLAST_FORM $line;
		}
	}
print BLAST_FORM "</FONT>\n";

if ($Multi_Select eq "yes")
{
	print BLAST_FORM <<EndOfOuput1;

Select the first <INPUT TYPE="text" NAME="FirstHomologs" VALUE="50" size=6> sequences <INPUT TYPE="button" NAME="button1" Value="Update Selection" onClick="SelectFirstSeq(document.Select_Blast_Hit.FirstHomologs.value)">

<input type="button" name="CheckAll" value="Check All" onClick="javascript:checkAll()"> <input type="button" name="UnCheckAll" value="Uncheck All" onClick="javascript:clearAll()">
<input type="submit" value="submit" name="submit" onclick="javascript:return VarifySelect_Box();">
</form>
</div>
<!--Select MSA Seq ends here-->
<SCRIPT LANGUAGE="JavaScript">
<!-- Begin
function checkAll() {
     var form=document.Select_Blast_Hit;
     var boxes = form.getElementsByTagName('input');

     for (var i = 0; i < boxes.length; i++) {
          if (boxes[i].type == 'checkbox'){
                boxes[i].checked = true;
          }
     }
}

function clearAll() {
     if(confirm("Are you sure you want to delete the selected items?"))
     {
        var form=document.Select_Blast_Hit;
     	var boxes = form.getElementsByTagName('input');
	
     	for (var i = 0; i < boxes.length; i++) {
        	  if (boxes[i].type == 'checkbox'){
                	boxes[i].checked = false;
			}
     	}
     }
}

function VarifySelect_Box() {
	var form=document.Select_Blast_Hit;
	var boxes = form.getElementsByTagName('input');
	var NoneChecked=1;
	var Num_of_Seq=0;
	
	for (var i=0; i<boxes.length;i++){
	   if (boxes[i].type == 'checkbox'){
	   	if (boxes[i].checked == true){
		    NoneCheacked=0;
		    Num_of_Seq++;
		    }
		else{
		    NoneCheacked=1;
		    }
	   }
	}
	Num_of_Seq=Num_of_Seq/2;
	if(Num_of_Seq < 1){
	alert("No sequences were selected, please choose sequences and continue;");
	return false;
	}
	if (Num_of_Seq<5)
	{
		alert("Only "+Num_of_Seq+" sequences were selected, please choose at least 5 sequences and continue;");
		return false;
	}
	if (Num_of_Seq<10)
	{
		confirm("Only "+Num_of_Seq+" sequences were selected. It is recomended to select at least 10 sequences for ConSurf analysis. Do you wnat to continue?");
	}
}
//function VarifySelect_Box() {
//	var form=document.Select_Blast_Hit;
//	var boxes = form.getElementsByTagName('input');
//	var NoneChecked=1;
//	
//	for (var i=0; i<boxes.length;i++){
//	   if (boxes[i].type == 'checkbox'){
//	   	if (boxes[i].checked == true){
//		    NoneCheacked=0;
//		    break;
//		    }
//		else{
//		    NoneCheacked=1;
//		    }
//	   }
//	}
//	if(NoneCheacked == 1){
//	alert("No sequences were selected, please choose sequences and continue;");
//	return false;
//	}		
//}


function VarifySelect_Radio() {
	var boxes = document.getElementById('Select_Blast_Hit').getElementsByTagName('input');
	var NoneChecked=1;
	
	for (var i=0; i<boxes.length;i++){
	   if (boxes[i].type == 'radio'){
	   	if (boxes[i].checked == true){
		    NoneCheacked=0;
		    break;
		    }
		else{
		    NoneCheacked=1;
		    }
	   }
	}
	if(NoneCheacked == 1){
	alert("No sequences were selected, please choose sequences and continue;");
	return false;
	}		
}
function SelectFirstSeq (First_Homologs) {
    var form=document.Select_Blast_Hit;
    var boxes = form.getElementsByTagName('input');
    var TotalNum_of_Seq=(boxes.length-5)/2;
    	for (var i = 0; i < First_Homologs; i++) {  
        	  if (boxes[i].type == 'checkbox'){ /select the e-value section checkbox/
                	boxes[i].checked = true;
			}
		  if (boxes[(i+TotalNum_of_Seq)].type == 'checkbox'){ /select the alignment section checkbox/
                	boxes[(i+TotalNum_of_Seq)].checked = true;
			} 	
     	}
}

//  End -->
</script>

<\?
} 
else 
	{
		\$fh = fopen("$Selected_Seq_File", 'w') or die("Can't open file");
		fwrite(\$fh, "The Selected Seq Are:\n");
		foreach (\$selected_seq as \$f) {
			fwrite (\$fh,"\$f\n");
		}
		fclose(\$fh);
		header('Location: $Redirect_Output');
	}
\?>
EndOfOuput1
}

else # Only One seq is selceted
{
	print BLAST_FORM <<EndOfOuput2;
<input type=\"submit\" value=\"submit\" name=\"submit\">
<\?
} 
else 
	{
		\$fh = fopen("$Selected_Seq_File", 'w') or die("Can't open file");
		fwrite (\$fh,"\$selected_seq");
		fclose(\$fh);
	}
\?>
EndOfOuput2
}

	close (BLAST_FORM);
	close (BLAST_HTML_RESULTS);
	return ("OK")
}

1;
