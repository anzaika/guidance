#!/usr/bin/perl -w

package ConSeq_gradesPE_and_Outputs;
use lib "/bioseq/ConSurf";
use CONSURF_CONSTANTS;

use Bio::SearchIO;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Align::AlignI;

my %consurf_html_colors=(	"1"=>"#10C8D1",		#less conserved
				"2"=>"#8CFFFF",
				"3"=>"#D7FFFF",
				"4"=>"#EAFFFF",	
				"5"=>"#FFFFFF",       #average
				"6"=>"#FCEDF4",
				"7"=>"#FAC9DE",
				"8"=>"#F07DAB",
                   		"9"=>"#A02560",		#most conserved
				"ISD"=>"yellow");		#InSufficient Data
    
#####################################################################
# Print the results: the colored sequence and the B/E information
#####################################################################
sub ConSeq_HTML_Output { 

    my $output = shift; # An array of Hashes that contains all the Scores data
    my $ref_Solv_Acc_Pred=shift; # An hash that holds all the solvent accessibility prediction data
    my $Out=shift;
    my $Just_Seq_Flag=shift; # Optional, if equal to yes just colored seq is printed, default is no
    unless ($Just_Seq_Flag eq "yes"){$Just_Seq_Flag="no";}
    #print LOG "\ncreate_output: cd $WorkingDir ; touch $colors_file ; chmod oug+wrx $colors_file\n";
    #system 'echo "(cd '.$WorkingDir.'; /bin/touch '.$colors_file.'; chmod oug+wrx '.$colors_file.')" | /bin/tcsh';
    unless (open COLORS, ">$Out"){
        return "ConSeq_gradesPE_and_Outputs::ConSeq_HTML_Output: Can\'t open the file $Out";}
    
    print COLORS "<html>\n<title>ConSeq Results</title>\n"; 
    print COLORS "<body bgcolor='white'>\n";
    print COLORS "<H1 align=center><u>ConSeq Results</u></H1>\n\n";

    print COLORS "\n<table border=0 width=750>\n";
    print COLORS "<tr><td>\n";

    ### print the colored sequence
    my $count = 1;
    my $letter_str = "";
    my $pred_str = "";
    my $func_str = "";
    foreach my $elem (@{$output}){

       # print the counter above the beginning of each 10 characters
       if ($count % 50 == 1){
           
           my $count_num = $count;
           for ($count_num ; $count_num < $count+50 ; $count_num += 10){

               if ($count_num <= @{$output}){

                   my $space_num = 11 - length($count_num);
                   my $spaces = "&nbsp;" x $space_num;

                   print COLORS "<font face='Courier New' color='black' size=+1>".$count_num.$spaces."</font>";
               }
            }
           print COLORS "<br>\n";
       }

       ### print the colored letters and 'e' for the exposed residues

       # after 50 characters - print newline
       if ($count % 50 == 0 or $count == @{$output})
       {
            if ($$elem{ISD}==1){ 	#INSUFFICIENT DATA
               $letter_str .= "<b><font face='Courier New' color='$consurf_html_colors{ISD}' size=+1><span style='background: $consurf_html_colors{$$elem{COLOR}};'>$$elem{SEQ}</span></font></b><br>";
            }
            else{
                if ($$elem{COLOR}==9){ # MOST CONSERVED
                    $letter_str .= "<b><font face='Courier New' color='white' size=+1><span style='background: $consurf_html_colors{$$elem{COLOR}};'>$$elem{SEQ}</span></font></b><br>\n";
                }
                else {
                    $letter_str .= "<b><font face='Courier New' color='black' size=+1><span style='background: $consurf_html_colors{$$elem{COLOR}};'>$$elem{SEQ}</span></font></b><br>\n";
                }
            }

	if ($Just_Seq_Flag eq "no"){
            # exposed
            if ($$ref_Solv_Acc_Pred{$$elem{POS}} eq "e"){
               $pred_str .= "<b><font face='Courier New' color='orange' size=+1>e</font></b><br>\n";
               
               # exposed & conserved = functional
               if ($$elem{COLOR} == 9 || $$elem{COLOR} == 8){

                   $func_str .= "<b><font face='Courier New' color='red' size=+1>f</font></b>\n";
               }
               # not functional
               else {

                   $func_str .= "<font face='Courier New' size=+1>&nbsp;</font>\n";
               }
            }
            # burried
            else {
                $pred_str .= "<b><font face='Courier New' color='#00cc00' size=+1>b</font></b><br>\n";

                # burried & conserved = structural
                if ($$elem{COLOR} == 9 ){
                    $func_str .= "<b><font face='Courier New' color='#000099' size=+1>s</font></b>\n";
                }
                # not structural
                else {
                    $func_str .= "<font face='Courier New' size=+1>&nbsp;</font>\n";
                }
            }
	  }
            print COLORS $letter_str;
            print COLORS $pred_str if ($Just_Seq_Flag eq "no"); 
            print COLORS $func_str if ($Just_Seq_Flag eq "no");
            print COLORS "</td></tr>\n";
            print COLORS "<tr><td>\n";
 
            $letter_str = "";
            $pred_str = "";
            $func_str = "";
        }
        # after 10 characters - print a space ('&nbsp;')
        elsif ($count % 10 == 0) {
            if ($$elem{ISD} eq '1'){
                $letter_str .= "<b><font face='Courier New' color='$consurf_html_colors{ISD}' size=+1><span style='background: $consurf_html_colors{$$elem{COLOR}};'>$$elem{SEQ}</span> </font></b>\n";
            }
            else{
                if ($consurf_html_colors{$$elem{COLOR}} eq "#A02560"){
                    $letter_str .= "<b><font face='Courier New' color='white' size=+1><span style='background: $consurf_html_colors{$$elem{COLOR}};'>$$elem{SEQ}</span> </font></b>\n";
                }
                else {
                $letter_str .= "<b><font face='Courier New' color='black' size=+1><span style='background: $consurf_html_colors{$$elem{COLOR}};'>$$elem{SEQ}</span> </font></b>\n";
                }
            }
            if ($Just_Seq_Flag eq "no"){
	    # exposed
            if ($$ref_Solv_Acc_Pred{$$elem{POS}} eq "e"){
                $pred_str .= "<b><font face='Courier New' color='orange' size=+1>e&nbsp;</font></b>";
               
                # exposed & conserved = functional
                if ($$elem{COLOR} == 9 || $$elem{COLOR} == 8){

                   $func_str .= "<b><font face='Courier New' color='red' size=+1>f&nbsp;</font></b>";
                }
                # not functional
                else {

                   $func_str .= "<font face='Courier New' size=+1>&nbsp;&nbsp;</font>";
                }
            }
            # burried
            else {
                $pred_str .= "<b><font face='Courier New' color='#00cc00' size=+1>b&nbsp;</font></b>";

                # burried & conserved = structural
                if ($$elem{COLOR} == 9) {

                   $func_str .= "<b><font face='Courier New' color='#000099' size=+1>s&nbsp;</font></b>";
                }
                 # not structural
                else {
                    $func_str .= "<font face='Courier New' size=+1>&nbsp;&nbsp;</font>";
                }
            }
        }
      }
        else {
            if ($$elem{ISD} eq '1'){
                $letter_str .= "<b><font face='Courier New' color='$consurf_html_colors{ISD}' size=+1><span style='background: $consurf_html_colors{$$elem{COLOR}};'>$$elem{SEQ}</span></font></b>";
            }
            else{
                if ($consurf_html_colors{$$elem{COLOR}} eq "#A02560"){
                    $letter_str .= "<b><font face='Courier New' color='white' size=+1><span style='background: $consurf_html_colors{$$elem{COLOR}};'>$$elem{SEQ}</span></font></b>";
                }
                else {
                    $letter_str .= "<b><font face='Courier New' color='black' size=+1><span style='background: $consurf_html_colors{$$elem{COLOR}};'>$$elem{SEQ}</span></font></b>";
                }
            }

       if ($Just_Seq_Flag eq "no"){
            # exposed
            if ($$ref_Solv_Acc_Pred{$$elem{POS}} eq "e"){
                $pred_str .= "<b><font face='Courier New' color='orange' size=+1>e</font></b>";

               # exposed & conserved = functional
               if ($$elem{COLOR} == 9 || $$elem{COLOR} == 8 ){

                   $func_str .= "<b><font face='Courier New' color='red' size=+1>f</font></b>";
               }
               # not functional
               else {

                   $func_str .= "<font face='Courier New' size=+1>&nbsp;</font>";
               }
            }
            #burried
            else {
                $pred_str .= "<b><font face='Courier New' color='#00cc00' size=+1>b</font></b>";

               # burried & conserved = structural
               if ($$elem{COLOR} == 9 ){

                   $func_str .= "<b><font face='Courier New' color='#000099' size=+1>s</font></b>";
               }
               # not structural
               else {

                   $func_str .= "<font face='Courier New' size=+1>&nbsp;</font>";
               }
            }
        }
      }
        $count++;
    }

    print COLORS "</td></tr>\n</table><br>\n";

    # print the color scale
    print COLORS "\n<br><b><u>Legend:</u><br><br>\nThe conservation scale:</b><br>\n<table border=0 cols=1 width=310>\n<tr><td align=center>\n<font face='Courier New' color='black' size=+1><center>\n";
    for (my $i=1 ; $i<=9 ; $i++){

       if ($i == 9){

           print COLORS "<font color='white'><span style='background: $consurf_html_colors{$i};'>&nbsp;$i&nbsp;</span></font>";
       }
       else {

           print COLORS "<span style='background: $consurf_html_colors{$i};'>&nbsp;$i&nbsp;</span>";
       }
    }
    print COLORS "</font></center>\n<center><table border=0 cols=3 width=310>\n<tr>\n<td align=left><b>Variable</b></td>\n<td align=center><b>Average</b></td>\n<td align=right><b>Conserved</b></td>\n</tr>\n</table></center>\n</td>\n</tr>\n</table>\n";
 if ($Just_Seq_Flag eq "no"){
    print COLORS "<b><font face='Courier New' color='orange' size=+1>e</font> - An exposed residue according to the neural-network algorithm.</b><br>\n";
    print COLORS "<b><font face='Courier New' color='#00cc00' size=+1>b</font> - A buried residue according to the neural-network algorithm.</b><br>\n";
    print COLORS "<b><font face='Courier New' color='red' size=+1>f</font> - A predicted functional residue (highly conserved and exposed).</b><br>\n";
    print COLORS "<b><font face='Courier New' color='#000099' size=+1>s</font> - A predicted structural residue (highly conserved and buried).</b><br>\n";
  }
   print COLORS "<b><font face='Courier New' color='$consurf_html_colors{ISD}' size=+1><span style='background: $consurf_html_colors{2};'>X</span></font> - Insufficient data - the calculation for this site was performed on less than 10% of the sequences.</font></b><br>\n"; 

    print COLORS "</body>\n</html>\n";

    close COLORS;
    return ("ok");
}


###################################################################################
# Print the results: the colored sequence of the msa acording to the query sequence
###################################################################################
sub print_msa_colors{ 
    
    my $output = shift;  # An array of Hashes that contains all the Scores data
    my $MSA=shift; # MSA in ClustalW format
    my $SeqName=shift; # The reference seq on the MSA
    my $Out_file=shift; # html File with the Colored MSA
#    my $TEST_QA="/bioseq/data/results/ConSurf/1251734866/qsub.QA";
#    open (QA,">$TEST_QA");
    my %msaColors = ();
    my %msaPrintColors = ();
    my $lineCounter = 0;
    my @line;
    my $key;
    my $query = 'query';
    my @querySequence;
    my $pdbIdLengthForPrinting;
    my $fontSize=2;
    my $sequenceLengthForDisplay=40;
    my @msaRightOrder=0;
    my $msaRightOrderCounter=0;
    my $tdWidth = 5;
    
    # counts how many times we print the whole seqtion(relevants to sequences longer than the sequenceLengthForDisplay) 
    my $msaCounter=0; 
    
    
    open (MSAFILE,$MSA);
    
    while (<MSAFILE>)#loop that matches the pdb_id and the sequence 
    {
    
      my $line=$_;  
       if ($line =~ /\s*(.+)[\s|\t]+([\w|\-]+)\s*/){ 
        unless ($msaColors{$1})  
        {
          if   ($1 !~ /\s*CLUSTAL\s*/)
          {   
#             print QA "MSA_SEQ:$1\n";
	     $msaRightOrder[$msaRightOrderCounter] = $1;
             $msaRightOrderCounter++;
          } 
        } 
          $msaColors{$1} .= $2 if $1 !~ /\s*CLUSTAL\s*/;            
        }    
    }
     
    
    foreach my $msaKey (%msaColors)
    {
        if ($msaKey =~  /\s*$SeqName\s*/){
            @line = split(//,$msaColors{$msaKey});    
        }
    }   
    
    # @line = split(//,$msaColors{$query});
    foreach my $elem (@{$output})
    {
    
       while ($line[$lineCounter] eq '-')
       {
               $querySequence[$lineCounter]{SEQ} = '-';
               $lineCounter++;
       }
          
      $querySequence[$lineCounter]{ISD} = $$elem{ISD};  #     
      $querySequence[$lineCounter]{COLOR} = $consurf_html_colors{$$elem{COLOR}};    #
      $querySequence[$lineCounter]{SEQ} = $$elem{SEQ};    
#    	print QA "SEQ:$$elem{SEQ}\tCOL:$querySequence[$lineCounter]{COLOR}\tISD:$$elem{ISD}\n";
      $lineCounter++;      
    }
    foreach $key (keys %msaColors)
    {
        @line = split(//,$key);
       if ($pdbIdLengthForPrinting <= $#line)
       {  
        $pdbIdLengthForPrinting = $#line; 
        }
    }
            
    $msaCounter = $#querySequence / $sequenceLengthForDisplay;
    $msaCounter =~ s/\..*//;
    foreach $key (keys %msaColors)#loop that runs over every sequence of the msa
    {
          $lineCounter = 0;
          foreach my $elem (@querySequence)#loop that runs over every amino acid of the sequence
          {
            @line = split(//,$msaColors{$key});
#	    foreach my $qa_key ( keys %{$elem}) {print QA " $elem->{COLOR}\t";}print QA "\n";
            if ((exists $elem->{COLOR} && $elem->{COLOR} eq '-') || (exists $line[$lineCounter] && $line[$lineCounter] eq '-'))
            {
                $msaPrintColors{$key} .= "<td width=$tdWidth><font face='Courier New' color='black' size=$fontSize><span style='background: white;'>$line[$lineCounter]</span></font></td>@\n";
            
            }  
            else
            {
               if (exists $$elem{ISD} && $$elem{ISD} eq '1')
               {                 
                  #$msaPrintColors{$key} .= "<td width=$tdWidth><b><font face='Courier New' color='yellow' size=$fontSize><span style='background: $consurf_html_colors{$$elem{COLOR}};'>$line[$lineCounter]</span></font></b></td>@\n";
		  $msaPrintColors{$key} .= "<td width=$tdWidth><b><font face='Courier New' color='yellow' size=$fontSize><span style='background: $elem->{COLOR};'>$line[$lineCounter]</span></font></b></td>@\n";
               }
               else
               {
                if (exists $elem->{COLOR}){
                    if($elem->{COLOR} eq "#A02560"){     
                        $msaPrintColors{$key} .= "<td width=$tdWidth><b><font face='Courier New' color='white' size=$fontSize><span style='background: $elem->{COLOR};'>$line[$lineCounter]</span></font></b></td>@\n";
                    }
                    else {                
                        $msaPrintColors{$key} .=  "<td width=$tdWidth><b><font face='Courier New' color='black' size=$fontSize><span style='background: $elem->{COLOR};'>$line[$lineCounter]</span></font></b></td>@\n";
                    }
                }
                else{
                    $msaPrintColors{$key} .=  "<td width=$tdWidth><b><font face='Courier New' color='black' size=$fontSize><span style='background: #FFFFFF;'>$line[$lineCounter]</span></font></b></td>@\n";
                }
               }
            }
          
            $lineCounter++;
          
          }    
        }
      #  print LOG "\nCreate colored msa file: cd $WorkingDir ; touch $SeqFile ; chmod oug+wrx $SeqFile\n";
      #  print LOG "\nSeqName : $SeqName \n";
    
      #  system 'echo "(cd '.$WorkingDir.'; /bin/touch '.$msa_colored_file.'; chmod oug+wrx '.$msa_colored_file.')" | /bin/tcsh';
        open (MSACOLOREDHTML, ">".$Out_file);
        print MSACOLOREDHTML "<html>\n<head>\n</head>\n<body>\n";
        print MSACOLOREDHTML "<H1 align=center><u>Color-Coded MSA</u></H1>\n\n";
        print MSACOLOREDHTML "<table border=0  CELLSPACING=1  CELLPADDING=0 >";
            
        for (my $j=0;$j<=$msaCounter;$j++)
        { 
            #     foreach $key (keys %msaPrintColors)#print the pdb_id and sequences
            for (my $k=0;$k<=$#msaRightOrder;$k++) 
            {
                $key = $msaRightOrder[$k]; 
                my $keyP=$key;        
                $keyP =~ s/\s*$//; 
                my  $query_spaces = ConSeq_gradesPE_and_Outputs::printSpaces($key,$pdbIdLengthForPrinting,$fontSize); 
                if ($key =~  /\s*$SeqName\s*/)
                {
                    print MSACOLOREDHTML "<tr><td><b><font face='Courier New' color='black' size=$fontSize><u>$keyP</u>" . $query_spaces . "    </font></b></td>";
                }
                else
                {
                    print MSACOLOREDHTML "<tr><td><b><font face='Courier New' color='black' size=$fontSize>$keyP" . $query_spaces . "    </font></b></td>"; 
                }
                 @line = split(/@/,$msaPrintColors{$key});
                for (my $i=$j*$sequenceLengthForDisplay;$i<($j+1)*$sequenceLengthForDisplay;$i++)
                {   
                   print MSACOLOREDHTML  $line[$i] if exists $line[$i];
                }
                print MSACOLOREDHTML '</tr><tr>';
          
            }
            print MSACOLOREDHTML '</tr><tr><td>&nbsp;</td></tr>'; 
        }
        print MSACOLOREDHTML "</table>";
	if ($MSA=~/([^\/]+)$/g)
	{
        	print MSACOLOREDHTML "<p><A HREF=$1 TARGET=MSA_window> Link to the Clustal formatted alignment</A></p>";
	}
        else
	{
		print MSACOLOREDHTML "<p><A HREF=$MSA TARGET=MSA_window> Link to the Clustal formatted alignment</A></p>";
	}
      # print the color scale
      print  MSACOLOREDHTML "\n<br><b><u>Legend:</u><br><br>\nThe conservation scale:</b><br>\n<table border=0 cols=1 width=310>\n<tr><td align=center>\n<font face='Courier New' color='black' size=+1><center>\n";
      for (my $i=0 ; $i<=9 ; $i++){
    
         if ($i == 9){
    
             print  MSACOLOREDHTML "<font color='white'><span style='background: $consurf_html_colors{$i} ;'>&nbsp;$i&nbsp;</span></font>";
         }
         else {
    
             print  MSACOLOREDHTML "<span style='background: $consurf_html_colors{$i};'>&nbsp;$i&nbsp;</span>";
         }
    }
      print  MSACOLOREDHTML "</font></center>\n<center><table border=0 cols=3 width=310>\n<tr>\n<td align=left><b>Variable</b></td>\n<td align=center><b>Average</b></td>\n<td align=right><b>Conserved</b></td>\n</tr>\n</table></center>\n</td>\n</tr>\n</table>\n";
    
     print  MSACOLOREDHTML "<b><font face='Courier New' color='yellow' size=+1><span style='background: $consurf_html_colors{9};'>X</span></font> - Insufficient data - the calculation for this site was performed on less than 10% of the sequences.</b><br>\n"; 
        print MSACOLOREDHTML "</body>\n<html>\n";
}
sub printSpaces 
{
 #check that the lenght (with apaces) of all the pdb_id are the same - for display  
    my $key = shift;
    my $pdbIdLengthForPrinting =shift;
    my $fontSize = shift;
    my $query_spaces = "<font face='Courier New'  size=$fontSize>&nbsp;";
    $key =~ s/\s*$//;
    my  @keyLine = split(//,$key);
    if ($#keyLine <= $pdbIdLengthForPrinting){
        for (my $i=1;$i <= $pdbIdLengthForPrinting  - $#keyLine;$i++)
        { 
            $query_spaces .= "&nbsp;"; 
        }
    }
    $query_spaces .= "</font>";         
    return  $query_spaces;
}

#++++++++++++++++++++++++++++++++++
sub ConSeq_NUC_PDF_Output { 
   
    my $Grades_PE = shift; # A GradesPE file
    my $Out=shift;
    my $Protein_Length=shift;
    my $ConSeq_PDF_Lib=CONSURF_CONSTANTS::CONSURF_PDF_LIB;
    #print LOG "\ncreate_output: cd $WorkingDir ; touch $colors_file ; chmod oug+wrx $colors_file\n";
    #system 'echo "(cd '.$WorkingDir.'; /bin/touch '.$colors_file.'; chmod oug+wrx '.$colors_file.')" | /bin/tcsh';
    unless (open PDF, ">$Out"){
        return "ConSeq_gradesPE_and_Outputs::ConSeq_PDF_Output: Can\'t open the file $Out";}
#    unless (open GRADES,"$Grades_PE"){ 
#    	return "ConSeq_gradesPE_and_Outputs::ConSeq_PDF_Output: Can\'t open the GradesPE file \'$Grades_PE\' $!";}
    my %ConSurf_Grades=();
    my @ans=ConSeq_gradesPE_and_Outputs::read_ConSeq_Grades_PE($Grades_PE,\%ConSurf_Grades);
    if ($ans[0] ne "ok") {return "ConSeq_gradesPE_and_Outputs::ConSeq_PDF_Output: \'".@ans."\'";}
    my $IS_THERE_FUNCT_RES=$ans[1];
    my $IS_THERE_STRUCT_RES=$ans[2];
    my $IS_THERE_INSUFFICIENT_DATA=$ans[3];
    
    print PDF <<EndOfHeader;	
<?php
require("$ConSeq_PDF_Lib");
\$pdf=new ConSurf_PDF();
\$pdf->AddPage();
\$pdf->SetFont('Times','BU',20);
\$pdf->Cell(0,0,'ConSurf Results',0,0,C);
\$pdf->SetY(\$pdf->GetY()+10);
EndOfHeader

my $Pos=1;
for ($Pos;$Pos<=$Protein_Length;$Pos++)
{
	
	if (($Pos-1) % 50 == 0)
	{
		print PDF "\$pdf->Ln();\n\$pdf->Ln();\n\$pdf->Ln();\n";
		print PDF "\$Rows_Pos=\$pdf->GetY();\n"; # Line Number
	}
#	elsif (($Pos-1) % 10 == 0)
#	{
#		print PDF "\$pdf->Print_ForegroundColor\("."''".",'B',10,"."0,4"."\);\n";
#	}
	if (($Pos-1) % 10 == 0)
	{
		print PDF "\$pdf->Print_2_Lines_Element\(\$Rows_Pos,".$ConSurf_Grades{$Pos}{'COLOR'}.",".$ConSurf_Grades{$Pos}{'AA'}.",".$Pos.","."10,".$ConSurf_Grades{$Pos}{'Insufficient_Data'}."\);\n";
	}
	else
	{
		print PDF "\$pdf->Print_2_Lines_Element\(\$Rows_Pos,".$ConSurf_Grades{$Pos}{'COLOR'}.",".$ConSurf_Grades{$Pos}{'AA'}.",'',"."10,".$ConSurf_Grades{$Pos}{'Insufficient_Data'}."\);\n";
	}
}

#print PDF <<EndBlock;
print PDF "\$pdf->Ln();\n";
print PDF "\$pdf->Ln();\n";
print PDF "\$pdf->Ln();\n";
print PDF "\$pdf->Ln();\n";
print PDF "\$pdf->Print_NEW_Legend_Nuc($IS_THERE_INSUFFICIENT_DATA);\n"; 
print PDF "\$pdf->Output();\n";
print PDF "?>\n";
print PDF " </span><font color='white'><span style='background: #A02560;'>&nbsp;&nbsp;</span></font></font></center>\n";

#EndBlock	
}
#++++++++++++++++++++++++++++++++++
sub ConSeq_PDF_Output { 
   
    my $Grades_PE = shift; # A GradesPE file
    my $Out=shift;
    my $Protein_Length=shift;
    my $ConSeq_PDF_Lib=CONSURF_CONSTANTS::CONSURF_PDF_LIB;
    #print LOG "\ncreate_output: cd $WorkingDir ; touch $colors_file ; chmod oug+wrx $colors_file\n";
    #system 'echo "(cd '.$WorkingDir.'; /bin/touch '.$colors_file.'; chmod oug+wrx '.$colors_file.')" | /bin/tcsh';
    unless (open PDF, ">$Out"){
        return "ConSeq_gradesPE_and_Outputs::ConSeq_PDF_Output: Can\'t open the file $Out";}
#    unless (open GRADES,"$Grades_PE"){ 
#    	return "ConSeq_gradesPE_and_Outputs::ConSeq_PDF_Output: Can\'t open the GradesPE file \'$Grades_PE\' $!";}
    my %ConSurf_Grades=();
    my @ans=ConSeq_gradesPE_and_Outputs::read_ConSeq_Grades_PE($Grades_PE,\%ConSurf_Grades);
    if ($ans[0] ne "ok") {return "ConSeq_gradesPE_and_Outputs::ConSeq_PDF_Output: \'".@ans."\'";}
    my $IS_THERE_FUNCT_RES=$ans[1];
    my $IS_THERE_STRUCT_RES=$ans[2];
    my $IS_THERE_INSUFFICIENT_DATA=$ans[3];
    
    print PDF <<EndOfHeader;	
<?php
require("$ConSeq_PDF_Lib");
\$pdf=new ConSurf_PDF();
\$pdf->AddPage();
\$pdf->SetFont('Times','BU',20);
\$pdf->Cell(0,0,'ConSurf Results',0,0,C);
\$pdf->SetY(\$pdf->GetY()+10);
EndOfHeader

my $Pos=1;
for ($Pos;$Pos<=$Protein_Length;$Pos++)
{
	
	if (($Pos-1) % 50 == 0)
	{
		print PDF "\$pdf->Ln();\n\$pdf->Ln();\n\$pdf->Ln();\n\$pdf->Ln();\n";
		print PDF "\$Rows_Pos=\$pdf->GetY();\n"; # Line Number
	}
	elsif (($Pos-1) % 10 == 0)
	{
		print PDF "\$pdf->Print_ForegroundColor\("."''".",'B',10,"."0,4"."\);\n";
	}
	if (($Pos-1) % 10 == 0)
	{
		print PDF "\$pdf->Print_4_Lines_Element\(\$Rows_Pos,".$ConSurf_Grades{$Pos}{'COLOR'}.",".$ConSurf_Grades{$Pos}{'AA'}.",".$Pos.",".$ConSurf_Grades{$Pos}{'B_E'}.",10,".$ConSurf_Grades{$Pos}{'Insufficient_Data'}.",'".$ConSurf_Grades{$Pos}{'Struct_Funct'}."'\);\n";
	}
	else
	{
		print PDF "\$pdf->Print_4_Lines_Element\(\$Rows_Pos,".$ConSurf_Grades{$Pos}{'COLOR'}.",".$ConSurf_Grades{$Pos}{'AA'}.",'',".$ConSurf_Grades{$Pos}{'B_E'}.",10,".$ConSurf_Grades{$Pos}{'Insufficient_Data'}.",'".$ConSurf_Grades{$Pos}{'Struct_Funct'}."'\);\n";
	}
}

#print PDF <<EndBlock;
print PDF "\$pdf->Ln();\n";
print PDF "\$pdf->Ln();\n";
print PDF "\$pdf->Ln();\n";
print PDF "\$pdf->Ln();\n";
print PDF "\$pdf->Print_NEW_Legend($IS_THERE_FUNCT_RES,$IS_THERE_STRUCT_RES,$IS_THERE_INSUFFICIENT_DATA);\n"; 
print PDF "\$pdf->Output();\n";
print PDF "?>\n";
print PDF " </span><font color='white'><span style='background: #A02560;'>&nbsp;&nbsp;</span></font></font></center>\n";

#EndBlock	
}
    
sub read_ConSeq_Grades_PE {
    my $Grades_PE = shift; # A GradesPE file
    my $ref_GradesPE_Hash=shift;
    
    my $IS_THERE_FUNCT_RES=0;
    my $IS_THERE_STRUCT_RES=0;
    my $IS_THERE_INSUFFICIENT_DATA=0;
    
    unless (open GRADES,"$Grades_PE"){ 
    	return "ConSeq_gradesPE_and_Outputs::ConSeq_PDF_Output: Can\'t open the GradesPE file \'$Grades_PE\' $!";}
    my $Protein_Length=1;
    while (my $line=<GRADES>)
    {
	chomp ($line);
	my ($pos,$AA,$Score,$Color,$B_E,$S_F)=0;
	if ((length $line)>=50)
	{
		my $pos=substr($line,0,4);
		$pos=trim($pos);
		my $AA=substr($line,5,4);
		my $Score=substr($line,9,7);
		my $Color=substr($line,20,2);
		my $B_E=substr($line,50,2);
		my $S_F=substr($line,59,1);
	
		if ($pos=~/[0-9]+/) # Score Line
		{
			$ref_GradesPE_Hash->{$pos}->{'AA'}=trim($AA);
			$ref_GradesPE_Hash->{$pos}->{'SCORE'}=trim($Score);
			$ref_GradesPE_Hash->{$pos}->{'COLOR'}=trim($Color);
			$ref_GradesPE_Hash->{$pos}->{'B_E'}=trim($B_E);
			$ref_GradesPE_Hash->{$pos}->{'Struct_Funct'}=trim($S_F);
			if ($ref_GradesPE_Hash->{$pos}->{'Struct_Funct'} eq "s")
			{
				$IS_THERE_STRUCT_RES=1;
			}
			elsif ($ref_GradesPE_Hash->{$pos}->{'Struct_Funct'} eq "f")
			{
				$IS_THERE_FUNCT_RES=1;
			}
			if ($Color=~/([0-9])[*]/)
			{		
				$ref_GradesPE_Hash->{$pos}->{'Insufficient_Data'}=1;  #Insufficient Data
				$ref_GradesPE_Hash->{$pos}->{'COLOR'}=$1;
				$IS_THERE_INSUFFICIENT_DATA=1;
			}
			else
			{
				$ref_GradesPE_Hash->{$pos}->{'Insufficient_Data'}=0; 
			}
			if (defined $line[0])
			{
				if ($line[0]>$Protein_Length) {$Protein_Length=$line[0];}
			}
			#print "POS:$pos*\tAA:$AA*\tScore:$Score*\tColor:$Color*\tB/E:$B_E*\tS_F:$S_F*\n";
			#print 	$pos,"\t",$ref_GradesPE_Hash->{$pos}->{'AA'},"\t",$ref_GradesPE_Hash->{$pos}->{'SCORE'},"\t",$ref_GradesPE_Hash->{$pos}->{'COLOR'},"\t",$ref_GradesPE_Hash->{$pos}->{'B_E'},"\t",$ref_GradesPE_Hash->{$pos}->{'Struct_Funct'},"\t",$ref_GradesPE_Hash->{$pos}->{'Insufficient_Data'},"\n"; 
		}
	}
    }
    return ("ok",$IS_THERE_FUNCT_RES,$IS_THERE_STRUCT_RES,$IS_THERE_INSUFFICIENT_DATA);
}

#++++++++++++++++++++++++++++++++++
sub Find_Good_Templates { 
	my $query=shift;
	my $Blast_vs_PDB=shift;
	my $Min_Overlap_Percent=shift;
	my $Min_ID=shift;
	my $Selected_Templates_Nems_ref=shift;
	my $Selected_Templates_Details_ref=shift;
	$Min_ID=$Min_ID*100;
	print "$query\t$Blast_vs_PDB\t$Min_Overlap_Percent\t$Min_ID\n";
	my %templates = ();
	# hash of templates names. each sequence name (unique) points to array of hashes.
    	# {seq_name} => [{e_val => <e_val1>, AAseq => <AAseq1>, beg => <S1_beg>, end => <S1_end>},...,{e_val => <e_valn>, AAseq => <AAseqn>, beg => <Sn_beg>, end => <Sn_end>}]

	my $query_seq_length;
	unless (open QUERY, $query){
        	return ("sys", "ConSeq_gradesPE_and_Outputs::Find_Good_Templates : can't open file $query for reading\n");
    	}
    ################################
    # Extracting Query information #
    ################################
    while (<QUERY>){
        if ($_ !~ m/>/){
            $_=~m/(\S+)/;
            $query_seq_length += length($1);
        }
    }
    close QUERY;

    ######################################################
    # defining the minimum length a homoloug should have #
    ######################################################
    $query_min_length = $query_seq_length*$Min_Overlap_Percent; #min length of overlap with the query
    
    my $searchio = new Bio::SearchIO(-format => 'blast',
                                     -file   => $Blast_vs_PDB);
    while( my $result = $searchio->next_result ) {  # $result is a Bio::Search::Result::ResultI object
        while( my $hit = $result->next_hit ) { # $hit is a Bio::Search::Hit::HitI object or undef if there are no more
            $s_description=$hit->description();
	    if ($s_description=~/\<UNP ([A-Za-z0-9]+_[A-Za-z0-9]+)\>/){$s_accession=$1;}
	    #$s_accession=join (" ",$hit->each_accession_number());
            $s_name= $hit->name();
            while( my $hsp = $hit->next_hsp ) { #hsp is the next available High Scoring Pair, Bio::Search::HSP::HSPI object or null if finished
                # extracting relevant details from the fragment
                ($s_beg, $s_end) = $hsp->range("sbjct");
		($q_beg, $q_end) = $hsp->range("query");
                $AAseq = $hsp->hit_string();
                $AAseq =~ s/-//g;
                $s_eval = $hsp->evalue();
                $s_eval =~ s/,//g;
                if ($s_eval =~ m/^e/) {$s_eval = "1".$s_eval;}
                $s_ident = $hsp->percent_identity();
		$s_similar = $hsp->frac_conserved()*100;
                $subject_seq_length=$hsp->subject->length;
		
		$Q_AlignmentLength=$q_end-$q_beg+1; # The length of Query Alignment
		$S_AlignmentLength=$s_end-$s_beg+1; # The Length of Subject Alignment
		# deciding if we take the template
		if ($Q_AlignmentLength/$query_seq_length < $Min_Overlap_Percent) # Query Alignment length cover the query less than the minimal coverage
		{
			print "REJECT: ".join(" ",$s_name)." - Q_Overlap:".$Q_AlignmentLength."/".$query_seq_length." < ".$Min_Overlap_Percent."\n";
		}	
		elsif ($S_AlignmentLength/$subject_seq_length < $Min_Overlap_Percent) # Subject Alignment length cover the subject less than the minimal coverage
		{
			print "REJECT: ".join(" ",$s_name)."- S_Overlap:".$S_AlignmentLength."/".$subject_seq_length." < ".$Min_Overlap_Percent."\n";
		}	
		elsif ($s_ident<$Min_ID) # share less than the minimal identitiy
		{
			print "REJECT: ".join(" ",$s_name)." - ".$s_ident."<".$Min_ID."\% ID (Min:95\%)\n";
			
		}
		else # Take the template
		{
			# print "SELECT: $s_name - ".join(" ",$s_description)."\n";
			push (@$Selected_Templates_Nems_ref,join(" ",$s_name));
			$Selected_Templates_Details_ref->{$s_name}->{'DESCR'} = join (" ",$s_description);
			$Selected_Templates_Details_ref->{$s_name}->{'E_VAL'} = $s_eval;
			$Selected_Templates_Details_ref->{$s_name}->{'IDENTITY'} = $s_ident;
			$Selected_Templates_Details_ref->{$s_name}->{'SIMILARITY'} = $s_similar;
			$Selected_Templates_Details_ref->{$s_name}->{'SEQ'}=$AAseq;
			$Selected_Templates_Details_ref->{$s_name}->{'Q_AlignmentLength'}=$Q_AlignmentLength;
			$Selected_Templates_Details_ref->{$s_name}->{'S_AlignmentLength'}=$S_AlignmentLength;
			$Selected_Templates_Details_ref->{$s_name}->{'UNIPROT_ID'}=$s_accession;
		}
	   }
	 }
	}
     return ("OK",$Selected_Templates_Nems_ref)         
}
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}   
1;
