#!/usr/bin/perl


use lib "/bioseq/bioSequence_scripts_and_constants";
use GENERAL_CONSTANTS;
use strict;

#****************************************
#command line arg
my $Prog = shift; # the calc program: 'consurf' or 'rate4s'
my $mode = shift; # 'msa' or 'fasta'
my $SeqName = shift; # in case of 'msa'
   $SeqName = '' if ($SeqName == -1);
my $bayesInterval = shift;
my $insufficientDataBool = shift; # 0 or 1
my $num_of_visit = shift; # 1 or 2
my $out_error_file = shift;
my $OutHtmlFile = shift;
my $pdbInsertionFile = shift;
my $aln_file = shift;
my $pdbdata = shift;
my $prog_out = shift;
my $final_out = shift;
my $pdb_file = shift;
my $msa = shift;
my $fasta_seq = shift;
my $insertionsListFile = shift;
my $gradefreq= shift;
my $gradefreqIsd= shift;
my $insertionPE = shift;
my $consurfSptFile = shift;
my $rasmolFile = shift;
my $sptFile = shift;
my  $insertionsList = ',';

##### INSERTION #############
my  $insertionsList;
open (INSERTIONSLISTFILE, $insertionsListFile);
while(<INSERTIONSLISTFILE>){
    $insertionsList = $_;
}

close INSERTIONSLISTFILE;


#******************************
#variables
my @Output = ();  # To hold all the information that should be printed in the output file
                  # (SEQRES, ATOM, grade, 3LATOM, color, reliability) 

my @PdbSeq = (); #to hold the PDB aligned to SEQRES (without the gaps in the seqres, if any)
my @SeqPdb = (); # to hold the SEQRES aligned to the PDB
my @PDB = ();    #to hold the PDB AND numbers ready to PE
my @gradesPE_data = (); #To hold data of updated gradesPE
my $ConsColorUnity; #unity of conservation to be colored
my $MaxCons;        #MaxCons Value
my @ColorLayers;    #color limit values between the layers
#calculate the number of informative sequences in each position (MSA DATA)
&calc_reliability;
my @colorstep; #color steps
$colorstep[10] = "[255,255,150]";
$colorstep[8] = "[16,200,209]";    #less conserved
$colorstep[7] = "[140,255,255]";
$colorstep[6] = "[215,255,255]";
$colorstep[5] = "[234,255,255]";
$colorstep[4] = "[255,255,255]";
$colorstep[3] = "[252,237,244]";
$colorstep[2] = "[250,201,222]";
$colorstep[1] = "[240,125,171]";
$colorstep[0] = "[160,37,96]";      #most conserved

my %ColorScale = (0 => 9,
		  1 => 8,
		  2 => 7,
		  3 => 6,
		  4 => 5,
		  5 => 4,
		  6 => 3,
		  7 => 2,
		  8 => 1);  
########### insertion #########
my %insertionsAtomNum;


### read the pairwise alignment
&Read_aln_file;
### read the pdbdata file that contains the ATOM sequence in 3 letter code, the residue number
### and the chain (for example: ARG34:A) 
&Read_pdb_data;
### read the scores file (the output of rate4site or consurf.pl)
&Read_grades;  
#calculate color unit layer thickness and max cons value
($MaxCons, $ConsColorUnity) = &Calc_layer_unit;
# calculate 1-9 color grades
@ColorLayers = &calc_color_layers ($MaxCons, $ConsColorUnity);

# set the color grade for each residue based on the conservation score with bayes algorithm
if ($Prog eq 'bayesRate'){   
   &Set_colors_bayes_interval;
}
# set the color grade for each residue based on the conservation score
&Set_colors_pdb;

#calculate the number of informative sequences in each position (MSA DATA)
&calc_reliability;
#calc alternative amino-acids (RESIDUE VARIETY)
&calc_residue_range; 
#print output file for PE protein coloring
&Print_consurf_spt;
#print file for PE sequence coloring (ConSurf Seq3D)
#&Print_consurf_js;  :  REMARK: removed this routine, as it creates a javascript which doesn't seem to be important
#print the 'amino acid conservation scores' file
&Print_PE;
#Insert the grades into the tempFactor column of the PDB file
# (not in use right now because for some reason it doesn't work with grasp)
&insert_grades;

&create_file_cgfPE;
&create_file_insertionPE;

#*** SUBROUTINES ****************************
#####################################
####################################

############################################################################
# Read the PDB-SEQRES (or PDB-MSA) alignment into two arrays
sub Read_aln_file {

    unless (open ALN2, "<$aln_file"){ #from SEQRES-PDB clustalw alignment
        &error_and_exit("sys", "r4s_res_to_gradesPE : Read_aln_file: Cannot open the file $aln_file $!");
    } 
    
    my $pdbseq; #to hold the ATOM sequence
    my $seqpdb; #to hold the SEQRES sequence

    while (<ALN2>) {        
        if (/^PDB_\w*\s+(\S+)/) {	    
            $pdbseq .= $1;
        }
        elsif (/^SEQRES_\w*\s+(\S+)/){
            $seqpdb .= $1;
        }
    	# in case of msa
        elsif (/^\Q$SeqName\E\s+(\S+)/){
        $seqpdb .= $1; 
        }
    }
    close ALN2;
  
    my @PdbSeq_temp = split (//, $pdbseq); 
    @SeqPdb = split (//, $seqpdb);
    # copy @PdbSeq_temp to @PdbSeq, 
    # without the places in which there are gaps in the SEQRES (if any)
    # or 'X' (unknown) amino acids if it is an MSA 
    my $i = 0;
    foreach my $AA (@SeqPdb){
        if ($AA ne "-" and $AA ne "X"){            
            push @PdbSeq, $PdbSeq_temp[$i];
        }
        $i++;
    }
}

################################################################################
#open and read file pdbdata into an array
sub Read_pdb_data {

    unless (open PDB, "<$pdbdata"){ #from SEQRES-PDB clustalw alignment
        &error_and_exit("sys", "r4s_res_to_gradesPE : Read_pdb_data: Cannot open the file $pdbdata $!");
    }
    my @PDB_temp = <PDB>;
    chomp @PDB_temp;
    close PDB;

    # copy @PDB_temp to @PDB, 
    # without the places in which there are gaps in the SEQRES (if any)
    # or non-sandard amino acids if it is an MSA 
    my $i = 0;
    foreach my $AA (@SeqPdb){
        if ($AA ne "-"){
	    
	    $PDB_temp[$i] =~ s/\s*$//;  # to remove the spaces
            push @PDB, $PDB_temp[$i];
        }
        $i++;
    }

}

#########################################################################################
# Read the grades file and put all the information in @output in the correct alignment!
sub Read_grades {
    # open the output file of 'consurf' or 'rate4s' for reading
    unless (open GRADES, "<$prog_out"){
        &error_and_exit("sys", "r4s_res_to_gradesPE : Read_grades: Cannot open the file $prog_out $!");
    }

    #open the output file of 'consurf insertions with atom number' for reading
    unless (open PDBDATAINSERTIONS, "<$pdbInsertionFile"){
        &error_and_exit("sys", "r4s_res_to_gradesPE : Read_grades: Cannot open the file $pdbInsertionFile $!");
    }

    while (<PDBDATAINSERTIONS>){
        if (/^\s*(.*)\s*\,\s*(.*)\s*,\s*(.*)\s*/)
        {
            my $key =$1;
            my $valueAtom = $2;
            my $valueAtomNo = $3;
            $key =~ s/\s*$//;
            my $atomNumber = $valueAtom;
            $atomNumber =~ s/\s*[A-Z]+(\d+)[A-Z]?\s*/$1/;
            if ($insertionsList =~ /\,$atomNumber\,/){
                $insertionsAtomNum{$key} = $valueAtomNo;
            }
            else{
                $insertionsAtomNum{$key} = $valueAtom;
            }
        }
    }

    # skip the title
    if ($Prog eq "consurf"){
        $_ = <GRADES>;
    }
    my $CountPdb = 0; # for the 3LATOM
    my $CountSeq = 0; # for the SEQ

    while (<GRADES>) {
        my $line = $_;
        chomp $line;
        if ($line =~ /^\s*(\d+)\s+(\S+)\s+(\S+)/){    
            my $AA = $2;
            my $grade = $3;
            # ignore lines with '*' - meaning gaps between the seqres and its homologues
            next if ($AA eq "*");
            $Output[$CountSeq]{SEQ} = $AA;
            $Output[$CountSeq]{GRADE} = $grade;
            $Output[$CountSeq]{ATOM} = $PdbSeq[$CountSeq];
            if ($PdbSeq[$CountSeq] ne "-") {
                $Output[$CountSeq]{PE} = $PDB[$CountPdb];
                $CountPdb++;
            }
            else {
                $Output[$CountSeq]{PE} = "-";
            }
            if ($line =~ /^\s*\d+\s+\S+\s+\S+\s+\[\s*([0-9\.\-]+)\,\s*([0-9\.\-]+)\]/){    
                $Output[$CountSeq]{INTERVALLOW}= $1;
                $Output[$CountSeq]{INTERVALHIGH}= $2;
             }
            $CountSeq++;    
        }
    }
    close GRADES;
}

#########################################################################
# calc the thickness of the color layers and the most conserved value
sub Calc_layer_unit {

    my $element;
    my $max_cons = $Output[0]{GRADE};
    my $ConsColorUnity; #unity of conservation to be colored
    foreach $element (@Output){
        if ($$element{GRADE} < $max_cons) {
            $max_cons = $$element{GRADE};
        }
    }
    $ConsColorUnity = $max_cons / 4.5 * -1; 
    if ($max_cons !~ /^\-/){
        $ConsColorUnity = $max_cons;      
    }
    return ($max_cons, $ConsColorUnity);
}

#################################################################
# calc the 1 - 9 color grades
sub calc_color_layers ($$) {

    my $MaxCons = $_[0];
    my $ConsColorUnity = $_[1];
    my $i;
    my $NoLayers = 9;
    my @ColorLayers;

    for ($i = 0; $i <= $NoLayers; $i++) {

	$ColorLayers[$i] = $MaxCons + ($i * $ConsColorUnity);
    }

    return @ColorLayers;

}

#############################################################
# Description : Calculates and finds the sites that have less
#               informative in the Confidence interval (Bayes).

sub Set_colors_bayes_interval
{
    my $i;
    my $element;
    my $Count = 0;
    foreach $element (@Output) {
        
        for ($i = 0; $i <= $#ColorLayers; $i++) {
            if ( ($i==$#ColorLayers) and !exists $$element{INTERVALLOWCOLOR}){
                $$element{INTERVALLOWCOLOR} = 8;             
            }
			elsif ($$element{INTERVALLOW} >= $ColorLayers[$i] and $$element{INTERVALLOW} < $ColorLayers[$i + 1]) {
                $$element{INTERVALLOWCOLOR} =$i;         
            } 
            elsif ( ($$element{INTERVALLOW} < $ColorLayers[$i]) and !exists $$element{INTERVALLOWCOLOR}){
                $$element{INTERVALLOWCOLOR} = 0;         
            }            
			if (($i == $#ColorLayers)  and !exists $$element{INTERVALHIGHCOLOR}){            
                $$element{INTERVALHIGHCOLOR} = 8;            
            } 
            elsif ($$element{INTERVALHIGH} >= $ColorLayers[$i] and $$element{INTERVALHIGH} < $ColorLayers[$i + 1]) {
                $$element{INTERVALHIGHCOLOR} =$i;
            }            
            elsif ( ($$element{INTERVALHIGH} < $ColorLayers[$i]) and ($$element{INTERVALHIGHCOLOR} eq '')){
                $$element{INTERVALHIGHCOLOR} = 0;
            } 
        }
    }
}

#Set the colors and the color grades for each AA
sub Set_colors_pdb {
    
    my $i;
    my $element;
    my $Count = 0;
   foreach $element (@Output)
    {
        for ($i = 0; $i <= $#ColorLayers; $i++) {
            if ($$element{GRADE} >= $ColorLayers[$i] && $$element{GRADE} < $ColorLayers[$i + 1]) {	
                $Output[$Count]{COLOR} = $colorstep[$i]; #color for chime
                $Output[$Count]{SEQ3DCOLOR} = $i;          #number of color for seq3d
                ($Output[$Count]{NUMBER}) = ($$element{PE}) =~ /(\d+)/; #res number
                $$element{reliability} =~ /(.*)\/(.*)/;
                my $msaDataNumerator = $1;
                my $msaDataDenominator = $2;				
                if ($Prog eq 'bayesRate'){					
                    $Output[$Count]{INTERVALLOWSEQ3DCOLOR} = $$element{INTERVALLOWCOLOR};
        		    $Output[$Count]{INTERVALHIGHSEQ3DCOLOR} = $$element{INTERVALHIGHCOLOR};
                    if (((($$element{INTERVALHIGHCOLOR} -  $$element{INTERVALLOWCOLOR}) > $bayesInterval ) or ($msaDataNumerator <= 5)) and ( $insufficientDataBool == 1 )){
                       $Output[$Count]{COLOR} = $colorstep[10];
                       $Output[$Count]{SEQ3DCOLOR} = 10;
                       $Output[$Count]{REALSEQ3DCOLOR} = $i;
                    }                  
                }
                else{               
                    if ( $msaDataDenominator != 0){	
                        if (($msaDataNumerator <= 5) and ( $insufficientDataBool == 1)){
                            $Output[$Count]{COLOR} = $colorstep[10];
                            $Output[$Count]{SEQ3DCOLOR} = 10;
                            $Output[$Count]{REALSEQ3DCOLOR} = $i; 
                        }
                    }
                }
                $Count++;
                last;
            }
            elsif ($i == 9) {
                $Output[$Count]{COLOR} = $colorstep[8];
                $Output[$Count]{SEQ3DCOLOR} = 8;          #color for seq3d
                ($Output[$Count]{NUMBER}) = ($$element{PE}) =~ /(\d+)/;
               
                $$element{reliability} =~ /(.*)\/(.*)/;
                my $msaDataNumerator = $1;
                my $msaDataDenominator = $2;  
                if ($Prog eq 'bayesRate'){
                    $Output[$Count]{INTERVALLOWSEQ3DCOLOR} = $$element{INTERVALLOWCOLOR};
                    $Output[$Count]{INTERVALHIGHSEQ3DCOLOR} = $$element{INTERVALHIGHCOLOR};
                   if (((($$element{INTERVALHIGHCOLOR} - $$element{INTERVALLOWCOLOR}) > $bayesInterval )  or ($msaDataNumerator <= 5)) and  ($insufficientDataBool == 1)){
                       $Output[$Count]{COLOR} = $colorstep[10];
                       $Output[$Count]{SEQ3DCOLOR} = 10; 
                       $Output[$Count]{REALSEQ3DCOLOR} = 8;
                   }                  
                }
                else{               
                    if (($msaDataNumerator <= 5) and ($insufficientDataBool == 1)){ 
                        $Output[$Count]{COLOR} = $colorstep[10];
                        $Output[$Count]{SEQ3DCOLOR} = 10;
                        $Output[$Count]{REALSEQ3DCOLOR} = 8;  
                    }
                }
                $Count++;
            }###elsif
        }###for
    }###foreach
}### sub


##########################################################################
# Print consurf.spt file (RasMol coloring script) for PE
sub Print_consurf_spt {

    my $i;
    my $element;
    my $color;
    my $LineMaxElem = 10;  #No of elements on each consurf.spt line
    my $ElemCount = 0;
    my $separator;
    my $chain = '';
    #my $chain;
    my $displayRsmlElem;
    my $insertionSimilarAAstr;

    unless (open SPT, ">$consurfSptFile"){
        &error_and_exit("sys", "r4s_res_to_gradesPE : Print_consurf_spt: Cannot open the file $consurfSptFile $!");
    }

#============================ RasMol (CON) =============
    unless (open RSML, ">$rasmolFile"){
        &error_and_exit("sys", "r4s_res_to_gradesPE : Print_consurf_spt: Cannot open the file $rasmolFile for writing $!");
    }
#========================================================

    # to solve the problem when chime runs out of colors.
    print SPT "select all \n";
    print SPT "color [200,200,200]\n\n";

#============================ RasMol (CON) =============
    print RSML "select all\n";
#=======================================================
#============================ RasMol (CON) =============
    print RSML "color [200,200,200]\n\n";
#=======================================================
    my $file_count = 9;
	mkdir "pdbspt";
    foreach $color (@colorstep) {
        my $color_elem_count = 0;
        unless (open COLOR, ">pdbspt/" . $file_count . $sptFile){
            &error_and_exit("sys", "r4s_res_to_gradesPE : Print_consurf_spt: Cannot open the file /pdbspt/" . $file_count . $sptFile." for writing: $!");
        }

        foreach $element (@Output) { 
            if (  $$element{PE} =~ /\s*.*\:([A-Z])\s*/ ){ 
                $chain = $1;
            }
            if (($$element{COLOR} eq $color) and ($$element{PE} ne "-" )) {
                if ($insertionsAtomNum{$$element{PE}} =~ /atomno/){
                    $insertionSimilarAAstr .= "\nselect selected or $insertionsAtomNum{$$element{PE}}";
                }
                else
                {
                    if ($ElemCount == 0){# and ($insertionsAtomNum{$$element{PE}} !~ /atomno/)) {
                        print SPT "\nselect  ";
                        print COLOR "\nselect  ";
#============================ RasMol (CON) =============
                        print RSML "\nselect  ";
#=======================================================		   
                    }
                    if ($ElemCount == 10){
                        print SPT "\nselect selected or $insertionsAtomNum{$$element{PE}}";
            		    print COLOR "\nselect selected or $insertionsAtomNum{$$element{PE}}";

#============================ RasMol (CON) =============
                        $displayRsmlElem = $$element{PE}; 
                        if ($chain =~ /[A-Z]/){
                            $displayRsmlElem =~ s/\:$chain//;
                        }
                        print RSML "\nselect selected or $insertionsAtomNum{$$element{PE}}";
#======================================================
		                $ElemCount = 0;
                    }
                    else {
                        print SPT "$separator $insertionsAtomNum{$$element{PE}}";
                        print COLOR "$separator $insertionsAtomNum{$$element{PE}}";
#============================ RasMol (CON) =============
                        $displayRsmlElem = $$element{PE};
                        if ($chain =~ /[A-Z]/){
                            $displayRsmlElem =~ s/\:$chain//;
                        }
                        print RSML "$separator $insertionsAtomNum{$$element{PE}}";
#======================================================
                    $separator = ",";
                    }
                    $ElemCount++;
                    $color_elem_count++;
                }
            }
        }
        if ($ElemCount > 0) {         
            print SPT  $insertionSimilarAAstr;
            print RSML $insertionSimilarAAstr;
            print COLOR $insertionSimilarAAstr;
            $insertionSimilarAAstr = '';
    
            print SPT "\nselect selected and:$chain\n";
            print SPT "\ncolor $color\nspacefill\n";
               #============================ RasMol (CON) =============
            if ($chain =~ /[A-Z]/){
               print RSML "\nselect selected and :$chain\n";
            }
            else{
                print RSML "\nselect selected and:\n";
            }
            print RSML "color $color\nspacefill\n";
            print RSML "\ndefine CON" . $file_count . " selected\n\n";
#=======================================================
            print COLOR "\nspacefill\n";
        }
        if ($color_elem_count == 0){
            print COLOR "javascript alert(\"No residues have this color\")";
        }
        $ElemCount = 0;
        $separator = " ";
        $file_count--;     
        if ($file_count == -1){
            $file_count = 10;
        }
        close COLOR;
    }    
    close SPT;
}

#################################################################################
# print consurf.js file for PE Seq3D
sub Print_consurf_js {

    my $i;
    my $element;
    my $Count = 0;
    my $consurfJSfile;
    my $total = scalar (@Output);

    if ( $insufficientDataBool == 1 ){
        $consurfJSfile = 'consurfisd.js';
    }
    else{
        $consurfJSfile = 'consurf.js';
    }   
    unless (open JS, ">$consurfJSfile"){
        &error_and_exit("sys", "r4s_res_to_gradesPE : Print_consurf_spt: Cannot open the file $consurfJSfile for writing$!");
    }
    print (JS "var csc = new Array(", $total +1 , ");\n\n");
    print (JS "for (i = 0; i <= " ,$total  , "; i++)\n");
    print JS "\tcsc[i] = -1;\n\n";
    foreach $element (@Output) {         
        if ($$element{PE} ne "-" ) {            
            print (JS "csc[" , $$element{NUMBER}, "] = ", $$element{SEQ3DCOLOR}, ";\n");       
        }
    }
}

############################################################################
# Print the 'amino acid conservation scores' file
sub Print_PE {
    
    # open the final output file for writing
    unless (open PE, ">$final_out"){
        &error_and_exit("sys", "r4s_res_to_gradesPE : Print_consurf_spt: Cannot open the file $final_out for writing $!");
    }
    print PE "\t Amino Acid Conservation Scores\n";
    print PE "\t===============================\n\n";
    print PE "- POS: The position of the AA in the SEQRES derived sequence.\n";
    print PE "- SEQ: The SEQRES derived sequence in one letter code.\n";
    print PE "- 3LATOM: The ATOM derived sequence in three letter code, including the AA's positions as they appear in the PDB file and the chain identifier.\n";
    print PE "- SCORE: The normalized conservation scores.\n";
    print PE "- COLOR: The color scale representing the conservation scores (9 - conserved, 1 - variable).\n";
    if ($prog_out eq 'bayesR4s.res'){
	print PE "- CONFIDENCE INTERVAL: When using the bayesian method for calculating rates, a confidence interval is assigned to each of the inferred evolutionary conservation scores.\n"; 
	print PE "- CONFIDENCE INTERVAL COLORS: When using the bayesian method for calculating rates. The color scale representing the lower and upper bounds of the confidence interval.\n"; 
    }
    print PE "- MSA DATA: The number of aligned sequences having an amino acid (non-gapped) from the overall number of sequences at each position.\n";
    print PE "- RESIDUE VARIETY: The residues variety at each position of the multiple sequence alignment.\n\n";

    if ($prog_out eq 'bayesR4s.res')
    {
        print PE " POS\t SEQ\t    3LATOM\tSCORE\t\tCOLOR\tCONFIDENCE INTERVAL\tCONFIDENCE INTERVAL COLORS\tMSA DATA\tRESIDUE VARIETY\n";
        print PE "    \t    \t        \t(normalized)\t        \t               \n";
    }     
    else
    {
        print PE " POS\t SEQ\t    3LATOM\tSCORE\t\tCOLOR\tMSA DATA\tRESIDUE VARIETY\n";
        print PE "    \t    \t        \t(normalized)\t        \t               \n";
    }

   my $pos = 1;
    foreach my $elem (@Output){

	printf (PE "%4d", "$pos");
	printf (PE "\t%4s", "$$elem{SEQ}");
	my $displayElement = $$elem{PE};
	printf (PE "\t%10s", "$displayElement");	  
	printf (PE "\t%6.3f", "$$elem{GRADE}");
	# in case it is undefiend: we add a '*' after the grade
          if  ($$elem{SEQ3DCOLOR} eq '10')
          {
            printf (PE "\t\t%3d", "$ColorScale{$$elem{REALSEQ3DCOLOR}}");
            printf (PE "%1s", "*");       
          }
          else
          {
	        printf (PE "\t\t%3d", "$ColorScale{$$elem{SEQ3DCOLOR}}");
          }
            if ($prog_out eq 'bayesR4s.res')
            {
              printf (PE "\t%6.3f", "$$elem{INTERVALLOW}");
              printf (PE "%1s", ",");
              printf (PE "%6.3f", "$$elem{INTERVALHIGH}"); 
			  printf (PE "\t\t\t%5d", "$ColorScale{$$elem{INTERVALLOWSEQ3DCOLOR}}");
              printf (PE "%1s", ","); 
              printf (PE "%1d\t\t", "$ColorScale{$$elem{INTERVALHIGHSEQ3DCOLOR}}");
              my $baysColorInterval =   $ColorScale{$$elem{INTERVALLOWSEQ3DCOLOR}} - $ColorScale{$$elem{INTERVALHIGHSEQ3DCOLOR}};
            } 
        printf (PE "\t%8s", "$$elem{reliability}");
        printf (PE "\t%-18s\n", "$$elem{res_range}"); # to align left

	$pos++;
    }
     if ($prog_out eq 'bayesR4s.res')
            {
               	print PE "\n\n*Below the confidence cut-off - The calculations for this site were performed on less than 6 non-gaped homologue sequences,\nor the confidence interval for the estimated score is equal to- or larger than- 4 color grades.\n";
            } 
            else
            {
	
		print PE "\n\n *Below the confidence cut-off - The calculations for this site were performed on less than 6 non-gaped homologue sequences. \n";
            }
    close PE;
}

###########################################################################
# Insert the scores into the Tempfactor column in the PDB file
sub insert_grades{

    my %grades = ();
    my $pe;
    my $chain;
    my $resnum;
    my $grade;

    # read the 3LATOM and GRADE fields from the output file
    open PE, "<$final_out";
    while (<PE>){

	if ($_ =~ /\s*\d+\s+\w+\s+(\S+)\s+(\S+)\s+\d+/){

	    $pe = $1;
	    $grade = $2;

	    # ignore gaps in the ATOM sequence
	    unless ($pe eq "-"){

		$pe =~ /\w\w\w(\S+):(\w)?/;
		$resnum = $1;
		$chain = $2;

		if ($chain eq ""){
		    
		    $chain = " ";
		}

		unless ($grade =~ /^-/){

		    # add a space so that the grades will be indented to the right
		    $grade = " " . $grade;
		}

		$grades{$resnum} = $grade;
	    }
	}
    }
    close PE;
    
    # call the function insert_Tempfactor to insert the grades to the PDB file
    # if temprature factor field was not found in some of the lines - report it to the user
    my $ans = &insert_tempFactor($pdb_file, $chain, \%grades);
    if ($ans eq "no_TF" && $num_of_visit==1){
        open OUTPUT, ">>$OutHtmlFile";
        print OUTPUT "\n<p><ul><li><font color='red'><b>Warning:</b></font> The \'Temperature Factor\' field was not found in some of your \'ATOM\' fields. Therefore it will not be fully updated with ConSurf grades.</li></ul></p>\n";
        close OUTPUT;
    }
}

#######################################################################
# calc the number of informative data per each position
# Description : This procedure Calculates and finds the sites that have less
#               informative data min 6 sites or 10% (ML).

sub calc_reliability{

    my @aln = ();

    &read_aln(\@aln); 
    
    # extract the target sequence name from the ".seq" file
    if ($mode eq "fasta"){
		open SEQ, "<$fasta_seq";
		while (<SEQ>){
			if ($_ =~ />(\S+)\s+/){
				$SeqName = $1;
			}
		}
		close SEQ;
    }
    
    # loop over the positions
    my $CountPos = 0;
    foreach my $pos (@aln){
		# calculate the reliability only for positions in which the target sequence is not "-" or "X"
		if ($$pos{$SeqName} ne "-" and $$pos{$SeqName} ne "X"){
			my $total = 0;
			my $inf = 0;
			# loop over the different values
			foreach my $aa (values %{$pos}){
				unless ($aa eq "-"){				
					$inf++;
				}
				$total++;
			}
			
			# insert the reliability field to the general array
			my $result = $inf . "/" . $total;
			$Output[$CountPos]{reliability} = $result;
			$CountPos++;
		}
    }
}

#####################################################################
# read the aln file into an array of hashes
# (for each position there is a hash, where the keys are the sequences names
# and the values are the residues)
#####################################################################
sub read_aln {

    my $aln = shift;

    # read the MSA file into array of hashes
    unless (open MSA, "<$msa"){
        &error_and_exit("sys", "r4s_res_to_gradesPE : read_aln: Cannot open the file $msa $!");
    }
    my $PosCount = 0;
    my $first_name = "";
    my $i;
    while (<MSA>){

	if ($_ =~ /^CLUSTAL/){
	    next;
	}

	elsif ($_ =~ /^(\S+)\s+(\S+)\s*/){

	    my $name = $1;
	    my $seq = $2;
	    my @seqs = split '', $seq;

	    # to set the first name
	    if ($first_name eq ""){

		$first_name = $name;
	    }
	    # a new block
	    elsif ($name eq $first_name){

		$PosCount = $i; 
	    }
	    

	    $i = $PosCount;
	    foreach my $cell (@seqs){

		$$aln[$i]{$name} = $cell;
		$i++;
	    }
	}
    }
    close MSA;
}

#####################################################################
# calc the residues variety per each position
sub calc_residue_range {

    my @aln = ();
    
    &read_aln(\@aln); 

     # extract the target sequence name from the ".seq" file
    if ($mode eq "fasta"){

	open SEQ, "<$fasta_seq";
	while (<SEQ>){

	    if ($_ =~ />(\S+)\s+/){

		$SeqName = $1;
	    }
	}
	close SEQ;
    }
    
    # loop over the positions
    my $CountPos = 0;
    foreach my $pos (@aln){

	# find the residues range only for positions in which the target sequence is not "-" or "X"
	if ($$pos{$SeqName} ne "-" and $$pos{$SeqName} ne "X"){

	    # sort the residues in alphabeical order
	    my @sorted = ();
	    push @sorted, (sort { $a cmp $b } values %{$pos});

	    $Output[$CountPos]{res_range} = "";

	    # loop over the different values
	    foreach my $aa (@sorted){

		# write the AA if it's not a gap and it hasn't appeared yet
		if ($aa ne "-" and $Output[$CountPos]{res_range} !~ /$aa/){

		    $Output[$CountPos]{res_range} .= "$aa,";
		}
	    }
	    
	    # delete the last comma
	    chop $Output[$CountPos]{res_range};
	    
	    $CountPos++;
	}
    }
}

##################### build a file in order to create consurf_grade_freqs for PE version 2.727

sub create_file_cgfPE{
    my %gradefreqHash;
    my %gradefreqIsdHash;  

    unless (open GRADESFREQ, ">$gradefreq"){
        &error_and_exit("sys", "r4s_res_to_gradesPE : create_file_cgfPE: Cannot open the file $gradefreq for writing$!");
    }
    unless (open GRADESFREQISD, ">$gradefreqIsd"){
        &error_and_exit("sys", "r4s_res_to_gradesPE : create_file_cgfPE: Cannot open the file $gradefreqIsd for writing$!");
    } 

 print  GRADESFREQISD "seq3d_grades_isd=";  
 print  GRADESFREQ     "seq3d_grades=";
 my $keyIsd; 
 my $key;

  $gradefreqHash{0}=0; 

    my $firstTimeBool = '1';
    my $counter = 1;
    my $counterIsd = 1; 
    my $delimiter = 70;
    my $lastResidueNumber;
    my $insertionCount = 0;
    my $withoutCoordinates =0;
    foreach my $elem (@Output){
	$key = $ColorScale{$$elem{SEQ3DCOLOR}}; 
        $keyIsd =  $ColorScale{$$elem{SEQ3DCOLOR}};
##########################################
	if (($firstTimeBool eq '1') && ($$elem{PE} =~ /(\d+)/)){
	    if ($$elem{PE} =~ /(\d+)/){
	       if ($1 > 1){
		  my $differ = $1 - 1;
		  for (my $i=1; $i<=$differ; $i++){
	
	               print  GRADESFREQISD '.'; 	 
		       print  GRADESFREQ '.';
		        	&print_delimiters_into_gradesfreq($counter,$delimiter,'0');
		       &print_delimiters_into_gradesfreq($counterIsd,$delimiter,'1');
		       $counter++;
		       $counterIsd++;
		       $withoutCoordinates++;
		      
		   }
	       }
	    }
	    $firstTimeBool = '0';
	}    
########################################################


    if ($$elem{PE} eq '-')
    {
	if ($lastResidueNumber eq ''){
	
	               print  GRADESFREQISD '.'; 	 
		       print  GRADESFREQ '.';
		       &print_delimiters_into_gradesfreq($counter,$delimiter,'0');
	&print_delimiters_into_gradesfreq($counterIsd,$delimiter,'1');
		       $counter++;
		       $counterIsd++;
		       $withoutCoordinates++;
         }
    }
    elsif ($$elem{PE} =~ /\d+\s*[A-Z]?\s*:?\s*/)
    {
     if ($$elem{PE} =~ /(\d+)\s*[A-Z]\s*:\s*/)
     {
	 $insertionCount++;
     }
	   $$elem{PE} =~ /(\d+)/;
	   my $numr=$1;
	   if ($lastResidueNumber ne ''){	  
	     my $dif = $numr - $lastResidueNumber;	   
	       if  ($dif>1){  
		   for (my $i = 1 ; $i<$dif; $i++){
		     
			    print  GRADESFREQISD '.'; 	 
			    print  GRADESFREQ '.';
			    &print_delimiters_into_gradesfreq($counter,$delimiter,'0');
                &print_delimiters_into_gradesfreq($counterIsd,$delimiter,'1');
                $counter++;
		       $counterIsd++;
		       $withoutCoordinates++;
		   }
	       }
	   }
	   $lastResidueNumber = $numr;

	   if  ($$elem{SEQ3DCOLOR} eq '10')
	   {
	        if ($$elem{PE} !~ /\d+\s*[A-Z]\s*:?\s*/)
			{			   
			    print  GRADESFREQISD '0';
			    print  GRADESFREQ $ColorScale{$$elem{REALSEQ3DCOLOR}};
			    &print_delimiters_into_gradesfreq($counter,$delimiter,'0');
                &print_delimiters_into_gradesfreq($counterIsd,$delimiter,'1');
			    $counter++;
			    $counterIsd++;
			}
	        $gradefreqIsdHash{0}++;
	        $gradefreqHash{$ColorScale{$$elem{REALSEQ3DCOLOR}}}++;	   
	   }   
	   else
	   {
	        if ($$elem{PE} !~ /\d+\s*[A-Z]\s*:?\s*/)
			{
			   
			    print  GRADESFREQISD $ColorScale{$$elem{SEQ3DCOLOR}}; 
			    print  GRADESFREQ $ColorScale{$$elem{SEQ3DCOLOR}};
                &print_delimiters_into_gradesfreq($counter,$delimiter,'0');
                &print_delimiters_into_gradesfreq($counterIsd,$delimiter,'1');
			    $counter++;
			    $counterIsd++; 
	        }
	       $gradefreqIsdHash{$ColorScale{$$elem{SEQ3DCOLOR}}}++;
	       $gradefreqHash{$ColorScale{$$elem{SEQ3DCOLOR}}}++;	    
	   }
    }
    }
    print GRADESFREQ "\n";
     print  GRADESFREQISD "\n"; 
  #  foreach $a (sort keys %gradefreqHash){
    my $sum;
    for (my $i=0;$i<10;$i++){
	
	if ($gradefreqHash{$i}){
	    print GRADESFREQ "$i,$gradefreqHash{$i}\n";
	    $sum = $sum + $gradefreqHash{$i};
	}
    else
	{
	    print GRADESFREQ "$i,0\n";
	    $sum = $sum + $gradefreqHash{$i};
	}    
    }
     my $sumIsd;
  #    foreach $a (sort keys %gradefreqIsdHash){
    for (my $i=0;$i<10;$i++){
	if ($gradefreqIsdHash{$i}){
	    print GRADESFREQISD "$i,$gradefreqIsdHash{$i}\n";
	    $sumIsd = $sumIsd + $gradefreqIsdHash{$i};

        }else
        {
	    print GRADESFREQISD "$i,0\n";
	     $sumIsd = $sumIsd + $gradefreqIsdHash{$i};
	} 
    }
 close GRADESFREQ;
 close GRADESFREQISD;
}
sub print_delimiters_into_gradesfreq{
    my $counterFreq = shift;
    my $delimiterFreq = shift;
   
    my $isd = shift;
    
    my $divFreq = $counterFreq / $delimiterFreq;

    	if (($divFreq !~ /\./)  && ($counterFreq ne '0'))# && ($counter ne '0'))
	{

	    if ($isd eq '0'){
	        print  GRADESFREQ   '"+"';
            }
   	elsif ($isd eq '1')
	 {
	   print GRADESFREQISD '"+"';
          }
	}
}


#########$$$$$$$$$ insertions file for pe 2.727 ##############
sub create_file_insertionPE{
  my %insertionsHash;

   unless (open INSERTIONPE, ">$insertionPE"){
    &error_and_exit("sys", "r4s_res_to_gradesPE : create_file_insertionPE: Cannot open the file $insertionPE for writing$!");
    }
  my $keyIsd; 
  my $key;
  my $CountPos;
  my $lastPEelemColorIsd;
  my $lastPEelemColor;
  my $firstTimeBool=1;

 foreach my $elem (@Output)
 {    
    if ($$elem{PE} =~ /\-?(\d+)\s*[A-Z]?\s*:\s*/)
    {
		my $insertPos = $1;
		if ($insertionsList =~ /\,$1\,/){
		   if  ($$elem{SEQ3DCOLOR} eq '10' )
		   {
			   $insertionsHash{$insertPos}{colors_isd} .= '0'; 
			   $insertionsHash{$insertPos}{colors} .= $ColorScale{$$elem{REALSEQ3DCOLOR}};
		   }   
		   else{
				$insertionsHash{$insertPos}{colors_isd} .= $ColorScale{$$elem{SEQ3DCOLOR}}; 
				$insertionsHash{$insertPos}{colors} .= $ColorScale{$$elem{SEQ3DCOLOR}};
		   }   
	   }
    }
 }
    foreach my $a (keys %insertionsHash){         
        my $line = "$a,\"$insertionsHash{$a}{colors_isd}\",\"$insertionsHash{$a}{colors}\"";
        $line =~ s/\"\"//;
	    $line =~ s/\,\,/,/;
        print  INSERTIONPE "$line\n";
    }
}
########################################
sub error_and_exit{
    my $h_err = shift;
    my $l_err = shift;
    open ERROR_OUT, ">".$out_error_file;
    print ERROR_OUT "HTML: $h_err\n";
    print ERROR_OUT "LOG: $l_err";
    close ERROR_OUT;
    chmod 0755, $out_error_file;
    exit;
}

#############################################################
# Description:  replace the tempFactor column in the PDB file
#               with another column. (for example: grades)
# Arguments: 1. PDB file
#            2. chain
#            3. reference to hash
#               key - the residue sequence number
#               value - the field to insert instead of the tempFactor
#############################################################
sub insert_tempFactor {

    my $PdbFile = shift;
    my $chain = shift;
    my $HashRef = shift;

    my $PDBchain;
    my $ResNum;
    my $iCode;
    my $ret = "nothing"; 
    
    # read the PDB file to an array
    open READPDB, "<$PdbFile";
    my @PDB = <READPDB>;
    close READPDB;

    # write the PDB file and replace the tempFactor column
    # with the new one.
    open WRITEPDB, ">$PdbFile";
    
    foreach my $line (@PDB) {
	if ($line =~ /^ATOM/){	 
	    $PDBchain = substr ($line, 21, 1);
            $ResNum = substr ($line, 22, 4);
	    $ResNum =~ s/\s*//g;
	    $iCode = substr ($line, 26, 1);
            if ($iCode =~ /^\s*$/){
                $iCode = '';
            }
            if ($PDBchain eq $chain){
                my $ResNumiCode=$ResNum . $iCode;
                $$HashRef{$ResNumiCode} =~ s/\s*$//;
                if (length($line)>60 && length($line)>=66)# to make sure that the temprature factor is written in the line
                    {substr ($line, 60, 6) = $$HashRef{$ResNumiCode};}
                else{
                    $ret = "no_TF";
                }
                print WRITEPDB $line;
            }
            else {
                print WRITEPDB $line;
            }
	}
	else {
	    print WRITEPDB $line;
	}
    }
    close WRITEPDB;
    return $ret;
}