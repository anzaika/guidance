#!/usr/local/bin/perl


#*** Modules ***#
#use lib '/data/users/consurf/WebConsurf/modules';  
#use PDB;
use Bio::SeqIO;
use strict;
use lib "/bioseq/bioSequence_scripts_and_constants"; 
use GENERAL_CONSTANTS;

#die "Usage: colorCoding.v2.pl pdbPrefix prog(ML or Bayesian) run_name referenceSeqName\n" if (scalar(@ARGV) < 4);

#****************************************
#command line arg
#my $pdbPrefix = $ARGV[0];  # The prefix of the file(pdb id)
my $Prog = shift;#$ARGV[1]; # the calc program: 'ML' or 'Bayesian'
my $SeqName = shift;#$ARGV[3]; 

# files names
my $WorkingDir = shift;
my $file_colorRes = shift; #$WorkingDir."selection4Site.txt";#file with list of residues and their color range
my $aln_file = shift; #$pdbPrefix."_PDB_MSA.aln";
my $pdbdata = shift; #$pdbPrefix . ".pdbdata";
my $prog_out = shift; #"kaks.res"; #sites and ka/ks values
my $final_out = shift; #$pdbPrefix . ".gradesPE"; #for protein explorer
my $pdb_file = shift; #$pdbPrefix . ".ent";
my $msa = shift; #$pdbPrefix ."AMINO". ".msa"; #amino msa file
#my $spt_file = "consurf.spt";  # REMARK: I DELETED ALL THE PRINTS TO THIS FILE, I DON'T THINK IT IS NESECERAY
my $rasmol_file = shift; #"consurf.rsml";
my $global_res = shift; #$WorkingDir."globalResult.txt";

#******************************
#variables

my @Output = ();  # To hold all the information that should be printed in the output file
# (SEQRES, ATOM, grade, pvalue, 3LATOM, color, reliability) 
my @Aln = ();		  # Holds all the aligment information: each position and its content
my @PdbSeq = (); #to hold the PDB aligned to SEQRES (without the gaps in the seqres, if any)
my @SeqPdb = (); # to hold the SEQRES aligned to the PDB
my @PDB = ();    #to hold the PDB AND numbers ready to PE
my @gradesPE_data = (); #To hold data of updated gradesPE
my $ConsColorUnity; #unity of conservation to be colored
my $MaxCons;        #MaxCons Value
my @ColorLayers;    #color limit values between the layers
my ($html_error, $log_error); # to hold the error lines that will be written in case an error was found

my @colorstep; #color steps
$colorstep[6] = "[255,190,0]";      # --> Ka/Ks>1 significant
$colorstep[5] = "[255,255,120]";    # --> Ka/Ks>1
$colorstep[4] = "[255,255,255]";    # --> Ka/Ks<1
$colorstep[3] = "[252,237,244]";    # --> Ka/Ks<1
$colorstep[2] = "[250,201,222]";    # --> Ka/Ks<1
$colorstep[1] = "[240,125,171]";    # --> Ka/Ks<1
$colorstep[0] = "[160,37,96]";      # --> Ka/Ks<1 significant


my %ColorScale = (0 => 7,
		  1 => 6,
		  2 => 5,
		  3 => 4,
		  4 => 3,
		  5 => 2,
		  6 => 1,
		  );

my %ColorScale_reverse = (7 => 0,
			  6 => 1,
			  5 => 2,
			  4 => 3,
			  3 => 4,
			  2 => 5,
			  1 => 6,
			  );

### read the pairwise alignment
&Read_aln_file;

### read the amino acid MSA and create a hash with all the relevant information
&read_aln();

### read the pdbdata file that contains the ATOM sequence in 3 letter code, the residue number
### and the chain (for example: ARG34:A) 
&Read_pdb_data;

### read the scores file (the output of kaks4site: kaks.res)
&Read_grades;  

### read the selection4site: color value for each residue
&Read_res_to_color;

#calculate the number of informative sequences in each position (MSA DATA)
#&calc_reliability;

#calc alternative amino-acids (RESIDUE VARIETY)
#&calc_residue_range; 

#print output file for PE protein coloring
&Print_consurf_spt;

#print file for PE sequence coloring (Seq3D)
#&Print_consurf_js;

#print the 'amino acid conservation scores' file
&Print_PE;

#Insert the grades into the tempFactor column of the PDB file
# (not in use right now because for some reason it doesn't work with grasp)
&insert_grades;



#*** SUBROUTINES ****************************
#####################################
####################################

############################################################################
# Read the PDB-SEQRES (or PDB-MSA) alignment into two arrays
sub Read_aln_file {

    my $path_to_aln_file="$WorkingDir" . "$aln_file"; #CHANGE - ADI - Selecton
    unless (open ALN2, $path_to_aln_file){ #from SEQRES-PDB clustalw alignment
        $html_error = "sys";
        $log_error = "colorCoding.v2.pl : Read_aln_file: Cannot open the file $aln_file $!";
        &print_error_and_exit;
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

    unless (open PDB, $pdbdata){ #from SEQRES-PDB clustalw alignment 
        $html_error = "sys";
        $log_error = "colorCoding.v2.pl : Read_pdb_data: Cannot open the file $pdbdata $!";
        &print_error_and_exit;
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
# also read pvalue on selecton server.
sub Read_grades {
    # check if $prog_out exists and if it contains any data
    if (!-e $prog_out or -z $prog_out){
        $html_error = "An error has occured during the calculation. We apologize for the inconvenience. Please contact us for support.";
        $log_error = "colorCoding.v2.pl : Read_grades : The file $prog_out does not exist or contains no data";
        &print_error_and_exit;
    }
    # open the output file of kaks4site for reading
    unless (open GRADES, "<$prog_out"){
        $html_error = "sys";
        $log_error = "colorCoding.v2.pl : Read_grades : Cannot open the file $prog_out for reading $!";
        &print_error_and_exit;
    }
    my $CountPdb = 0; # for the 3LATOM
    my $CountSeq = 0; # for the SEQ
    while (<GRADES>) {
        my $line = $_;
        chomp $line;
        my $AA;
	my $grade;
	my $p_value; # used only for ML
	if ($Prog eq "Bayesian"){
	    $line =~ /^\s*(\d+)\s+(\S+)\s+(\S+)/;
	    $AA = $2;
	    $grade = $3;	
	}
	if ($Prog eq "ML"){
	    $line =~ /^\s*(\d+)\s+(\S+)\s+(\S+)\s+(\S+)/;
	    $AA = $2;
	    $grade = $3;
	    $p_value = $4;
	}
	if ($AA eq "*" || $AA eq "_" || $AA eq "-" ||$AA eq ""){ #gaps
	    next;
	}
	$Output[$CountSeq]{SEQ} = $AA;
	$Output[$CountSeq]{GRADE} = $grade;
	$Output[$CountSeq]{ATOM} = $PdbSeq[$CountSeq];
	if ($Prog eq "ML") {
	    $Output[$CountSeq]{PVALUE} = $p_value;  
	}
	if ($PdbSeq[$CountSeq] ne "-") {
	    $Output[$CountSeq]{PE} = $PDB[$CountPdb];
	    $CountPdb++;
	}
	else {
	    $Output[$CountSeq]{PE} = "-";
	}
	$CountSeq++;  
    }
    close GRADES;
}

################################################################################
#reading the file which contains a list of the residues and what category of color they fall in (selection4site.txt):
################################################################################
sub Read_res_to_color{

    if (!-e $file_colorRes or -z $file_colorRes){
        $html_error = "An error has occured during the calculation. Please try running the program again.";
        $log_error = "colorCoding.v2.pl : Read_res_to_color : The file $file_colorRes does not exist or contains no data";
        &print_error_and_exit;
    }
    # open the output file of kaks4site for reading
    unless (open RES_2_COLORS, "<$file_colorRes"){
        $html_error = "sys";
        $log_error = "colorCoding.v2.pl : Read_res_to_color : Cannot open the file $file_colorRes $!";
        &print_error_and_exit;
    }   
    my $CountSeq = 0; # for the SEQ
    while (<RES_2_COLORS>) {
        my $line = $_;
        chomp $line;
        if ($line =~ /^\s*(\d+)\s+(\S+)\s+(\d+)/){ 
            my $ID = $1;
            my $AA = $2;
            my $color = $3;		
            if ($AA eq "*" || $AA eq "_" || $AA eq "-"){ #gap			
                next;
            }
    #	    if ($Output[$CountSeq]{SEQ} ne $AA) {
    #            open LOG, ">>$OutLogFile";
    #            print LOG "\ncolorCoding.v2.pl:\Read_res_to_color: at position $CountSeq AA is : $Output[$CountSeq]{SEQ} different in color file and grade file: $AA\n";
    #            close LOG;
    #	    }
            $Output[$CountSeq]{COLOR} = $colorstep[$ColorScale_reverse{$color}]; 
            $Output[$CountSeq]{SEQ3DCOLOR} = $ColorScale_reverse{$color};
            $Output[$CountSeq]{NUMBER} = $ID;
            $CountSeq++;
        }
    }
    close RES_2_COLORS;
}

##########################################################################
# Print consurf.spt file (RasMol coloring script) for PE
# in case the residue has an additional letter after its name, we remove it, since rastop cannot
# color correctly in this case.
sub Print_consurf_spt {
    my $i;
    my $element;
    my $color;
    my $LineMaxElem = 10;  #No of elements on each consurf.spt line
    my $ElemCount = 0;
    my $separator;
    my $chain;

#============================ RasMol (CON) =============
    unless (open RSML, ">".$rasmol_file){
        $html_error = "sys";
        $log_error = "colorCoding.v2.pl : Print_consurf_rsml : Cannot open the file $rasmol_file for writing $!";
        &print_error_and_exit;
    }    
#========================================================
    # to solve the problem when chime runs out of colors.
    #print SPT "select all\n";
#============================ RasMol (CON) =============
    print RSML "select all\n";
#=======================================================

    #print SPT "color [200,200,200]\n\n";
#============================ RasMol (CON) =============
    print RSML "color [200,200,200]\n\n";
#=======================================================

    my $file_count = 9;
    foreach $color (@colorstep) {        
        my $color_elem_count = 0;
        #open COLOR, ">./pdbspt/" . $file_count . ".spt";

        foreach $element (@Output) {
            my $elem_to_print = "";
            if ($$element{COLOR} eq $color and $$element{PE} ne "-" ) {                
                if ($$element{PE} =~ /^(\w\w\w\d+)\D+(:.*)/){
                    $elem_to_print = $1.$2;
                }
                else{
                    $elem_to_print = $$element{PE};
                }
                if ($ElemCount == 0) {
                #print SPT "\nselect ";
                #print COLOR "\nselect ";
#============================ RasMol (CON) =============
                print RSML "\nselect ";
#=======================================================
                }
                if ($ElemCount == 10) {
		    #       $$element{PE} =~ s/(\:)([\d|\w])//;
		    #     $chain=$2;
                    #print SPT "\nselect selected or $$element{PE}";
                    #print COLOR "\nselect selected or $$element{PE}";
#============================ RasMol (CON) =============
                    print RSML "\nselect selected or $elem_to_print";
#======================================================
                    $ElemCount = 0;
                }
                else {
		    #                    $$element{PE} =~ s/(\:)([\d|\w])//;
		    #                   $chain=$2;
                    #print SPT "$separator $$element{PE}";
                    #print COLOR "$separator $$element{PE}";
#============================ RasMol (CON) =============
        		    print RSML "$separator $elem_to_print";
#======================================================
                    $separator = ",";
                }
                $ElemCount++;
                $color_elem_count++;
            }
        }
    	if ($ElemCount > 0) {
            #print SPT "\nselect selected and:$chain\n";
            #print SPT "color $color\nspacefill\n";    
    #============================ RasMol (CON) =============
            print RSML "\nselect selected and:$chain\n";
            print RSML "color $color\nspacefill\n";
            # 10/10/07. according to a check we did, seems that this line is not relevant for Selecton
            #print RSML "\ndefine CON" . $file_count . " selected\n\n";
    #=======================================================
            #print COLOR "\nspacefill\n";
        }
        if ($color_elem_count == 0){
            #print COLOR "javascript alert(\"No residues have this color\")";
        }
        $ElemCount = 0;
        $separator = " ";
        $file_count--;        
        #close COLOR;
    }    
    #close SPT;
    close RSML;
}

#################################################################################
# print consurf.js file for PE Seq3D
# NOT RELEVANT FOR FGIJ !!
#sub Print_consurf_js {
#    my $i;
#    my $element;
#    my $Count = 0;
#    my $total = scalar (@Output);
#    unless (open JS, ">consurf.js"){
#	open OUTPUT, ">>$OutHtmlFile";
#	print OUTPUT $SysErrorDef;
#	print OUTPUT $ContactDef;
#	close OUTPUT;
#	open LOG, ">>$OutLogFile";
#	print LOG "\ncolorCoding.v2.pl:\nPrint_consurf_js: Can\'t open the fileconsurf.js\n ";
#	close LOG;
#	print "Can\'t open the fileconsurf.js\n ";
#	open ERROR, ">error";
#	close ERROR;
#	exit;
#    }
#    print (JS "var csc = new Array(", $total +1 , ");\n\n");
#    print (JS "for (i = 0; i <= " ,$total  , "; i++)\n");
#    print JS "\tcsc[i] = -1;\n\n";
#    foreach $element (@Output) { 
#        if ($$element{PE} ne "-" ) {
#            print (JS "csc[" , $$element{NUMBER}, "] = ", $$element{SEQ3DCOLOR}, ";\n");       
#        }
#    }
#}

############################################################################
# Print the 'amino acid conservation scores' file
sub Print_PE {	
    # open the final output file for writing
    unless (open PE, ">$final_out"){
        $html_error = "sys";
        $log_error = "colorCoding.v2.pl : Print_PE : Cannot open the file $final_out for writing $!";
        &print_error_and_exit;
    }
    my $global_kaks;
    if ($Prog eq "ML"){
        $global_kaks = &Get_global_KaKs;
        print PE "\t Ka/Ks Scores\n";
        print PE "\t===============================\n\n";
        print PE "- POS: The position of the AA in the SEQRES derived sequence.\n";
        print PE "- SEQ: The SEQRES derived sequence in one letter code.\n";
        print PE "- KaKs: The non-synonymous to synonymous ratio.\n";
        print PE "- P-Value: The pvalue of the AA in the SEQRES derived sequence.\n";
        print PE "- 3LATOM: The ATOM derived sequence in three letter code, including the AA's positions as they appear in the PDB file and the chain identifier.\n";
        print PE "- COLOR: The color scale representing Ka/Ks (7 - conserved, 1 - positive selection).\n";
        print PE "- MSA DATA: The number of aligned sequences having an amino acid (non-gapped) from the overall number of sequences at each position.\n";
        print PE "- RESIDUE VARIETY: The residues variety at each position of the multiple sequence alignment.\n\n";
        print PE " The Ka/Ks score of the all protein is: $global_kaks. \n\n";
        print PE " POS\t SEQ\t    3LATOM\t  Ka/Ks\t      P-Value\tCOLOR\tMSA DATA\tRESIDUE VARIETY\n";
        
        my $pos = 1;
        foreach my $elem (@Output){
            printf (PE "%4d", "$pos");
            printf (PE "\t%4s", "$$elem{SEQ}");
            printf (PE "\t%10s", "$$elem{PE}");
            printf (PE "\t%8.3f", "$$elem{GRADE}");
            printf (PE "\t%2.3f", "$$elem{PVALUE}");
            printf (PE "\t%4d", "$ColorScale{$$elem{SEQ3DCOLOR}}");
            printf (PE "\t%8s", "$$elem{reliability}");
            printf (PE "\t%-s\n", "$$elem{res_range}"); # to align left
            $pos++;
        }
    } #Bayes
    else {
        print PE "\t Ka/Ks Scores\n";
        print PE "\t===============================\n\n";
        print PE "- POS: The position of the AA in the SEQRES derived sequence.\n";
        print PE "- SEQ: The SEQRES derived sequence in one letter code.\n";
        print PE "- KaKs: The non-synonymous to synonymous ratio.\n";
        print PE "- 3LATOM: The ATOM derived sequence in three letter code, including the AA's positions as they appear in the PDB file and the chain identifier.\n";
        print PE "- COLOR: The color scale representing Ka/Ks (7 - conserved, 1 - positive selection).\n";
        print PE "- MSA DATA: The number of aligned sequences having an amino acid (non-gapped) from the overall number of sequences at each position.\n";
        print PE "- RESIDUE VARIETY: The residues variety at each position of the multiple sequence alignment.\n\n";
        
        print PE " POS\t SEQ\t    3LATOM\t  Ka/Ks\t      \tCOLOR\tMSA DATA\tRESIDUE VARIETY\n";
        my $pos = 1;
        foreach my $elem (@Output){	
            printf (PE "%4d", "$pos");
            printf (PE "\t%4s", "$$elem{SEQ}");
            printf (PE "\t%10s", "$$elem{PE}");
            printf (PE "\t%8.3f", "$$elem{GRADE}");
            printf (PE "\t%4d", "$ColorScale{$$elem{SEQ3DCOLOR}}");
            printf (PE "\t%8s", "$$elem{reliability}");
            printf (PE "\t%-s\n", "$$elem{res_range}"); # to align left
            $pos++;
        }
    }
    close PE;
}

############################################################################
# Read global result 
sub Get_global_KaKs {	    
    unless (open GLOBAL, "<$global_res"){
        $html_error = "sys";
        $log_error = "colorCoding.v2.pl : Get_global_KaKs: Cannot open the file $global_res for reading $!";
        &print_error_and_exit;
    }
    while (<GLOBAL>){	
        my $line=$_;
        my @global   = split(/\s/,$line);
            if ($global[1] =~ m/ka/){
                return $global[2];
        }
    }
    # if we arrived here, it is because we didn't return from the while loop
    $html_error = "sys";
    $log_error = "colorCoding.v2.pl : Get_global_KaKs : Cannot find the global kaks in file $global_res";
    &print_error_and_exit;
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
    unless (open PE, "<$final_out"){
        $html_error = "sys";
        $log_error = "colorCoding.v2.pl : insert_grades : Cannot open the file $final_out for reading $!";
        &print_error_and_exit;
    }
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
    #&PDB::insert_tempFactor($pdb_file, $chain, \%grades);
}

#######################################################################
# calc the number of informative data per each position
sub calc_reliability{
    my @aln = ();
    &read_aln(\@aln);   
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
# read the MSA file into array of hashes
    unless (open MSA, "<$msa"){
        $html_error = "sys";
        $log_error = "colorCoding.v2.pl : calc_reliability : Cannot open the file $msa for reading $!";
        &print_error_and_exit;
    }        
    my $inFile  = Bio::SeqIO->new('-file' => "$msa" , '-format' => 'Fasta');
    my $numberOfSeqs=0;
    my $seqLength;
    
    #### reading alignment data into an array of a hash 
    while ( my $seqObj = $inFile->next_seq() ) {
        if ($numberOfSeqs==0){
            $seqLength=length($seqObj->seq());
        }
	my $seq = $seqObj->seq();
	my @seq=split("",$seq);
	my $name = $seqObj->display_id();
	if ($seqObj->desc() ne ""){
	    $name .= " ".$seqObj->desc();
	}   
	for (my $counter=0; $counter < $seqLength; $counter++) {
	    $Aln[$counter]{$name} = $seq[$counter];
	}
	$numberOfSeqs++;
    }
    #### calculating the reliability at each position: i,.e how much sequence data (non-gaps) exists at each position
    my $aaInRefSeq=0;
    for (my $j=0; $j < $seqLength; $j++) {
        if ($Aln[$j]{$SeqName} ne "-" and $Aln[$j]{$SeqName} ne "X"){	#not a gap or unknown in the reference sequence
            my $inf = 0;
            # loop over the different values
            foreach my $aa (values %{$Aln[$j]}){
                unless ($aa eq "-"){	    
                    $inf++;
                }
            } 
            # insert the reliability field to the general array
            my $result = $inf . "/" . $numberOfSeqs;
            $Output[$aaInRefSeq]{reliability}=$result;
            $aaInRefSeq++;
        }	
    }
    
    # calculating the amino acid content at each position
    my $CountPos = 0;
    foreach my $pos (@Aln){
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
    #my $PosCount = 0;
    #my $first_name = "";
    #my $i;
    #while (<MSA>){
    #	if ($_ =~ /^CLUSTAL/){
    #	    next;
    #	}
    #	elsif ($_ =~ /^(\S+)\s+(\S+)\s*/){
    #	    my $name = $1;
    #	    my $seq = $2;
    #	    my @seqs = split '', $seq;
    #	    # to set the first name
    #	    if ($first_name eq ""){
    #			$first_name = $name;
    #	    }
    #	    # a new block
    #	    elsif ($name eq $first_name){
    #			$PosCount = $i; 
    #	    }
    #	    $i = $PosCount;
    #	    foreach my $cell (@seqs){
    #			$$aln[$i]{$name} = $cell;
    #			$i++;
    #	    }
    #	}
    #}
    close MSA;
}

#####################################################################
# calc the residues variety per each position
sub calc_residue_range {
    my @aln = ();  
    &read_aln(\@aln); 
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
#####################################################################
sub print_error_and_exit{
    open ERROR, ">".$WorkingDir."error";
    print ERROR "HTML: $html_error\n";
    print ERROR "LOG: $log_error\n";
    close ERROR;
    exit;
}


