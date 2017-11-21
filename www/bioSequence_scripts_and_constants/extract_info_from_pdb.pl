#!/usr/bin/perl

#######################################################################
# The script extracts the SEQRES and ATOM sequences from the PDB file,
# and write them in files
# Script was edited by Elana Erez, July 2007, to match the needs of the new ConSurf
#######################################################################

#*** Modules ***#
use lib "/bioseq/bioSequence_scripts_and_constants"; #"/bioseq/scripts_for_servers";  
use lib "/bioseq/ConSurf"; 

use GENERAL_CONSTANTS;
use CONSURF_CONSTANTS;
use strict;


my $file = shift;  #file to read
my $chain = shift; #chain to read
my $pdb = shift;
my $mode = shift; # 'msa' or 'fasta'
my $run_name = shift;
my $OutputHtml = shift;
my $pdb_data = shift;
my $pdb_data_insertion = shift;
my $insertionsListFile = shift;
my $out_error_file = shift;   
my $atom_length_file = shift;  # "atom.length"
my $seqres_length_file = shift;  # "seqres.length"
my $aa_seq_file = shift; # was $NamePref. ".seq" , should be renamed to "protein_aa.seq"
my $seqres_pdb_strings_file = shift; # was $NamePref. ".fasta2"
my $pdb_in_fasta_file = shift; # was $NamePref. ".pdbfasta"
my $pdb_title_file = shift; # was $NamePref. ".title"
my $QuickHelp = shift;
my $server_name = shift;  # currently : selecton or consurf

my %AAdata = ( 
    "ALA" => {"name"    => "A", "area" => "222" }, 
    "ARG" => {"name"    => "R", "area" => "366" }, 
    "ASN" => {"name"    => "N", "area" => "274" }, 
    "ASP" => {"name"    => "D", "area" => "267" }, 
    "CYS" => {"name"    => "C", "area" => "249" }, 
    "GLN" => {"name"    => "Q", "area" => "300" }, 
    "GLU" => {"name"    => "E", "area" => "294" }, 
    "GLY" => {"name"    => "G", "area" => "193" }, 
    "HIS" => {"name"    => "H", "area" => "309" }, 
    "ILE" => {"name"    => "I", "area" => "301" }, 
    "LEU" => {"name"    => "L", "area" => "301" }, 
    "LYS" => {"name"    => "K", "area" => "336" }, 
    "MET" => {"name"    => "M", "area" => "314" }, 
    "MSE" => {"name"    => "M", "area" => "314" },#SELENIOMETHIONINE
    "PHE" => {"name"    => "F", "area" => "329" }, 
    "PRO" => {"name"    => "P", "area" => "257" }, 
    "SER" => {"name"    => "S", "area" => "235" }, 
    "THR" => {"name"    => "T", "area" => "260" }, 
    "TRP" => {"name"    => "W", "area" => "372" }, 
    "TYR" => {"name"    => "Y", "area" => "349" }, 
    "VAL" => {"name"    => "V", "area" => "274" },                 
    ); 
my %Mod_AAdata = (
    "MSE" => {"name"	=> "MET"},
    "ASN" => {"name"	=> "ASN"},
    "MLY" => {"name"	=> "LYS"},
    "HYP" => {"name"	=> "PRO"},
    "SER" => {"name"	=> "SER"},
    "THR" => {"name"	=> "THR"},
    "CME" => {"name"	=> "CYS"},    	
    "CGU" => {"name"	=> "GLU"},
    "SEP" => {"name"	=> "SER"},
    "KCX" => {"name"	=> "LYS"},
    "MLE" => {"name"	=> "LEU"},
    "TPO" => {"name"	=> "THR"},
    "CSO" => {"name"	=> "CYS"},
    "PTR" => {"name"	=> "TYR"},
    "DLE" => {"name"	=> "LEU"},
    "LLP" => {"name"	=> "LYS"},    
    "DVA" => {"name"	=> "VAL"},
    "TYS" => {"name"	=> "TYR"},
    "AIB" => {"name"	=> "ALA"},
    "DVA" => {"name"	=> "VAL"},
    "TYS" => {"name"	=> "TYR"},
    "AIB" => {"name"	=> "ALA"},
    "OCS" => {"name"	=> "CYS"},
    "NLE" => {"name"	=> "LEU"},
    "MVA" => {"name"	=> "VAL"},
    );
my $modified_residue_counter = 0;
my @sequence = (); # to hold the SEQRES sequence
my @residues = (); # to hold the ATOM sequence
my $seq_string = "";  # the SEQRES string
my $atom_string = ""; # the ATOM string
my ($seq_length, $atom_length, $title, $ans, $log_error, $html_error);  
########## insertions ##########
my  $insertionsList = ',';
###################################

#---- only for Selecton: we check if it is NMR and output a message accordingly
# --- added also for consurf, as the old server doesn' check it anywhere else
if ($server_name eq "selecton" or $server_name eq "consurf"){
    open PDBFILE, "<$file";
    my @orig_pdb = <PDBFILE>;
    close PDBFILE;
    # check if it is NMR and there is more than one model
    my $many_models = 0;
    foreach my $i (@orig_pdb){
        if ($i =~ /^MODEL\s+2/){
            $many_models = 1;
            last;
        }
    }
    # cut all the other models and print a warning
    if ($many_models == 1){
        open PDBFILE, ">".$file;
        foreach my $line (@orig_pdb){
            print PDBFILE $line;
            last if ($line =~ /^ENDMDL/);
        }
        close PDBFILE;        
        open OUTPUT, ">>$OutputHtml";
        print OUTPUT "\n<p><ul><li><font color='red'><b>Warning:</font></b> The PDB file contains more than one NMR structure. The first model is arbitrarily chosen in order to display the conservation grades using FirstGlance in Jmol. (If you wish to display the results on a different model, create your own PDB file and cut all the ATOM records of the other models).</li></ul></p>\n";
        close OUTPUT;
    }    
} # NMR for selecton

my $printwarning = 0; # in order to print the warning only once
### read the SEQRES and ATOM sequences (use the PDB module)
# in case of MSA, there is no need to read the SEQRES
if ($server_name eq "consurf" and $mode eq "fasta"){#--- ONLY CONSURF---#
    unless(open PDB, $file){
        &error_and_exit("sys", "extract_info_from_pdb.pl : cannot open file $file for reading $!");
    }
    @sequence = &read_SEQRES_pdb ("consurf");
}
unless(open PDB, $file){
    &error_and_exit("sys", "extract_info_from_pdb.pl : cannot open file $file for reading $!");
}
@residues = &read_complex;
$chain = " " if ($chain eq "NONE"); 
open (PDBDATA, "> " .$pdb_data);   # to PE 
open (PDBDATAINSERTION, "> " .$pdb_data_insertion) if ($server_name eq "consurf");   # to PE         #--- ONLY CONSURF---#
# check the ATOM sequence
my @ans = &check_sequence("ATOM", \@residues, $OutputHtml, $QuickHelp);
close PDBDATAINSERTION if ($server_name eq "consurf"); #--- ONLY CONSURF---#
close PDBDATA;

if ($ans[0] eq "ok"){
    $atom_string = $ans[1];
    $atom_length = $ans[2];
}
else{
    &error_and_exit($ans[1], $ans[2]);
}
#--- ONLY CONSURF---#
if ($server_name eq "consurf"){
    open (INSERTIONSLISTFILE, ">$insertionsListFile");
    print INSERTIONSLISTFILE "$insertionsList";
    close INSERTIONSLISTFILE;
    #############
    #for PE New Version
    open (ATOMLENGTH, ">".$atom_length_file);   # to PE 
    print ATOMLENGTH $atom_length;
    close ATOMLENGTH;
    
    my @sequencePE;
    my $seq_string_PE;
    my $seq_length_PE;
    
    unless(open PDB, $file){
        &error_and_exit("sys", "extract_info_from_pdb.pl : cannot open file $file for reading $!");
    }
    @sequencePE = &read_SEQRES_pdb ("consurf_PE");
    @ans = &check_sequence("SEQRES", \@sequencePE, $OutputHtml);
    if($ans[0] eq "ok"){
        $seq_string_PE = $ans[1];
        $seq_length_PE = $ans[2];
    }
    else{
        &error_and_exit($ans[1], $ans[2]);
    }
    
    open (SEQRESLENGTH, ">".$seqres_length_file);
    print SEQRESLENGTH $seq_length_PE;
    close SEQRESLENGTH;    
    
    ##############
    # check the SEQRES sequence
    if ($mode eq "fasta"){  
    @ans = &check_sequence("SEQRES", \@sequence, $OutputHtml);
    if($ans[0] eq "ok"){
        $seq_string = $ans[1];
        $seq_length = $ans[2];
    }
    else{
        &error_and_exit($ans[1], $ans[2]);  
    }
        # in case there is no SEQRES
        if ($seq_string eq "NO_SEQRES"){
            $seq_string = $atom_string;
            $seq_length = $atom_length;
        }
    }
    ### print the sequences in the related files    
    if ($mode eq "fasta"){    
        # the SEQRES file
        open SEQ, ">" .$aa_seq_file;
        print SEQ ">" .$pdb . "_" . $chain . "  " . $seq_length . "\n" . $seq_string . "\n" ;
        close SEQ;
    
        # the SEQRES and ATOM file
        open OUT, ">" .$seqres_pdb_strings_file;
        print OUT ">SEQRES_$chain\n";
        print OUT "$seq_string\n\n";
        print OUT ">PDB_$chain\n";
        print OUT "$atom_string\n";
        close OUT;
    }
}
#--- ONLY CONSURF finish---#

# the ATOM file
open PDBFASTA, ">" .$pdb_in_fasta_file;
print PDBFASTA ">PDB_$chain\n";
print PDBFASTA "$atom_string\n";
close PDBFASTA;

# the title file
open TITLE, ">" .$pdb_title_file;   
print TITLE $title;
print TITLE "\nCHAIN: $chain\n";
close TITLE;


################################################################
# check if the SEQRES or ATOM sequence exists and contain standard amino-acids, 
# and write the file prefix.pdbdata
# returns $atom_string, $atom_length
sub check_sequence {

    my $field = shift; # 'ATOM' or 'SEQRES'
    my $array = shift;
    my $OutputHtml  = shift;

    my $seq = "";
    my $anyAA = 0;  # to check if there is a sequence at all
    my $CountAA = 0; # to count the AAs of the given chain
    my $StandardAA = 0;
    my %uniqueInsertion;
    my $html_error;
 
    foreach my $AA (@{$array}) {
	$anyAA++ unless ($$AA{name3L} =~ /^\s*$/);
        #for chain $chain only
        if ($$AA{chain} eq $chain) { 	    
            $CountAA++;	    
	    if ($$AA{name1L} ne ""){ # a standard AA
		$seq .= $$AA{name1L}; #AA to print
		$StandardAA++;
		if ($field eq "ATOM"){
###############insertion begin #####################
		    print PDBDATA "$$AA{name3L}$$AA{number}:$$AA{chain}\n";#PDBDATA was opened by the function who called the check_sequence procedure
		    #--- ONLY CONSURF---#
                    if ($server_name eq "consurf"){
                        if ($$AA{number} =~ /[A-Z]$/){
                            my $insertionNumber = $$AA{number};
                            $insertionNumber =~ s/[A-Z]$//;
                            $insertionsList .= "$insertionNumber\," if ($insertionsList !~ /,$insertionNumber\,/);		
                        }	
                        print PDBDATAINSERTION "$$AA{name3L}$$AA{number}:$$AA{chain},$$AA{name3L}$$AA{number},$$AA{atomNo}\n";
                    }  #--- ONLY CONSURF finish---#
############### insertion end ######################		   
		}
	    }
        }
    }
    ### There is no sequence at all
    if ($anyAA == 0){
	# no SEQRES
	if ($field eq "SEQRES"){
	    $seq = "NO_SEQRES";            
            if ($printwarning == 0){
                open OUTPUT, ">>$OutputHtml" or die("cannot open $OutputHtml $!");
                print OUTPUT "\n<p><ul><li>Warning: There is no SEQRES derived information in the PDB file. The calculation will be based on the ATOM derived sequence. If this sequence is incomplete, we recommend to run $server_name again using an external multiple sequence alignment file, which is based on the complete protein sequence.</li></ul></p>\n";
                close OUTPUT;
                $printwarning =1;
            }
	    return "ok", $seq, $StandardAA;
	}
	# no ATOM
	else {
            $html_error = "<b>There is no $field derived information in the PDB file.<br>Please refer to the OVERVIEW for detailed information about the PDB format.</b></p>\n";
            $log_error = "pdb_to_fasta_insertion.pl: exit - There is no $field field.\n";

            return "error", $html_error, $log_error;
	}
    }

    my $chain_link = $QuickHelp . "#note2";
    ### The given chain is not correct
    if ($CountAA == 0){

	# the given chain was "none" 
	if ($chain eq " "){
	    if ($pdb eq "FILE"){  # an uploaded PDB file
		$html_error ="<b>The <a href=$chain_link>chain identifier</a> column in the $field field of the <a href=$file target=PDB>PDB file</a> is not empty.<br>Please check your file or run $server_name again with a specific chain identifier.</b></p>\n";
		
	    }
	    else {
		$html_error ="<b>The <a href=$chain_link>chain identifier</a> column in the <a href=$file target=PDB>PDB file</a> is not empty.<br>Please run $server_name again with a specific chain identifier (A, B, C, etc.)</b></p>\n";
	    }
	}	
	else {
	    if ($pdb eq "FILE"){  # an uploaded PDB file

		$html_error = "<b>Chain \"$chain\" does not exist in the $field field of the <a href=$file target=PDB>PDB file</a>.<br> Please check your file or run $server_name again with another chain. If there is no <a href=$chain_link>chain identifier</a>, set the 'Chain Identifier' field to \"none\".</b></p>\n";
	    }
	    else {
		$html_error = "<b>Chain \"$chain\" does not exist in the <a href=$file target=PDB>PDB file</a>.<br> Please run $server_name again with another chain. If there is no <a href=$chain_link>chain identifier</a>, set the 'Chain Identifier' field to \"none\".</b></p>\n";
	    }
	}
	$log_error = "pdb_to_fasta_insertion.pl: exit - the chain is not correct\n";
        return "error", $html_error, $log_error;
    }
    ### Too much nodified residues
    
    if ($modified_residue_counter>0 and $CountAA>0)
	{
		if ($modified_residue_counter/$CountAA>CONSURF_CONSTANTS::MAXIMUM_MODIFIED_PERCENT)
		{
	        	&error_and_exit("Too much modified or not standard amino acid recognized by the $server_name server<br />", "read_SEQRES_pdb : found $modified_residue_counter irregular char 'in SEQRES line\n");
		}
	}
    ### The chain exist but it is not a protein
    if ($StandardAA == 0){
	if ($pdb eq "FILE"){  # an uploaded PDB file
	    $html_error = "<b>Chain \"$chain\" in the $field field of the <a href=$file target=PDB>PDB file</a> does not contain any standard amino acids.<br />Please check your file or run $server_name again with another chain.</b></p>\n";
	}
	else {
	    $html_error = "<b>Chain \"$chain\" in the $field field of the <a href=$file target=PDB>PDB file</a> does not contain any standard amino acids.<br />Please run $server_name again with another chain.</b></p>\n";
	}	
	$log_error = "pdb_to_fasta_insertion.pl: exit - the chain is not a protein.\n";
        return "error", $html_error, $log_error;
    }
    return "ok", $seq, $StandardAA;
}
#********************************************* #
# sub read_sequence_pdb($) 
# ------------------
# reads a pdb sequence file containg
# into @residues array, on a 1LETTER format
# This array will contain an entry per each 
# amino acid, containing all the relevant data
#********************************************* #
sub read_SEQRES_pdb {
    
    my $calling_server=shift; # Who Calls the function - Considering modified residues only in ConSurf/ConSurf_PE)
    my $_chain; #the chain
    my $ResList; #list of residue names.

    my @names3L = (); #here are written the names in 3L   
    my $resCount = 0; 
    my $AA3L;
    my @residues = ();
    my @sequence = (); #main data to be returned
    
    if (!defined $calling_server)
    {
    	$calling_server="";
    }
    while (<PDB>) {        
        if (/^SEQRES/){            
            $_chain = substr ($_, 11, 1);
            $ResList = substr ($_, 19, 51);    
            @names3L = split (/\s+/, $ResList);
      
            #*************************
            #make data structure like in read_complex
 
            foreach $AA3L (@names3L) {
	    	if($AAdata{$AA3L}{name} eq "" && $AA3L ne ""){ #MODIFIED RESIDUE
                    # count this modified residue
                    $modified_residue_counter++;

		    if ($Mod_AAdata{$AA3L}{name} ne "" and ($calling_server eq "consurf" or $calling_server eq "consurf_PE")) # Known Modified Residue
		    {
		    	
			$sequence[$resCount]{name1L} = $AAdata{$Mod_AAdata{$AA3L}{name}}{name};   #residue name 1L
                    	$sequence[$resCount]{name3L} = $Mod_AAdata{$AA3L}{name}; #residue 3L name
                    	$sequence[$resCount]{chain} = $_chain;  #chain        
                    	$resCount++;
		    	if ($chain eq $_chain)
			{
				open OUTPUT, ">>$OutputHtml";
				print OUTPUT "<p><ul><li><font color='red'><b>Warning:</font></b> Please note: The modified residue \'$AA3L\' found in your pdb file in the line:<br />$_<br />, was converted back to the original before the analysis took place</li></ul></p>";
				close OUTPUT;
			}
		    }
		    elsif ($Mod_AAdata{$AA3L}{name} eq "" and ($calling_server eq "consurf" or $calling_server eq "consurf_PE"))
		    {
		    	$sequence[$resCount]{name1L} = "X";   #residue name 1L
			$sequence[$resCount]{name3L} = $AA3L; #residue 3L name
			$sequence[$resCount]{chain} = $_chain;  #chain        
			$resCount++;
			if ($chain eq $_chain and ($calling_server eq "consurf"))
			{
				chomp($_);
                    		open OUTPUT, ">>$OutputHtml";
				print OUTPUT "<p><ul><li><font color='red'><b>Warning:</font></b> The residue \'$AA3L\' found in your pdb file in the line:<br />$_<br />This residue is not a standard 3 letters code for an amino acid recognized by the $server_name server and thus replaced by 'X'</li></ul></p>";
				close OUTPUT;
			}
		   }
		   else
		   {
		   	chomp($_);
			&error_and_exit("The residue \'$AA3L\' found in your pdb file in the line:<br />$_<br />This residue is not a standard 3 letters code for an amino acid recognized by the $server_name server<br />", "read_SEQRES_pdb : found irregular char \'$AA3L\'in SEQRES line\n") if ($chain eq $_chain);
		   }
                }
                else{
                    $sequence[$resCount]{name1L} = $AAdata{$AA3L}{name};   #residue name 1L
                    $sequence[$resCount]{name3L} = $AA3L; #residue 3L name
                    $sequence[$resCount]{chain} = $_chain;  #chain        
                    $resCount++;
                }
            }
        }
    }
    close (PDB);    
    return (@sequence); #returns entire array with data
}

##############################################
##############################################
#********************************************* #
# sub read_complex($) 
# ------------------
# reads a pdb file containg
# a pdb complex (more than one protein chain
# into @residues array.
# This array will contain an entry per each 
# amino acid, containing all the relevant data
#********************************************* #

sub read_complex {
    
    my $resCount = 0; 
    my @residues = ();
    my $i;
    my $Res3L; #residue before 3letter validation (AARG/ BARG)
    my $_chain;
    my $number;
    my $iCode; # insertion code
    my $insCode='';
    my $atomNo;
    my $chainInsCode = '';
    $residues[0]{number} = -300; #to initialize this scalar
    
    while (<PDB>) {        
        if (/^ATOM/){
  ############# insertion begin ###########
	    $atomNo .= substr ($_, 7, 5);
	############# insertion end ###########
            $Res3L = substr ($_, 17, 3);
            $_chain = substr ($_, 21, 1);
            $number = substr ($_, 22, 4);
	    $number =~  s/\s+//g;
	    $iCode = substr ($_, 26, 2);
	    if ($iCode !~ /[A-Z]/i){
		$chainInsCode = '';
	    }
########################$$$$$$$$$$$$$$$$$$###########################
	    $number =~ s/\s+//;
	    if (($iCode =~ /[A-Z]/i) && ($chainInsCode !~ /$iCode/))
	    {
		$number =~  s/\s+//g;		
	        $insCode .= $Res3L . "," . $_chain . "," . $number . "," . $iCode  . ";";	
		$chainInsCode =  $insCode;
		$chainInsCode =~ s/.*,(.)/$1/;
	    }
            ($number) =~ s/\s+//g;
  ############# insertion begin ###########
	    ($atomNo) =~ s/\s+//g;
             $atomNo .=','; 
	############# insertion end ###########
	    $iCode =~ s/\s*//g;
            my $numberAndInsertion = $number . $iCode; 
	    if ($numberAndInsertion ne  ($residues[$resCount - 1]{number})) { #create residue hash
        ############# insertion begin ###########
                $residues[$resCount-1]{atomNo} = $atomNo;
                $residues[$resCount-1]{atomNo} =~ /^(\d+),.*\,(\d+)\,$/;
                my $leftAtomNo = $1-1; #?
                my $rightAtomNo = $2-1; #?
                $residues[$resCount-1]{atomNo} = '(atomno >= ' . $leftAtomNo . ' and atomno <= ' . $rightAtomNo . ')';
	############# insertion end ###########
           
        ##############insertion begin #######################
		if ($iCode =~ /[A-Z]/i){    
                    $insCode =~ s/\;$//;
		    $residues[$resCount]{name1L} = $AAdata{$Res3L}{name};   #residue name 1L           
		    $residues[$resCount]{name3L} = $Res3L;   #residue name 3L
		    $residues[$resCount]{chain} = $_chain;  #residue chain
		    $residues[$resCount]{number} = $number . $iCode;
                    $resCount++;
		    $insCode='';
		    $atomNo ='';
	        }
                else
		{
		    $residues[$resCount]{name1L} = $AAdata{$Res3L}{name};   #residue name 1L           
		    $residues[$resCount]{name3L} = $Res3L;   #residue name 3L
                    $residues[$resCount]{chain} = $_chain;  #residue chain
		    $residues[$resCount]{number} = $number; #residue number
		    $resCount++;
		    $atomNo ='';
		}
                ############### insertion end #####################
            }
        }
        # READ TITLES
        elsif (/^TITLE(.+)/){
            $title .= $1;            
        }
        # for NMR files - read only one model
        elsif (/^ENDMDL/){
            last;
        }
    }
    close (PDB);    
    return (@residues); #returns entire array with data
}

##########3
sub error_and_exit{
    my $h_err = shift;
    my $l_err = shift;
    if ($server_name eq "consurf"){
        $out_error_file = GENERAL_CONSTANTS::SERVERS_RESULTS_DIR."ConSurf/$run_name"."/$out_error_file";
    }
    elsif ($server_name eq "selecton"){
        $out_error_file = GENERAL_CONSTANTS::SERVERS_RESULTS_DIR."Selecton/$run_name"."/$out_error_file";
    }
    open ERROR_OUT, ">".$out_error_file;
    print ERROR_OUT "HTML: $h_err\n";
    print ERROR_OUT "LOG: $l_err";
    close ERROR_OUT;
    chmod 0755, $out_error_file;
    exit;
}

