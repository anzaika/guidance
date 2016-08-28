#!/usr/local/bin/perl

use strict;
use Storable;
use lib "/bioseq/bioSequence_scripts_and_constants/"; #"/db1/System/bioseq/scripts_for_servers";  
use GENERAL_CONSTANTS;
use BIOSEQUENCE_FUNCTIONS;
use SELECTON_CONSTANTS;
use lib "/bioseq/Selecton/external_scripts";
use pipe_for_selecton;

my $WorkingDir = shift;
my $runCalcInput = $WorkingDir.shift;
my $formInput = $WorkingDir.shift;

# hashes to retrieve the runs' input from the storable element
my %FORM = ();
my %run_data = ();

# VARS FROM INPUT FILE
# vars read from form:
my ($epsilonPrecision, $query_seq_name_to_run, $upload_TREE_file, $optimizeBL);

# General vars
my ($OutLogFile, $querySeqFoundinMSA, $tree_faq, $PdbPrefix, $method, $run_name, $PdbFileNameUnc, $was_Pdb_uploaded, $proc_comm, $estimated_run_time, $begin_Q_runtime);

# HTML related vars
my ($SysErrorDef, $ContactDef, $WWWdir, $ErrorDef, $OutputURL);

# files which were created by the cgi
my ($fileDna_aligned, $treeUpload, $OutHtmlFile, $cgi_log_file, $sequences_names_file, $fileName_amino_aligned, $upload_unaligned_file_dna, $upload_MSA_file_dna, $clustal_aligned_file, $pdb_data);

my @sequences_names = (); # This array will hold the names of the sequences, as it appears in the input file.

my $qsub_ans_file = $WorkingDir."qsub_ans.txt"; # The flag file which is read by the daemon process
my $finish_flag = $WorkingDir."END_OK"; # The file which denotes that the run was finished.

&readInput();

my $http_path = GENERAL_CONSTANTS::SELECTON_URL;
### Sending e-mail from the cluster
my $smtp_server = GENERAL_CONSTANTS::SMTP_SERVER;
my $userName = GENERAL_CONSTANTS::ADMIN_USER_NAME;
my $userPass = GENERAL_CONSTANTS::ADMIN_PASSWORD;
my $mail = "mailto:".GENERAL_CONSTANTS::ADMIN_EMAIL."?subject=Selecton%20Run%20No:%20$run_name";
my $email_subject;
my $email_message;
my $email_system_return;
my $send_email_dir = GENERAL_CONSTANTS::SEND_EMAIL_DIR;

# VARS DEFINED IN THE SCRIPT:
# scripts or executables
my $biocluster_external_scripts_path = "/bioseq/Selecton/external_scripts/";
#my $selecton = "/bioseq/pupkoSVN/trunk/programs/selecton/selecton"; #"/d/bioinfo/users/adist/pupkoSVN/trunk/programs/selecton/selecton";
my $selecton = "/bioseq/Selecton/selecton.exe";
#my $selecton = $biocluster_external_scripts_path."srcSelecton/srcV2.2/selecton";
#my $mecSelecton = "/bioseq/pupkoSVN/trunk/programs/mec/mec";
my $mecSelecton = "/bioseq/Selecton/mec.exe";
#"/d/bioinfo/users/adist/pupkoSVN/trunk/programs/mec/mec"; # the exe of the kaks using the mec model is added by adid (29.1.07)
my $colorCoding = $biocluster_external_scripts_path . "colorCoding.v2.pl";
my $colorCodingLinear = $biocluster_external_scripts_path . "colorCodingLinear.pl";
my $statTest = "/cgi-bin/statTest.cgi"; # REMARK : SHOULD BE RE-WRITTEN AND RUN FROM THE CLUSTER
# vars
my $pdbUpload = $WorkingDir . $PdbFileNameUnc; #name for an uploaded PDB FILE
my $significance_test_faq = "/overview.html#meth5";
# files which will be created by this script
my $rsml = "rasmol.txt";#"colors.txt";
my $areTherePositiveSites = $WorkingDir."areSitesPositive.txt";
my $colorsLinear = "colors.html"; #results color-coded onto the linear sequence
my $outputScoreFile="kaks.res";
my $selection4Site_file = $WorkingDir."selection4Site.txt";
my $kaks_file = $WorkingDir."kaks.res";
my $params = "globalResult.txt";
my $log =  "kaks4site.log";
my $final_out = $PdbPrefix . ".gradesPE";
my $finalDNAFile_seq_names = "DNA.names.msa";
my $finalAminoFile_seq_names = "AMINO.names.msa";
my $tree_out = "kaks4site.tree";
my $rasmol_file = $PdbPrefix .".rsml";
my $statistics_file = SELECTON_CONSTANTS::STATISTICS_FILE;

#***************************** 
### FGiJ path and related files
my $FGiJ_path = "/fgij/";
my $FGiJ_pipe_pdb = $FORM{pdb_ID} . "_selecton" . $run_name . "_pipe.pdb";
my $FGiJ_link = $FGiJ_path . "fg.htm?mol=/results/" . $run_name . "/" . $FGiJ_pipe_pdb;
my $pipe_error = "pipe.error";

#---------------------------------------------------------------------------
#                          C A L C U L A T I O N
#---------------------------------------------------------------------------
$begin_Q_runtime = &printTime();
open LOG, ">".$OutLogFile;
print LOG &printTime();
print LOG "\nEntered selecton_run_calc.pl\nUpdating the HTML with running status\n";
# updating the HTML status to "running"
my $ans = &GENERAL_CONSTANTS::print_Q_status_in_html($OutHtmlFile, "Running", "no", $estimated_run_time);
print LOG $ans if ($ans ne "OK");


&run_calc();
        
### changing back the names in the DNA and AMINO files, so it will hold the original sequnces names.
&return_seq_names_to_files($fileName_amino_aligned, $finalAminoFile_seq_names, \@sequences_names);
&return_seq_names_to_files($fileDna_aligned, $finalDNAFile_seq_names, \@sequences_names);
#if a tree was created, we output it to the users, than first have to write to it original names
#if ((-e $WorkingDir.$tree_out) && !(-z $WorkingDir.$tree_out)){
#      print LOG "renaming numbered tree output $tree_out to : no_tree.txt\n";
#      my $cmd = "mv ".$WorkingDir.$tree_out." ".$WorkingDir."no_tree.txt";
#      my $out = `$cmd`;
#      print LOG "moving returned: $out\n";
#      chmod 0744, $WorkingDir."no_tree.txt";
#      print LOG "going to run change_tree_file_names\n";
#      #&change_tree_file_names(\@sequences_names, $WorkingDir."no_tree.txt", $WorkingDir.$tree_out);
#      print LOG "after change_tree_file_names\n";
#}
#------------
### create a pdb file with pipe in its header, for FGiJ to read
if ($was_Pdb_uploaded eq "yes"){
    print LOG "Going to prepare pipe file\n";
      &prepare_pipe(\@sequences_names);
}
else{
    print LOG "Not creating pipe, since the value of \$was_Pdb_uploaded is: \"$was_Pdb_uploaded\"\n";
}
###  PRINTING output FINAL NOTES AND LINKS

open OUTPUT, ">>" .$OutHtmlFile;
flock OUTPUT,2;
print OUTPUT "\n<H1><center><a name=finish>Selecton calculation is <font color=\"red\">FINISHED </font></a></center></H1>\n";	
print OUTPUT "<h3><i>Final Result:</i></h3>";
if ($was_Pdb_uploaded eq "yes") { #3D
    print OUTPUT "\n<p><b><A HREF='".$FGiJ_link."' TARGET=_blank>Graphical display of Selecton results</b></A> with FirstGlance in Jmol<br></p>\n";
}
else {
    print OUTPUT "\n<p><b><A HREF= $colorsLinear> View Color Coded Selecton Results </A></b></p>\n";
}
# check if there are P.S. sites. If yes - a button for statistical testing
open POS, "<$areTherePositiveSites";
my @pos = <POS>;
chomp @pos;
close POS;
if ($pos[0] eq "yes"){
    print OUTPUT "<font face=Verdana size=2>Positively selected sites found.</font><br>\n";
###### here is a link for statistical testing for the M8 and MEC models
    if (($FORM{MODEL} eq "M8") || ($FORM{MODEL} eq "MEC")) {
        my $ans = &add_data_to_input_file;
        # in case the data was not added - we don't create a submit button, we let the user know he can contact us
        if ($ans eq "OK"){            
            print OUTPUT "\n<FORM ENCTYPE=\"multipart/form-data\" ACTION=\"".$statTest."\" METHOD=\"POST\">\n";
            print OUTPUT "<INPUT TYPE=\"submit\" VALUE=\"Test statistical significance of positive selection\"><br>\n";
            print OUTPUT "<INPUT TYPE=hidden NAME=\"run_number\" VALUE=\"".$run_name ."\">\n";
            print OUTPUT "<INPUT TYPE=hidden NAME=\"input_file\" VALUE=\"".$runCalcInput."\">\n";
            print OUTPUT "<INPUT TYPE=hidden NAME=\"Q_log\" VALUE=\"".$OutLogFile."\">\n";
            print OUTPUT "<INPUT TYPE=hidden NAME=\"receipent\" VALUE=\"".$FORM{email_address}."\">\n";
            print OUTPUT "<font face=Verdana size=1><a href = \"".$significance_test_faq."\">This will run your data with a null model of evolution</a></font></FORM><br>\n";
        }
        else{
            print OUTPUT "<font face=Verdana size=2 color=\"green\">Please <a href=\"$mail\">contact us</a> if you wish to run a statistical test for your results, and mention this number: $run_name</font><br>";
        }
    }
}
else {
    print OUTPUT "<font face=Verdana size=2>No positively selected sites found in the protein.</font><br>\n";
}
print OUTPUT "<h3>Output Files:</h3>\n\n";
### print the links
if ($was_Pdb_uploaded eq "yes") { #if user ran selecton with a PDB struct.: Ka/Ks scores 2gether with color coding
    print OUTPUT "<p><A HREF= $colorsLinear TARGET=colors_window> Codon Ka/Ks scores color-coded on the linear sequence</A></p>\n";
}
print OUTPUT "<p><b><A HREF= $outputScoreFile> Codon Ka/Ks scores (numerical values)";

# add the gapped output only in case there was gapped output (if the files are identical, the only difference will be in the line 'Displayed on sequence 1< including gaps>'
if ((-s $WorkingDir."kaks.res.gaps") - (-s $WorkingDir.$outputScoreFile) > 15){
    print OUTPUT " - Reference sequence only</A></b></p>\n";
    print OUTPUT "<p><img SRC=\"/New.gif\" BORDER=0 height=30 width=40><b><A HREF= \"kaks.res.gaps\"> Codon Ka/Ks scores (numerical values) - For each position in the MSA, including gaps</A></b><img SRC=\"/New.gif\" BORDER=0 height=30 width=40></p>\n";
}
else{
    print OUTPUT "</A></b></p>\n";
}
if ( $upload_unaligned_file_dna ne "no"){ # print amino aln only if the user supplied a non-aligned file (then we run codon-align & produce an amino aln)
    print OUTPUT "<p><A HREF= $finalAminoFile_seq_names TARGET=MSA_window> Amino Acid Multiple Sequence Alignment (in Fasta format)</A>\n";
}
print OUTPUT "<p><A HREF= $finalDNAFile_seq_names TARGET=codonMSA_window> Codon Multiple Sequence Alignment (in Fasta format)</A>\n";
if ($upload_TREE_file eq "NOT_GIVEN") { # no user tree supplied - NJ tree created
    print OUTPUT "<p><A HREF= $tree_out TARGET=tree_window> Phylogenetic Tree</A>\n";
}
print OUTPUT "<p><A HREF = \"$params\">Likelihood and parameters of Selecton run</A>\n";
print OUTPUT "<p><A HREF = \"$log\">Log-file of Selecton run</A>\n";
if ($was_Pdb_uploaded eq "yes") { #if user ran selecton with a PDB struct.
    print OUTPUT "<p><A HREF= $rsml TARGET=spt_window>RasMol coloring script source</A></p>\n";
    print OUTPUT "<p><A HREF= $FGiJ_pipe_pdb TARGET=pipe_pdb>PDB file updated with Selecton results in its header</A></p>\n";
}
print OUTPUT "\n<br><br><p><center>Please <a href= \"$mail\">report any problem</a> in case of need.</center></p>\n";
flock OUTPUT,8;
close OUTPUT;

### write to the output that the job has finished
open OUTPUT, "<$OutHtmlFile";
flock OUTPUT,2;
my @output = <OUTPUT>;
flock OUTPUT,8;
close OUTPUT;   
open OUTPUT, ">$OutHtmlFile";
flock OUTPUT,2;
foreach my $line (@output){
	  if ($line =~ /Selecton Job Status Page/i){ #finds the phrase "Selecton" job status page, case-insensitive
            print OUTPUT "<H1 align=center>Selecton Job Status Page - <font color='red'>FINISHED</font></h1>\n";
		    print OUTPUT "<a href=#finish><H2 align=center>Go to the results</font></H2></a>\n";
	  }
	  else {
		   print OUTPUT $line; 
	  }
}
flock OUTPUT,8;
close OUTPUT;   
### stop the automatic reload
system 'echo "(cd '.$WorkingDir.' ; chmod -R og=rx * )" | /bin/tcsh';
&stop_reload;
# reporting the statistics file on a succssful ending
my $total_runtime = BIOSEQUENCE_FUNCTIONS::subtract_time_from_now($begin_Q_runtime);
open STATISTICS, ">>".$statistics_file;
flock STATISTICS, 2;
print STATISTICS "$run_name total runTime: $total_runtime\n";
flock STATISTICS, 8;
close STATISTICS;

$email_subject = "Your Selecton results for run number $run_name are ready";
$email_message = "Selecton finished calculation. Please click on the following link to view the results:\n$WWWdir"."output.html\nPlease note: the results will be kept on the server for three months.";
open LOG, ">>$OutLogFile";
print LOG "\nSending mail to user.\n";
GENERAL_CONSTANTS::send_mail("Selecton", $FORM{email_address}, $run_name, $email_subject, $email_message);
 
if (-e $WorkingDir."core"){
    print LOG "remove core file from working directory\n";
    unlink $WorkingDir."core";
}
	
print LOG "\nSelecton run completed successfully!";
print LOG "\n************** END OF LOG FILE *****************\n";
close LOG;

exit;  

#----------------------------------------------------------------------------------------
#                 S U B   R O U T I N E S 
#----------------------------------------------------------------------------------------
sub readInput{
# ------ storable on -------
      my $input_data = retrieve($runCalcInput);
      %run_data = %$input_data;
      my $FORM_data = retrieve($formInput);
      %FORM = %$FORM_data;
      $run_name = $run_data{run_name}; $WorkingDir = $run_data{WorkingDir}; $WWWdir = $run_data{WWWdir}; $epsilonPrecision = $run_data{epsilonPrecision}; $query_seq_name_to_run = $run_data{query_seq_name_to_run}; $optimizeBL = $run_data{optimizeBL}; $querySeqFoundinMSA = $run_data{querySeqFoundinMSA}; $tree_faq = $run_data{tree_faq}; $method = $run_data{method}; $SysErrorDef = $run_data{SysErrorDef}; $ErrorDef = $run_data{ErrorDef}; $ContactDef = $run_data{ContactDef}; $fileDna_aligned = $run_data{fileDna_aligned}; $treeUpload = $run_data{treeUpload}; $OutHtmlFile = $run_data{OutHtmlFile}; $cgi_log_file = $run_data{cgi_log_file}; $OutLogFile = $run_data{OutLogFile}; $fileName_amino_aligned = $run_data{fileName_amino_aligned}; $OutputURL = $run_data{OutputURL}; $estimated_run_time = $run_data{estimated_run_time}; $sequences_names_file = $run_data{sequences_names_file}; $was_Pdb_uploaded = $run_data{was_Pdb_uploaded};
      
      ($run_data{upload_TREE_file} eq "NOT_GIVEN") ? $upload_TREE_file = "" : $upload_TREE_file = $run_data{upload_TREE_file};
      ($run_data{PdbFileNameUnc} eq "NOT_GIVEN") ? $PdbFileNameUnc = "" : $PdbFileNameUnc = $run_data{PdbFileNameUnc};
($run_data{PdbPrefix} eq "NOT_GIVEN") ? $PdbPrefix = "" : $PdbPrefix = $run_data{PdbPrefix};
($run_data{pdb_data} eq "NOT_GIVEN") ? $pdb_data = "": $pdb_data = $run_data{pdb_data};
($run_data{clustal_aligned_file} eq "NOT_GIVEN") ? $clustal_aligned_file = "" : $clustal_aligned_file = $run_data{clustal_aligned_file};
      
      $upload_MSA_file_dna = $run_data{upload_MSA_file_dna};
      $upload_unaligned_file_dna = $run_data{upload_unaligned_file_dna};
# ------ storable on -------
# ------ storable off -------      
            #unless(open INPUT, $runCalcInput){
            #      open ANS, ">".$qsub_ans_file;
            #      print ANS "NOT_OK";
            #      close ANS;
            #      chmod 0755, $qsub_ans_file;
            #      exit;                  
            #}
            #while(<INPUT>){
            #            chomp;
                        #if(/RUN NAME: (.+)/) {$run_name = $1;}
                        #elsif(/WORKING DIR: (.+)/) {$WorkingDir = $1;}
                        #elsif(/WWW DIR: (.+)/){$WWWdir = $1;}
                        #elsif(/PRECISION LEVEL: (.+)/) {$epsilonPrecision = $1;}
                        #elsif(/EVOLUTONARY MODEL: (.+)/){$FORM{MODEL} = $1;}
                        #elsif(/EMPIRICAL MATRIX: (.+)/) {#ONLY IF IT IS MEC MODEL
                              #($1 eq "NOT_GIVEN") ? $FORM{EMPIRICAL_MATRIX} = "" : $FORM{EMPIRICAL_MATRIX} = $1;} 
                        #elsif(/QUERY NAME TO RUN: (.+)/){$query_seq_name_to_run = $1;}
                        #elsif(/DISTRIBUTE CATEGORIES: (.+)/){$FORM{CATEGORIES} = $1;}
                        #elsif(/TREE_WAS_UPLOADED\?: (.+)/){#REMARK: CHANGE THIS VAR'S CONTENT IN THE REST OF THE SCRIPT TO TRUE/FALSE.
                              #($1 eq "NOT_GIVEN") ? $upload_TREE_file = "" : $upload_TREE_file = $1;}
                              #$upload_TREE_file = $1;}
                        #elsif(/OPTIMIZE BRANCH LENGTH\? (.+)/){$optimizeBL = $1;}
                        #elsif(/GENETIC CODE: (.+)/){$FORM{GENCODE} = $1;}
                        #elsif(/GIVEN QUERY NAME: (.+)/){$FORM{msa_SEQNAME} = $1;}
                        #elsif(/PDB ID: (.+)/){#if given
                        #      ($1 eq "NOT_GIVEN") ? $FORM{pdb_ID} = "" : $FORM{pdb_ID} = $1;}
                        #elsif(/PDB NAME: (.+)/){
                        #      ($1 eq "NOT_GIVEN") ? $PdbFileNameUnc = "" : $PdbFileNameUnc = $1;}
                        ##elsif(/PDB CHAIN: (.+)/){
                        ##      ($1 eq "NOT_GIVEN") ? $FORM{chain} = "" : $FORM{chain} = $1;} 
                        #elsif(/PDB PREFIX: (.+)/){
                        #      ($1 eq "NOT_GIVEN") ? $PdbPrefix = "" : $PdbPrefix = $1;}
                        #elsif(/PDB DATA FILE: (.+)/){
                        #      ($1 eq "NOT_GIVEN") ? $pdb_data = "" : $pdb_data = $1;}
                        #elsif(/CLUSTAL ALN: (.+)/){
                        #      ($1 eq "NOT_GIVEN") ? $clustal_aligned_file = "" : $clustal_aligned_file = $1;}                        #elsif(/FOUND QUERY IN MSA\?: (.+)/){$querySeqFoundinMSA = $1;}
                        #elsif(/TREE FAQ: (.+)/){$tree_faq = $1;}
                        #elsif(/METHOD: (.+)/){$method = $1;}                        
                        #elsif(/USER EMAIL: (.+)/){
                        #      ($1 eq "NOT_GIVEN") ?  $FORM{email_address} = "" : $FORM{email_address} = $1;}                        
                        #elsif(/SYS ERROR: (.+)/){$SysErrorDef = $1;}
                        #elsif(/ERROR DEF: (.+)/){$ErrorDef = $1;}
                        #elsif(/CONTACT DEFINITION: (.+)/){$ContactDef = $1;}                        
                        #elsif(/WAS PDB UPLOADED\?: (.+)/){$was_Pdb_uploaded = $1;}
                        #elsif(/DNA FILE NAME: (.+)/){$fileDna_aligned = $1;}
                        #elsif(/UPLOADED TREE PATH: (.+)/){$treeUpload = $1;}
                        #elsif(/OUTPUT HTML PATH: (.+)/){$OutHtmlFile = $1;}
                        #elsif(/LOG PATH: (.+)/){$cgi_log_file = $1;}
                        #elsif(/QSUB LOG: (.+)/){$OutLogFile = $1;}
                        #elsif(/SEQ NAMES FILE: (.+)/){$sequences_names_file = $1;}
                        #elsif(/AMINO FILE NAME: (.+)/){$fileName_amino_aligned = $1;}
                        #elsif(/WAS A DNA UNALIGNED FILE UPLOADED\?: (.+)/){$upload_unaligned_file_dna = $1;}
                        #elsif(/WAS A DNA ALIGNED FILE UPLOADED\?: (.+)/){$upload_MSA_file_dna = $1;}                        #elsif(/URL OUTPUT: (.+)/){$OutputURL = $1;}
                        #elsif(/ESTIM RUNTIME: (.+)/){$estimated_run_time= $1;}
            #}
            #close INPUT;
            
# ------ storable off -------
            # if reading was OK, we report it, for the daemon
            open ANS, ">".$qsub_ans_file;
            print ANS "OK";
            close ANS;
            chmod 0755, $qsub_ans_file;
            
            # recreating the sequences array of the DNA sequences names
            unless (open SEQ_NAMES, $WorkingDir.$sequences_names_file){
                  &sys_error_exit("cannot open the file ".$WorkingDir.$sequences_names_file." for reading $!\n");}
            while(<SEQ_NAMES>){
                  chomp;
                  $sequences_names[0] = "";
                  if(/(\d+) (.+)/){
                        $sequences_names[$1] = $2;
                  }
            }
            close SEQ_NAMES;
}
#########################################################################################
# CALCULATION AND POST-PROCESSING
sub run_calc {
    &run_kaks4site;
	# The next routine assumes these files were created, therefore first we check that it was created
    unless ((-e $selection4Site_file) && !(-z $selection4Site_file)
            && (-e $kaks_file) && !(-z $kaks_file))
    {
      &sys_error_exit("run_calc: The file $selection4Site_file or $kaks_file was not created (or contains no data) during the run of selecton.v2.2. Cannot process outputs");
    }
    #print LOG "run_calc : going to run routine change_colors_if_significant\n";
    #&change_colors_if_significant;
    if ($was_Pdb_uploaded eq "yes") {  #if user ran selecton with a PDB struct.
	    ## run the script colorCoding.pl to produce the output files
	    print LOG "run_calc : touch $final_out\n";
		$proc_comm = "perl $colorCoding $method \'$query_seq_name_to_run\' $WorkingDir $selection4Site_file $clustal_aligned_file $pdb_data $kaks_file $final_out $PdbFileNameUnc $fileName_amino_aligned $rsml $params";
        print LOG "run_calc: running $proc_comm\n";
        system 'echo "(cd '.$WorkingDir.';touch '.$final_out.'; chmod oug+rx '.$final_out.')" | /bin/tcsh';
		system 'echo "(cd '.$WorkingDir.'; '.$proc_comm.')" | /bin/tcsh';	    
	    # check if the script $colorCoding found an error
	    if (-e $WorkingDir."error"){
            &read_colors_error_and_exit;
        }		    
        # REMARK: I don't think it is necesseray in the new server, since we wont use PE
	    ##### copy final PE files to pdbspt dir and compress the PDB file 
	    #my $string1 = "cd $WorkingDir";
	    #my $string2 = "cp consurf.spt  pdbspt/consurf.spt";
	    #my $string3 = "gzip -c $PdbFileNameUnc > pdbspt/pdbfile.ent";
	    #my $string4 = "mv consurf.spt colors.txt";
	    #my $string6 = "mv $rasmol_file rasmol.txt"; 
	    #my $string5 = "chmod  ogu+rx  pdbspt/*";
	    #
	    #print LOG "\nrun_calc: Copy the final PE files to pdbspt dir\n";
	    #
	    #system 'echo "('.$string1.'; '.$string2.'; '.$string3.'; '.$string4.'; '.$string6.';' .$string5.';)" | /bin/tcsh';
	}

	$proc_comm = "perl $colorCodingLinear $run_name $WorkingDir $colorsLinear $selection4Site_file $areTherePositiveSites";    
	print LOG "\nrun_calc: running $proc_comm \n";	
	system 'echo "(cd '.$WorkingDir.'; '.$proc_comm.')" | /bin/tcsh';
    
    # check if the script $colorCodingLinear found an error
    if (-e $WorkingDir."error"){
            &read_colors_error_and_exit;
    }	
}

######################################################################################
# run kaks4site.exe 
sub run_kaks4site {
###### should add verification that two refuting arguments aren't given here... 
#  	my $selecton_comm="$selecton -c \'$fileDna_aligned\'"; #default run:  bayesian, beta+w>1 (M8) , w=8, ref=1st seq, std nuc.code, NJ tree
    
my $selecton_comm="$selecton -i $WorkingDir" . $fileDna_aligned . " -e" . $epsilonPrecision; #default run:  bayesian, beta+w>1 (M8) , w=8, ref=1st seq, std nuc.code, NJ treeepsilon by default set to 0.1
    
    if ($FORM{MODEL} eq "MEC") { #if the model is MEC call another exe
                        $selecton_comm  = "$mecSelecton -i $WorkingDir" . "$fileDna_aligned"; 
	   	if ($FORM{EMPIRICAL_MATRIX} eq "JTT"){
	   		$selecton_comm .= " -z 0";
		}
		if ($FORM{EMPIRICAL_MATRIX} eq "WAG"){
	   		$selecton_comm .= " -z 1";
		}
		if ($FORM{EMPIRICAL_MATRIX} eq "mtREV24"){
	   		$selecton_comm .= " -z 2";
		}
		if ($FORM{EMPIRICAL_MATRIX} eq "cpREV45"){
	   		$selecton_comm .= " -z 3"; # if no z is given mecSelecton will run by default with JTT
		}
	}
    if ($querySeqFoundinMSA eq "yes") {
	 	$selecton_comm .= " -q \'$query_seq_name_to_run\'";  
    }
   	if ($FORM{MODEL} eq "M7") { # beta no additional omega1. prob(beta) set to 1
	   	$selecton_comm .= " -p1 -Fp";
   	}
   	if ($FORM{MODEL} eq "M8a") { #beta + w = 1
	   	$selecton_comm .= " -w1 -Fw"; #do not optimize omega (omega set to 1) (M8a)
   	}
   	if ($FORM{MODEL} eq "M5") {
	   	$selecton_comm .= " -dg";
   	}
   	if ($FORM{CATEGORIES} ne "") { 
	   	my $catAdd=" -n ".$FORM{CATEGORIES};
	   	$selecton_comm .= $catAdd;
   	}
   	if ($upload_TREE_file ne "NOT_GIVEN") {
    	$selecton_comm .= " -u \'$treeUpload\'";
	}  	
   	if ($optimizeBL eq "n"){
	   $selecton_comm .= " -bn";	
   	}
   	if ($FORM{GENCODE} != 0) {
	   	my $genAdd=" -g ".$FORM{GENCODE};
	   	$selecton_comm .= $genAdd;	
   	}

 
    print LOG "\nrun_kaks4site: running $selecton_comm\n";
    print LOG "run_kaks4site: SeqName = ***$query_seq_name_to_run***\n";
    
    system "cd $WorkingDir; $selecton_comm; chmod  ogu+rx *"; # The program can't run from rsh
    
    #check for user errors in kaks4ite.log
   my $kaksLogFile = $WorkingDir.$log; 
   open OUTPUT, ">>$OutHtmlFile";
   unless (open (LOGFILE,"$kaksLogFile")) {
	        close (LOGFILE);
      &sys_error_exit("Error in run_kaks4site, $kaksLogFile does not exist");
   }
   while (<LOGFILE>){
       my $line=$_;
       if ($line =~ /\S+/){
            my @userError   = split(/\s/,$line);
            if ($userError[0] eq "USER"){
                 $line =~ s/USER ERROR://; ## query sequence not found
                 my $first_line  = $line;
                 while (<LOGFILE>){
                     $line = $_;
                     $first_line = "<br>".$first_line."<br>".$line;
                 }
                 print OUTPUT "\n<p><ul><li><font color='red'>Warning:</b></font> The query sequence name \'$FORM{msa_SEQNAME}\' is not found in the <A HREF=$fileDna_aligned TARGET=MSA_window>MSA file</A>.<br>The calculation continues. The first sequence in MSA is used as a query.</li></ul></p>\n";
                 close OUTPUT;
                 print LOG "\nrun_kaks4site: query sequence not found. Calculation continues with 1st sequence in MSA.\n";
            }
            if (($line =~ /found in the tree file but not found in the sequence file/) || ($line =~ /Error reading tree file/)) { #mismatch between MSA names and tree names
                 my $err=$line;
                 my $line1 = <LOGFILE>;
                 my $line2 = <LOGFILE>;
                 $err .= "$line1 "."$line2";
                 close (LOGFILE);
                 &print_to_output_and_exit("Error in tree file:<br>$err Please check that all the names in the sequence file are identical to all the names in the tree.", "run_kaks4site: $err");                      		
            }
            elsif ($line =~ /Bad format in tree file/){
                 close (LOGFILE);
                 &print_to_output_and_exit("<b>Bad format in tree file.</b><br>Please correct your tree file according to <a href=\"$tree_faq\">Selecton accepted format</a> and re-submit your query.", "run_kaks4site: $line");        
            }
            elsif ($line =~ /The nucleotide sequences contained the character: (.*)/) {
            my $illegal = $1;
            close (LOGFILE);
            &print_to_output_and_exit("<b>The nucleotide sequences file contained an illegal character: $1. Only the following characters are accpted: A,C,G,T,-. Please correct your file and re-submit it to Selecton.","run_kaks4site: $line");
            }
            elsif($line =~ /Unable to read file. It is required that each line is no longer than/){
                  close (LOGFILE);
                  &print_to_output_and_exit("<b>Selecton does not accept DNA sequences which are longer than ".GENERAL_CONSTANTS::SELECTON_MAX_NUCLEOTIDE." nucleotides.</b>", "run_kaks4site: $line");
            }
        }
   }	
    close OUTPUT;
    close (LOGFILE);
}
######################################################################
# creating new files to hold DNASeqNames in DNA file and AMINO file
sub return_seq_names_to_files{
   
   my $current_file = shift;
   my $new_file = shift;
   my $ref_seq_name_arr = shift; # reference to the sequence names array
   
   unless (open IN, $WorkingDir.$current_file){
      print LOG "could not open file $WorkingDir"."$current_file for reading. Names of files will be displayed as numbers.\n" ;
   }
   else {
      unless (open OUT, ">".$WorkingDir .$new_file){
         print LOG "could not open file $WorkingDir"."$new_file for writing. Names of files will be displayed as numbers.\n";
      }
      else{
         while (<IN>){
            if(/>(\d+)/){
               print OUT ">".$ref_seq_name_arr->[$1]."\n";
            }
            else{
               print OUT $_;
            }
         }
         close OUT;
      }
      close IN;      
   }   
}
######################################################################################
sub change_tree_file_names{
   
   my $ref_sequences_names = shift;
   my $input_tree = shift;
   my $output_tree = shift;
   my ($tree, $err);
   
   print LOG "change_tree_file_names : going to change numbers from tree $input_tree to names in file $output_tree\n";
   
   unless (open TREE, $input_tree){      
         print LOG "change_tree_file_names : could not open the file $input_tree for reading. the tree file will be presented with numbers\n";
         return;
   }
   #check validity of input tree file
   $tree = <TREE>;   
   close TREE;
   
   
   my @tree_arr = split(/\(/, $tree);
   my @sub_tree = ();
   my @temp_arr;
   my $sub_counter = 0;
   
   # building the array @sub_tree, so that each cell will hold maximum one sequence name
   for(my $i=0; $i<@tree_arr; $i++){
      if ($tree_arr[$i] ne ""){
         $tree_arr[$i] = "(".$tree_arr[$i];
      }
      if ($tree_arr[$i] =~ m/.*,.+/){
         @temp_arr = split(/,/, $tree_arr[$i]);
         foreach (@temp_arr){
            $sub_tree[$sub_counter] = $_.",";
            $sub_counter++;
         }
      }
      else{
         $sub_tree[$sub_counter] = $tree_arr[$i];
         $sub_counter++;        
      }
   }

   # rebuilding the tree, this time replacing the sequences names with the names found in the DNA input file
   my $final_tree = "";
   my ($exp, $rest_of_exp, $new_rest_exp);
   my $seq_found = "no";
   for (my $k=1; $k<@sub_tree; $k++){
      #in this part we wish to split the expression to 2 parts; left part : (?seq_name ; right part: all the rest       
      if ($sub_tree[$k] ne ""){
         if ($sub_tree[$k] =~ m/(.+)(:.+)/){
            $exp = $sub_tree[$k];
            $rest_of_exp = "";
            while ($exp =~ m/(.+)(:.+)/){
               $exp = $1;
               $rest_of_exp = $2.$rest_of_exp;
            }
         }
           # in case the expression is of format:  seq_name:distance,
         elsif($sub_tree[$k] =~ m/(.+)(\);.+)/){
            $exp = $1;
            $rest_of_exp = $2;
            while ($exp =~ m/(.+)(\))/){
               $exp = $1;
               $rest_of_exp = $2.$rest_of_exp;
            }
         }
         #  in case the expression is of format:  seq_name)*,
         elsif($sub_tree[$k] =~ m/(.+)(\)?.+)/){
            $exp = $1;
            $rest_of_exp = $2;
            while ($exp =~ m/(.+)(\))/){
               $exp = $1;
               $rest_of_exp = $2.$rest_of_exp;
            }            
         }
         # if the length (value after the ":") is equal to zero, we replace it with a very small value,
         # because the selecton.exe cannot calculate trees with zeros
         $new_rest_exp = "";
         while($rest_of_exp =~ m/(.?:)(\d\.?\d*)(.+)/){
            if(!($2>0) && !($2<0)){
               $rest_of_exp = $3;
               $new_rest_exp .= $1."0.000000001";
            }
            else{
               $rest_of_exp = $3;
               $new_rest_exp .= $1.$2;
            }
         }
         $new_rest_exp .=$rest_of_exp;
         $rest_of_exp = $new_rest_exp;
         
         $exp =~ m/(\(?)(.+)/;                    
         $final_tree.= $1.$ref_sequences_names->[$2].$rest_of_exp;         
      }       
       #an empty cell stands for a "(" sign
      else{
         $final_tree.= "(";
      }
   }
   if ($final_tree =~ m/,$/){
      chop $final_tree;
   }
   
   unless (open NEW_TREE, ">".$output_tree){
      &sys_error_exit("change_tree_file_names:: cannot open file $output_tree for writing.");
   }
   print LOG "change_tree_file_names : printing edited tree to file $output_tree and chmod it.\n";
   print NEW_TREE $final_tree;
   close NEW_TREE;
   chmod 0755, $output_tree;   
}
######################################################################
# prepare the variables to be sent to the "pipe" script. than calls the script with all needed vars.
sub prepare_pipe{
   
   my $sequences_names = shift;
   
   print LOG "\nEentered prepare_pipe()\n";
   # vars that should be edits before sent to the pipe script:
   my ($pipe_pdb_id, $pipe_chain, @pipe_Model, $model_ref, $pipe_distribute_categories, $pipe_optimizeBL, @pipe_genetic, $genetic_ref, $pipe_dna_input, $pipe_query, $pipe_precision, $pipe_empirical);
   
   ($FORM{pdb_ID} eq "") ? $pipe_pdb_id = "UPLOADED" : $pipe_pdb_id = $FORM{pdb_ID};
   ($FORM{chain} eq "") ? $pipe_chain = "none" : $pipe_chain = $FORM{chain};
   
   if ($FORM{MODEL} eq "M8") {$pipe_Model[0] = 'Positive selection enabled (M8, beta + w >= 1)' ;}
   elsif ($FORM{MODEL} eq "M8a") {$pipe_Model[0] = 'Null model: no positive selection(M8a, beta + w = 1)' ;}
   elsif ($FORM{MODEL} eq "M7") {$pipe_Model[0] = 'Null model: no positive selection(M7, beta)' ;}
   elsif ($FORM{MODEL} eq "M5") {$pipe_Model[0] = 'Positive selection enabled(M5, gamma)' ;}
   elsif ($FORM{MODEL} eq "MEC") {$pipe_Model[0] = 'Mechanistic Empirical Combination Model (MEC)' ;}
   else {$pipe_Model[0] = "IGNORED";}
   $model_ref = \@pipe_Model;
   
   ($FORM{MODEL} eq "MEC") ? $pipe_empirical = $FORM{EMPIRICAL_MATRIX} : $pipe_empirical = "IRRELEVANT";
       
   if ($epsilonPrecision == 0.1) {$pipe_precision = "Intermediate precision";}
   elsif ($epsilonPrecision == 1) {$pipe_precision = "Low precision- faster run";}
   elsif ($epsilonPrecision == 0.01) {$pipe_precision = "High precision- slower run";}
   else {$pipe_precision = "DEFAULT";}
       
   ($FORM{CATEGORIES} eq "") ? $pipe_distribute_categories = 8 : $pipe_distribute_categories = $FORM{CATEGORIES};
   ($optimizeBL eq "y") ? $pipe_optimizeBL = "True" : $pipe_optimizeBL = "False";
   
   if ($FORM{GENCODE}==0) {$pipe_genetic[0] = "Nuclear Standard" ;}
   elsif ($FORM{GENCODE}==1) {$pipe_genetic[0] = "Nuclear Blepharisma";}
   elsif ($FORM{GENCODE}==2) {$pipe_genetic[0] = "Nuclear Ciliate";}
   elsif ($FORM{GENCODE}==3) {$pipe_genetic[0] = "Nuclear Euplotid";}
   elsif ($FORM{GENCODE}==4) {$pipe_genetic[0] = "Mitochondria Vertebrate";}
   elsif ($FORM{GENCODE}==5) {$pipe_genetic[0] = "Mitochondria Invertebrate";}
   elsif ($FORM{GENCODE}==6) {$pipe_genetic[0] = "Mitochondria Yeast";}
   elsif ($FORM{GENCODE}==7) {$pipe_genetic[0] = "Mitochondria Ascidian";}
   elsif ($FORM{GENCODE}==8) {$pipe_genetic[0] = "Mitochondria Echinoderm";}
   elsif ($FORM{GENCODE}==9) {$pipe_genetic[0] = "Mitochondria Flatworm";}
   elsif ($FORM{GENCODE}==10) {$pipe_genetic[0] = "Mitochondria Protozoan";}
   else {$pipe_genetic[0] = "IGNORED";}
   $genetic_ref = \@pipe_genetic;
      
   ($sequences_names->[$query_seq_name_to_run] eq "") ? $pipe_query = "\"\"" : $pipe_query = $sequences_names->[$query_seq_name_to_run];
   
   #since there is only 1 input file, there will be only var sent to the pipe script. In order that the pipe script will know what kind of input it is - a short string is added at the beginning.
   if ( $upload_unaligned_file_dna ne "no"){
      $pipe_dna_input = "SELECTON_UN".$upload_unaligned_file_dna;
   }
   elsif ($upload_MSA_file_dna ne "no"){
      $pipe_dna_input = "SELECTON_MSA".$upload_MSA_file_dna;
   }
   else{
      $pipe_dna_input = "no_dna_input";
   }
   
   open ERROR, $WorkingDir.$pipe_error;
   close ERROR;
   chmod 0755, $WorkingDir.$pipe_error;
   
   print LOG "running pipe_for_selecton::create_pipe with parameters:\n";
   print LOG "$run_name $cgi_log_file $WorkingDir $rsml $FGiJ_pipe_pdb $pdbUpload $pipe_pdb_id $pipe_chain $pipe_dna_input $pipe_query $model_ref $pipe_distribute_categories $pipe_optimizeBL $genetic_ref $pipe_error $pipe_empirical $pipe_precision\n";
   
   pipe_for_selecton::create_pipe($run_name, $cgi_log_file, $WorkingDir, $rsml, $FGiJ_pipe_pdb, $pdbUpload, $pipe_pdb_id, $pipe_chain, $pipe_dna_input, $pipe_query, $model_ref, $pipe_distribute_categories, $pipe_optimizeBL, $genetic_ref, $pipe_error, $pipe_precision, $pipe_empirical);
   
   # checking if there was an error and the pipe file was not created properly
   if (-e $WorkingDir.$pipe_error && !(-z $WorkingDir.$pipe_error)){
      unless (open ERROR, $WorkingDir.$pipe_error){
         &sys_error_exit("An error was found when trying to create the pipe file for Selecton.\nThe Error message should be written to file $WorkingDir"."$pipe_error, however this file could not be opened.\n");
      }
      &sys_error_exit("An error was found while trying to create the pipe file : ".<ERROR>);
   }   
}
######################################################################
sub sys_error_exit{
   my $err = shift;
   
   open OUTPUT, ">>$OutHtmlFile";
   print OUTPUT $SysErrorDef;
   print OUTPUT $ContactDef;
   close OUTPUT;
   print LOG "\n$err\n";
   &send_mail();
   &send_mailSelecton("SYSTEM ERROR\n".$err);
   &stop_reload;
   exit;   
}

##########################################################################################
# Stops the reload of the output page
sub stop_reload {
   sleep 5;
   open OUTPUT, "<$OutHtmlFile";
   flock OUTPUT,2;
   my @output = <OUTPUT>;
   flock OUTPUT,8;
   close OUTPUT;   
   open OUTPUT, ">$OutHtmlFile";
   flock OUTPUT,2;
   foreach my $line (@output){  # we remove the refresh lines and the button which codes for Selecton cancelled job
      unless ($line =~ /REFRESH/ or $line =~ /NO-CACHE/ or $line =~ /ACTION=\"cgi.+kill/ or
      $line =~ /VALUE="Cancel Selecton Job"/ or $line =~/TYPE=hidden NAME="Qstat_file"/ or $line =~/TYPE=hidden NAME="selecton_http"/ or $line =~ /TYPE=hidden NAME="run_no"/ or $line =~ /TYPE=hidden NAME="cgi_pid"/ or $line =~ /Estimated run time is:/ or $line =~ /kill_process.cgi/ or $line =~ /<!--job_/){
            print OUTPUT $line;
      }
   }
   flock OUTPUT,8;
   close OUTPUT;
   print LOG "\n\nEnd time: ";
   print LOG &printTime();
   close LOG;
   
    # remove the job from the running jobs list
    &BIOSEQUENCE_FUNCTIONS::remove_job_from_running_log("Selecton", $run_name);
        
    open FINISH, ">".$finish_flag;
    close FINISH;
    unlink $qsub_ans_file if (-e $qsub_ans_file);        
    
    chmod 0755, $OutHtmlFile;
    chmod 0600, $WorkingDir."user_email.txt";
    chmod 0711, $WorkingDir;
    chmod 0600, $runCalcInput;
#    if (-e $WorkingDir."core"){
#	print LOG "remove core file from working directory\n";
#	unlink $WorkingDir."core";
#    }
}

#########################################################################################
# Sends an automatic mail when there are errors
sub send_mail { # to user
   if ($FORM{email_address} ne ""){
      $email_subject = "Error in Selecton running";
      $email_message = "Hello!\n\nUnfortunately there was an error while running Selecton.\nPlease click on the following link to see more details\nWe apologize for the inconvenience\n\n$OutputURL";   
      print LOG "send_mail: sending system error to user:\n";
      my $mail_line = 'perl sendEmail.pl -f \''.GENERAL_CONSTANTS::ADMIN_EMAIL.'\' -t '.$FORM{email_address}.' -u \''.$email_subject.'\' -xu '.$userName.' -xp '.$userPass.' -s '.$smtp_server.' -m \''.$email_message.'\'';
      chdir $send_email_dir;
      $email_system_return = `$mail_line`;
      unless ($email_system_return =~ /successfully/)    {
         print LOG "The message was not sent successfully. system returned: $email_system_return\n";
      }
   }
}
#########################################################################################
sub send_mailSelecton{ # to selecton administrator
   my $email_message = shift;
   if ($FORM{email_address} eq ""){$FORM{email_address} = "NOT_GIVEN";}
   $email_subject = "Error in Selecton running $run_name";
   print LOG "send_mailSelecton: send error message to admin:\n";   
   my $mail_line = 'perl sendEmail.pl -f \''.GENERAL_CONSTANTS::ADMIN_EMAIL.'\' -t '.GENERAL_CONSTANTS::ADMIN_EMAIL.' -u \''.$email_subject.'\' -xu '.$userName.' -xp '.$userPass.' -s '.$smtp_server.' -m \''.$email_message.'\nUser email is:  '.$FORM{email_address}.'\'';
   chdir $send_email_dir;
   $email_system_return = `$mail_line`;
   unless ($email_system_return =~ /successfully/)    {
      print LOG "The message was not sent successfully. system returned: $email_system_return\n";
   }
}
#########################################################################################
sub printTime {   
   my $theTime = BIOSEQUENCE_FUNCTIONS::printTime;
   return $theTime;
}
#########################################################################################
sub print_to_output_and_exit{
   my $html_err = shift;
   my $log_err = shift;

   open OUTPUT, ">>$OutHtmlFile";
   print OUTPUT "\n<p>$ErrorDef<br>$html_err</p>\n";
   print OUTPUT $ContactDef;
   close OUTPUT;	    
   print LOG "$log_err";
   &send_mail();
   &stop_reload;
   exit;
}
#########################################################################################
# in case an error file was found during the run of either $colorCoding or $colorCodingLinear
sub read_colors_error_and_exit{
            unless (open ERROR, $WorkingDir."error"){
                        &sys_error_exit("run_calc : An error was found after running $colorCoding, but the error file $WorkingDir"."error could not be opened.");
            }
            my ($h_err, $l_err);
            while (<ERROR>){
                  chomp;
                  if (/HTML: (.+)/){
                        $h_err = $1;}
                  elsif(/LOG: (.+)/){
                        $l_err = $1;}                        
            }
            close ERROR;
            ($h_err =~ m/^sys$/) ? &sys_error_exit("run_calc : $l_err") : &print_to_output_and_exit($h_err, $l_err);        
}
#########################################################################################
sub add_data_to_input_file{
# ---- storable on ----
      $run_data{params} = $params;
      $run_data{log} = $log;
      $run_data{outputScoreFile} = $outputScoreFile;
      $run_data{FORMInput} = $formInput;
      # recreating the storable element
      unlink $runCalcInput;
      store \%run_data, $runCalcInput;
      if (!(-e $runCalcInput) or (-z $runCalcInput)){
            &sys_error_exit("Seems that the store filed. the file $runCalcInput doesn't exists or is of size zero\n");
      }
      chmod 0600, $runCalcInput;
# ---- storable on ----
# ---- storable off ----
      #unless (open INPUT, ">>".$runCalcInput){
      #      open LOG, ">>$OutLogFile";
      #      print LOG "add_data_to_input_file : cannot open file $runCalcInput for writing $!. The statistical testing button cannot be active\n";
      #      close LOG;
      #      return "no";                        
      #}
      #print INPUT "GLOBAL RESULTS FILE: $params\n";
      #print INPUT "KAKS LOG FILE: $log\n";
      #print INPUT "SCORE RESULTS FILE: $outputScoreFile\n";
      #close INPUT;
# ---- storable off ----      
      return "OK";
}
######################################################################
sub extract_file_name{
   my $file = shift;
   $file =~ s/\s/_/;
   return $file;   
}
