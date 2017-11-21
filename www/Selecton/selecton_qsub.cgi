#!/usr/bin/perl 


use CGI;
use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);

use strict;
use Storable;
use Bio::SeqIO;
#use Bio::TreeIO;
use Bio::Tree::NodeI;
use Bio::Root::Root;
use Bio::Tree::TreeI;
use lib "/bioseq/bioSequence_scripts_and_constants";
use GENERAL_CONSTANTS;
use SELECTON_CONSTANTS;
use BIOSEQUENCE_FUNCTIONS;
use TREE_parser;
use lib "/bioseq/Selecton/external_scripts";
use codonAlign;



###### READING DATA FROM FORM - FIRST TO DO

my $queryForm = new CGI;
my %FORM;        #hash with form information


#***********************************
### input file, ref seq. and email address

my $upload_unaligned_file_dna = $queryForm->param('userFILEunaligned');
my $upload_MSA_file_dna=$queryForm->param('userFILEaligned');
$FORM{msa_SEQNAME} = $queryForm->param("msa_SEQNAME");
$FORM{email_address} = $queryForm->param("email_add");
my $recipient = $FORM{email_address};

# check if the user which runs this run has not exeeded its maximal number of runs
my $user_ip = $ENV{'REMOTE_ADDR'};
BIOSEQUENCE_FUNCTIONS::check_if_user_is_allowed("selecton",$user_ip, $recipient);

#***************************************
### advanced options
$FORM{pdb_ID} = $queryForm->param("pdb_ID");
my $upload_PDB_file = $queryForm->param("pdb_FILE");
my $pdbUploadName = "FILE";
$FORM{chain} = $queryForm->param("chain");
$FORM{MODEL} = $queryForm->param("MODEL");
$FORM{EMPIRICAL_MATRIX} = $queryForm->param("EMPIRICAL_MATRIX"); #adid added the mec model
my $empiricalMatrix = $FORM{EMPIRICAL_MATRIX};
my $epsilonPrecision = $queryForm->param("PRECISION"); # added v2.2  # NOTE: if more genetic codes are added (or changed), the info should be edit in the "prepare_pipe" routine.
my $model=$FORM{MODEL}; # NOTE: if more models are added (or changed), the info should be edit in 2 more places, as text: 1. in the OUTPUT file, in "start_output_html" routine, 2. in the "prepare_pipe" routine.
my $method="Bayesian"; ## ugly patch, since we removed ML: DO NOT REMOVE THIS LINE, it is necessary for colorCoding
#$FORM{DISTRIBUTION} = $queryForm->param("DISTRIBUTION");
$FORM{CATEGORIES} = $queryForm->param("CATEGORIES");
my $upload_TREE_file = $queryForm->param("tree_FILE");
$FORM{BL} = $queryForm->param("BL");
my $optimizeBL="y"; #true if checked
if ($FORM{BL} eq ""){
   $optimizeBL="n";
}
$FORM{GENCODE} = $queryForm->param("GENCODE"); # NOTE: if more genetic codes are added (or changed), the info should be edit in the "prepare_pipe" routine.

#****************************************
#### TRANSLATE to LOWER case the PDB file and UPPER case the CHAIN

$FORM{pdb_ID} =~ tr/[A-Z]/[a-z]/;
$FORM{chain} =~ tr/[a-z]/[A-Z]/;


#****************************************
#### set pdb_ID to FILE if uploading a file

if ($upload_PDB_file ne "") {   
   $FORM{pdb_ID} = $pdbUploadName;
}

#**************************************
##### variables #######
my $querySeqFoundinMSA="no";
my $query=$FORM{msa_SEQNAME};

#**********************************
##### GENERAL PATHS 
my $run_name = $^T;  #the running dir NAME old $$
my $ibis_external_scripts_path = "/bioseq/bioSequence_scripts_and_constants/";
my $WorkingDir = GENERAL_CONSTANTS::SERVERS_RESULTS_DIR."Selecton/" . $run_name . "/";
my $http_path = GENERAL_CONSTANTS::SELECTON_URL;
my $job_canceled_page = $http_path.'/cancel_page.html';
my $WWWdir = $http_path."/results/" . $run_name . "/";
my $PdbPath = GENERAL_CONSTANTS::PDB_DIVIDED;

#***************************************
### PROGRAMS OR PERL OR EXECUTABLES 
my $extractPDBinfo = $ibis_external_scripts_path . "extract_info_from_pdb.pl";
my $kill_job_script = "/cgi-bin/kill_process.cgi";
my $qsub_script = $WorkingDir."qsub.sh";
my $runClac_inQ = "/bioseq/Selecton/selecton_run_calc.pl";
my $clustalw= 'ssh bioseq@biocluster clustalw'; #'/usr/local/bin/clustalw';
my $muscle = 'ssh bioseq@biocluster muscle'; #'/usr/local/bin/muscle';


#***************************************
### Sending e-mail from ibis
my $send_email_dir = GENERAL_CONSTANTS::SEND_EMAIL_DIR_IBIS;
my $smtp_server = GENERAL_CONSTANTS::SMTP_SERVER;
my $userName = GENERAL_CONSTANTS::ADMIN_USER_NAME;
my $userPass = GENERAL_CONSTANTS::ADMIN_PASSWORD;
my $email_subject;
my $email_message;
my $email_system_return;

#***************************************
### output related paths
my $InpSeqFile = $WorkingDir . "path.txt"; #file containing ls result of pdb 
my $OutputURL = $WWWdir  ."output.html"; #link to output file
my $OutHtmlFile = $WorkingDir . "output.html"; #OUTPUT to the user
my $Logs_dir = GENERAL_CONSTANTS::SERVERS_LOGS_DIR."Selecton/";
my $OutLogFile = $Logs_dir.$run_name.".log";
my $QsubLogFile = $Logs_dir.$run_name."_Q.log";

#**************************************
##### SPECIFIC TO THIS RUN variables
my $pid;  #to set the pid of the child
my $PdbFileDir = ""; #pdb file specific dir on path (2 letters) to check
my ($cgi_pid, $cmd);
my @sequences_names = (); # This array will hold the names of the sequences, as it appears in the input file.
my $Qstat_No_file = "QSTAT_NO";
my $estimated_run_time = "none";

#***************************** 
### making absolute path file for ATEN pdb1ed5.ent.Z or pdbFILE.ent.Z if file was uploaded
my $PdbFileName = "pdb" .$FORM{pdb_ID}. ".ent.gz"; 
my $PdbFileNameUnc = "pdb" .$FORM{pdb_ID}. ".ent"; 
my $PdbPrefix = "pdb" .$FORM{pdb_ID};
if ($FORM{pdb_ID} eq "") {
	$PdbPrefix="";
}
# file names in use only in case of PDB reading, for the use of the script $extractPDBinfo
my ($pdb_data) = ($PdbFileNameUnc) =~ /(\w+)/; #extracting prefix, put it in first var
my $title_file = $pdb_data.".title";  
my $pdb_fasta = $pdb_data.".pdbfasta";
$pdb_data.=".pdbdata";
my $pdb_to_fasta_error = "pdb_to_fasta.error";
my $pdb_msa = $PdbPrefix . "_PDB_MSA.pdbfasta";
my $clustal_outFile = $PdbPrefix . "_PDB_MSA.out"; # the outfile for clusalw
my $clustal_aligned_file = $PdbPrefix."_PDB_MSA.aln";

#***************************** 
### FILE NAMES IN USE
my $dnaMSAprefix= $PdbPrefix ."DNA";
my $aminoMSAprefix= $PdbPrefix ."AMINO";
my $fileUploadName_dna_unaligned=$dnaMSAprefix."_unaligned".".txt";
my $fileDna_aligned = $dnaMSAprefix.".msa";
my $fileName_amino_aligned = $aminoMSAprefix.".msa";
my $fileUploadPath_dna_unaligned=$WorkingDir . $fileUploadName_dna_unaligned; #uploaded file from user - unaligned
my $copied_dna_unaligned = $WorkingDir . "COPY_".$fileUploadName_dna_unaligned;;  # a copy of the same file
my $wwwfileUploadPath_dna_unaligned = $WWWdir . $fileUploadName_dna_unaligned;#www path : uploaded file from user - unaligned
my $fileDnaPath_aligned = $WorkingDir . $fileDna_aligned;  #dna file after alignment
my $fileAminoPath_aligned = $WorkingDir . $fileName_amino_aligned;  #amino file after alignment
($PdbFileDir) = ($FORM{pdb_ID}) =~ /\w(\w{2})/;
my $PdbFilePath =  $PdbPath . $PdbFileDir . "/" . $PdbFileName; 
my $codonAlignLogFile=$WorkingDir."codonAlign.log";
my $treeUpload = $WorkingDir . "userTree.txt";
my $copied_treeUpload = $WorkingDir . "COPY_userTree.txt";  # a copy of the same file
my $userTree = $WWWdir . "userTree.txt";
my $userTree_copy = $WWWdir . "COPY_userTree.txt";
my $PdbFile = $PdbPrefix . ".ent";
my $sequences_names_file = "sequences.names";
my $runCalcInput = "runCalcInput.txt";
my $FormInput = "formInput.txt";
my $qsub_ans = $WorkingDir."qsub_ans.txt";
my $statistics_file = SELECTON_CONSTANTS::STATISTICS_FILE;

#***************************** 
### pdb related links
my $pdbUpload = $WorkingDir . $PdbFileNameUnc; #name for an uploaded PDB FILE

#**********************************
### HTML definitions
my $ErrorDef = "<font size=+3 color='red'>ERROR!  Selecton session has been terminated: </font>\n";
my $SysErrorDef = "<p><font size=+3 color='red'>SYSTEM ERROR - Selecton session has been terminated!</font><br><b>Please verify that there are no errors in your input file/s, and try to run Selecton again. Specifically, make sure your file is in the <a href=\"$http_path/faq.html#q4\">correct format</a>, and refer to the <a href=\"$http_path/faq.html\">FAQ</a> for further assistance.</b></p>\n";
my $tree_faq = $http_path.'/faq.html#q8';
my $SystemError = "<b>A system error occured during the calculation. Please try to run Selecton again in a few minutes.</b>\n";
my $mail = "mailto:".GENERAL_CONSTANTS::ADMIN_EMAIL."?subject=Selecton%20Run%20No:%20$run_name";
my $ContactDef = "<H3><center>For assistance please <a href=\"$mail\">contact us</a> and mention this number: $run_name</H3>\n";
my $QuickHelp = 'http://consurf.tau.ac.il/quick_helpver3.html#chain';  # REMARK : ugly, should change
my $status_faq = $http_path."/faq.html#q13";

my $reload_interval = 30;




#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#***************************************************
#### MAIN CGI THAT FORKS AND USES THE FILES TO RUN
#***************************************************

#********** FATHER **********
##### TO BE DONE BEFORE LONG PROCESS
if ($pid = fork) {
   exit;
} 


#******************************
#### CHILD ###
elsif (defined $pid) {
   # create the WorkingDir
   mkdir  $WorkingDir;
   system 'echo "(touch '.$OutLogFile.'; chmod 750 '.$OutLogFile.')" | /bin/tcsh';
   open LOG, ">$OutLogFile";   
   print LOG "\n************** LOG FILE *****************\n\n";
   print LOG "\nBegin time: ";
   my $curr_time = &printTime();   
   print LOG "\nuser's email is: $recipient\n";
   
   # add this run to the log
    open LIST, ">>".GENERAL_CONSTANTS::SELECTON_LOG;
    flock LIST, 2;
    print LIST $curr_time." $run_name $user_ip $recipient\n";
    flock LIST, 8;
    close LIST;
# print the user's ip and run name in the running jobs list

    open RUN_LIST, ">>".GENERAL_CONSTANTS::SELECTON_RUNNING_JOBS;
    flock RUN_LIST, 2;
    print RUN_LIST "$run_name $user_ip $recipient\n";
    flock RUN_LIST, 8;
    close RUN_LIST;      

   
   #********************************************
   ##### CREATE AND START file for output UPDATE
   &start_output_html;
   
   ### Move directly to the output file
   print "Location: $OutputURL\n\n";
   # unmark comment if you want to print to the screen a message, instead of a new location:
   #print header;
   # print start_html;
   # print h3("Due to a temporal error in Selecton server, the results page can not be viewed via the browser. Please save the link to your results, as soon as the error is fixed - you will be able to view it in this url:<br> $OutputURL<br><br>We deeply apologize for the inconvenience.<br><br><br><br>");
   # print $ContactDef;
   # print end_html;
   #**************************************
   ##### disconnecting CHILD flush buffers  
   close(STDIN); 
   close(STDOUT); 
   close(STDERR);
   ### if PDB ID: checking if PDB FILE exists on aten pdb DB
   # send user's e-mail to a special file
   open USER_MAIL, ">$WorkingDir" . "user_email.txt";
   print USER_MAIL $FORM{email_address};
   close USER_MAIL;
   chmod 0600, $WorkingDir. "user_email.txt";      
   
   if ($FORM{pdb_ID} ne "" and $FORM{pdb_ID} ne $pdbUploadName) {
      &check_copy_uncomp_pdbID;
   }
   ### PDB file: UPLOAD PDB file
   elsif ($upload_PDB_file ne "") {
      &upload_file ($pdbUpload, $upload_PDB_file);
   }
   ### Verifying that user supplied a file:
   if ($upload_unaligned_file_dna eq "" and $upload_MSA_file_dna eq "") {      
      my $err1="No file provided by the user. Aborting, this is supposed file name:"." $upload_unaligned_file_dna";
      &print_to_output_and_exit($err1,$err1);
   }
   ### to extract SEQRES and ATOMS fasta format from pdb aa
   if ($upload_PDB_file ne "" or $FORM{pdb_ID} ne "") {  #if user ran selecton with a PDB struct.
      &extract_PDB_info();
   }
   ### update page-a.htm file with title and ...
   #&update_page_a_html;  # REMARK: i think that this is not neceserray if not having PE
   ### UPLOAD MSA file
   if ($upload_unaligned_file_dna ne ""){ # user uploaded an un-aligned file
      print LOG "\nUser uploaded UNALIGNED DNA file. seqname: $FORM{msa_SEQNAME}\n";
      &upload_file ($fileUploadPath_dna_unaligned, $upload_unaligned_file_dna); # the first is the UNIX path, the second is the name user gave
   }
   else{ # user uploaded an aligned file
      print LOG "\nUser uploaded ALIGNED DNA file. seqname: $FORM{msa_SEQNAME}\n";
      &upload_file ($fileUploadPath_dna_unaligned, $upload_MSA_file_dna); # the first is the UNIX path, the second is the name user gave
   }
   ### UPLOAD TREE file
   if ($upload_TREE_file ne "") {
      &upload_file ($treeUpload, $upload_TREE_file);  
      &convertNewline($treeUpload) ;
      &removeBPvalues($treeUpload);
      &check_validity_tree_file($treeUpload);         
   }
   # converts Mac or DOS newline to Unix newline
   &convertNewline($fileUploadPath_dna_unaligned);
      
   # create a copy to the original DNA file, since the names of the sequences will be replaced with numbers.
   system("mv $fileUploadPath_dna_unaligned $copied_dna_unaligned");
   
   # if the user supplies an unaligned file - aligns according to codons, and then translate.
   # if the user supplies a codon-aligned file - verifies that there are no internal stop codons and that ORF/3=whole and translates
   &codonAlign(\@sequences_names);
   
   # in case there is a tree file: changes its names accordingly, saving the original file under the name $copied_treeUpload
   if ($upload_TREE_file ne "") {
      system("mv $treeUpload $copied_treeUpload");
      &change_tree_file_names(\@sequences_names, "change_to_numbers", $copied_treeUpload, $treeUpload);
   }
   ###**** SET CORRECT FORMAT MSA file - extract query seq and verify no. of seqs >=5
   my $number_of_sequences_in_msa = &msa_format_extract(\@sequences_names);  

   
   # create a file that will hold all the relevant data for this run.
   &print_data_files;

   system 'echo "(touch '.$WorkingDir.'std.out; chmod 755 '.$WorkingDir.'std.out)" | /bin/tcsh';

   unless (open QSUB_SH, ">$qsub_script")
      {&sys_error_exit("cannot open the file $qsub_script for writing $!\n");}
   print QSUB_SH '#!/bin/sh';
   print QSUB_SH "\nperl $runClac_inQ $WorkingDir $runCalcInput $FormInput > $WorkingDir"."std.out";
   close QSUB_SH;

   chmod 0755, $qsub_script;
   print LOG "submitting qsub job \"SELECTON_$run_name\"\n";

   system "touch $QsubLogFile";
   chmod 0744, $QsubLogFile;
   #chmod 077, $QsubLogFile;
    
   #submitting the qsub_script using a qsub command.
   my $qsub_command = "ssh bioseq\@biocluster qsub -q bioseq -e $WorkingDir -o $WorkingDir -N SELECTON_$run_name $qsub_script";
   if ($estimated_run_time ne "none"){
      if ($model eq "M8" and $number_of_sequences_in_msa>48){ # the estimation for this case is not very accurate, so it is better
         # to give it the maximum run time
         $qsub_command.=" -l walltime=".GENERAL_CONSTANTS::MAX_WALLTIME;
      }
      else{
         $qsub_command.=" -l walltime=$estimated_run_time";
      }
      
   }
   else{
      $qsub_command.=" -l walltime=".GENERAL_CONSTANTS::MAX_WALLTIME;
   }
   print LOG "going to run qsub command: ".$qsub_command."\n";

   my $qsub_job_no = `$qsub_command`;
   chomp($qsub_job_no); # the command returns an extra new line character

   # writing the job number to a file.
   open JOB_NO, ">".$WorkingDir.$Qstat_No_file;
   print JOB_NO $qsub_job_no;
   close JOB_NO;
   chmod 0644, $WorkingDir.$Qstat_No_file;

   $qsub_job_no =~ /(\d+)/;
   $qsub_job_no = $1; # extracting only the process (Q) number
   print LOG "\nafter submitting qsub job. JOB NUMBER: $qsub_job_no\n";
   
   my $ans = BIOSEQUENCE_FUNCTIONS::enqueue_job($qsub_job_no, "Selecton", $run_name);
   print LOG $ans if ($ans ne "ok");


   #unless (open LIST, ">>".GENERAL_CONSTANTS::QUEUING_JOBS){
   #   print LOG "Could not open file ".GENERAL_CONSTANTS::QUEUING_JOBS.". Reason: $!\nThe job was not listed in the queuing_jobs list.\n";
   #   &printTime();
   #}
   #else{
   #   flock LIST, 2; # locks the list, so no other process will write to it. On the same time - if the list is currently locked by another process - it waits until the list file is realeased. The "2" and "8" are the operation symbols for "lock" and "unlock".
   #   print LIST "$qsub_job_no Selecton $run_name ".&printTime()."\n";
   #   flock LIST, 8;
   #   close LIST;
   #}

   print LOG "\n\nCGI ended successfully.";
   close LOG;   
}
else{
   die "Can not fork the process, please contact ".GENERAL_CONSTANTS::ADMIN_EMAIL."\n";
}
exit;

###################################################################################
#                               SUB ROUTINES                                      #
###################################################################################
# Start writing the output web page of Selecton
sub start_output_html {

   system 'echo "(touch '.$OutHtmlFile.'; chmod 755 '.$OutHtmlFile.')" | /bin/tcsh';
   unless (open OUTPUT, ">$OutHtmlFile"){
      print LOG "\nstart_output_html: Cannot open the output file $OutHtmlFile\n";
      exit;
   }
   print LOG "\nstart_output_html: Opening the file $OutHtmlFile, and change the permissions of the WorkingDir\n";
   print OUTPUT <<EndOfHTML;    

<HTML>
<HEAD> <META HTTP-EQUIV="REFRESH" CONTENT=$reload_interval> </HEAD>
<HEAD> <META HTTP-EQUIV="PRAGMA" CONTENT="NO-CACHE"> </HEAD>


<TITLE>Selecton Results $run_name</TITLE>

<style type="text/css">
#menu {

text-decoration: none;
	color: white;
font-size: 12px;
font-weight: 700;
}

</style>




</HEAD>
<BODY bgcolor="#FFF5EE">

<table  width=100%  bgcolor="#400080"> <tr><td>

<table  border=0 cols=6 width=450 bgcolor="#400080" cellpadding=1 cellspacing=0 >
<tr>
<td align=center><a href="/index.html" id=menu target=_top>HOME</a></td>
<td align=center><a href="/overview.html"  id=menu target=_top>OVERVIEW</a></td>
<td align=center><a href="/gallery.html"  id=menu target=_top>GALLERY</a></td>
<td align=center><a href="/faq.html"  id=menu target=_top>FAQ</a></td>
<td align=center><a href="/credits.html"  id=menu target=_top>CREDITS</a></td>
</tr>
</table>

</td><tr></table>


<H1 align=center>Selecton Job Status Page</h1>

<blockquote>

<p><font face=Verdana size=2>
<br><b>Selecton is now processing your request.<br><font size=+1>
<!--job_stat--><a href ="$status_faq" target=stat_faq>Your job status is:</a> Queued<br>
<!--job_pass-->The time that passed since submitting the query is: 00:00<br>
<!--job_time--Estimated run time is: -->
</font>
Please note this may be a lengthy process and an email will be sent to the address you supplied once the calculation is finished.</b><br>
This page will be automatically updated every 30 seconds. You can also reload it manually.<br> 
Once the job has finished, several links to the output files will appear below. 
<br><br>
If you wish to view these results at a later time without recalculating
them, please bookmark this page. The results will be kept on the server for three months.

</font></p>

<h4><font face=Verdana><u>Running Parameters:</u></h4>

EndOfHTML

   print OUTPUT "<p><font face=Verdana size=2>\n";
   print OUTPUT "PDB ID = $FORM{pdb_ID} <br>\n" if ($FORM{pdb_ID} ne "");
   if ($upload_PDB_file ne ""){
      if ($upload_PDB_file =~ m/^.*(\\|\/)(.*)/) {print OUTPUT "PDB file = $2 <br>\n";}
      else {print OUTPUT "PDB file = $upload_PDB_file <br>\n";}
   }
   print OUTPUT "Chain identifier = $FORM{chain} <br>\n" if ($FORM{chain} ne "");
   if ($upload_unaligned_file_dna ne ""){
      if ($upload_unaligned_file_dna =~ m/^.*(\\|\/)(.*)/) {print OUTPUT "DNA unaligned file = $2 <br>\n" ;}
      else {print OUTPUT "DNA unaligned file = $upload_unaligned_file_dna <br>\n" ;}
   }
   if ($upload_MSA_file_dna ne ""){
      if ($upload_MSA_file_dna =~ m/^.*(\\|\/)(.*)/) {print OUTPUT "DNA MSA file = $2 <br>\n";}
      else {print OUTPUT "DNA MSA file = $upload_MSA_file_dna <br>\n" ;}
   }
   print OUTPUT "Query sequence name in MSA file = $FORM{msa_SEQNAME} <br>\n" if ($upload_unaligned_file_dna ne "");
   print OUTPUT "Model = Positive selection enabled (M8, beta + w >= 1) <br>\n" if ($model eq "M8");
   print OUTPUT "Model = Null model: no positive selection(M8a, beta + w = 1) <br>\n" if ($model eq "M8a");
   print OUTPUT "Model = Null model: no positive selection(M7, beta)<br>\n" if ($model eq "M7");
   print OUTPUT "Model = Positive selection enabled(M5, gamma) <br>\n" if ($model eq "M5");
   print OUTPUT "Model = Mechanistic Empirical Combination Model (MEC) <br>\n" if ($model eq "MEC");
   if ($model eq "MEC"){    
      print OUTPUT "Amino-Acid empirical matrix to be expanded =  ";
      print OUTPUT $empiricalMatrix;
      print OUTPUT "<br>\n";    	
   }
   print OUTPUT "Number of categories = $FORM{CATEGORIES} <br>\n" if ($FORM{CATEGORIES} ne "");
   if ($upload_TREE_file ne ""){
      if ($upload_TREE_file =~ m/^.*(\\|\/)(.*)/) {print OUTPUT "User tree file: $2 <br>\n" ;}
      else {print OUTPUT "User tree file: $upload_TREE_file <br>\n" ;}      
   }
   print OUTPUT "</p></font><br>\n\n";
   $cgi_pid = $$;
   chomp($cgi_pid);

   print OUTPUT "\n<FORM ENCTYPE=\"multipart/form-data\" ACTION=\"$kill_job_script\" METHOD=\"POST\">\n";
   print OUTPUT "<INPUT TYPE=\"submit\" VALUE=\"Cancel Selecton Job\"><br>\n";
   print OUTPUT "<INPUT TYPE=hidden NAME=\"Qstat_file\" VALUE=\"".$Qstat_No_file."\">\n";
   print OUTPUT "<INPUT TYPE=hidden NAME=\"selecton_http\" VALUE=\"".$job_canceled_page."\">\n";
   print OUTPUT "<INPUT TYPE=hidden NAME=\"run_no\" VALUE=\"".$WorkingDir."\">\n";
   print OUTPUT "<INPUT TYPE=hidden NAME=\"cgi_pid\" VALUE=\"".$cgi_pid."\">\n";
   print OUTPUT "<h4><u>Messages regarding input and calculation:</u></h4>\n";
   close OUTPUT;   
}
###################################################################################
# checking if PDB FILE exists on aten pdb DB
sub check_copy_uncomp_pdbID {

   #saving the ls output on a file to check file existence
   system 'echo "(ls '.$PdbFilePath.'  > '.$InpSeqFile.';)" | /bin/tcsh';  
   my $file_exist = `cat "$InpSeqFile"`;
   #checking if file pdb exists 
   if ($file_exist eq "") {
      print LOG "\ncheck_copy_uncomp_pdbID:ls $PdbFilePath  > $InpSeqFile";
      &print_to_output_and_exit("The PDB file with the ID \'$FORM{pdb_ID}\' is not found on our database, or not available right now.","\ncheck_copy_uncomp_pdbID: The pdb file $PdbFilePath does not exist on our DB\n");
   }

   #####if pdb exists continues calc
   print LOG "\ncheck_copy_uncomp_pdbID: pdb file exists !\n";
   #***************************************************
   #####copy pdb.Z file locally to work on it freely an uncompress
   
   #copy file to working dir
   print LOG "\ncheck_copy_uncomp_pdbID: Copy $PdbFilePath to $WorkingDir\n";
   system 'echo "(cp '.$PdbFilePath.' '.$WorkingDir.' ;)" | /bin/tcsh';
#    system 'echo "(cd '.$WorkingDir.' ; chmod -R ogu+rx * )" | /bin/tcsh';

   # check if the file exists in the WorkingDir
   my $pdb_file = $WorkingDir . $PdbFileName;
   unless (-e $pdb_file){       
      &sys_error_exit("check_copy_uncomp_pdbID: The PDB file was not copied to the directory $WorkingDir");
   }
     
   # change the permissions of the file and uncompress it
   print LOG "\ncheck_copy_uncomp_pdbID: Change the permissions of $PdbFileName and gunzip it\n";    
   system 'echo "(cd '.$WorkingDir.'; chmod +wxr  '.$PdbFileName.'; gunzip '.$PdbFileName.'; )" | /bin/tcsh';
}
####################################################################################
# UPLOAD file
sub upload_file {
   my $full_path = shift; #NAME to be saved    
   my $file_full_name = shift; #FULL name of file   
   # strip the remote path and keep the filename
   $file_full_name =~ m/^.*(\\|\/)(.*)/; 
   my $name = $2;
   print LOG "\n upload_file : FILE_NAME = $full_path \n";  
   system 'echo "(touch '.$full_path.'; chmod oug+w '.$full_path.')" | /bin/tcsh';
   # if the upload didn't work
   unless (open(UPLOADFILE, ">$full_path")){ 
      &sys_error_exit("upload_file: Can\'t open the file $full_path");
   }
   print LOG "\nupload_file: Upload the file $file_full_name and save it as $full_path\n";   
   while (<$file_full_name>)  {    
      print UPLOADFILE;  
   }
   close UPLOADFILE;

   # verify that the size of the file is not zero
   if (-z $full_path){
      my $err = "upload_file: Cannot upload the file \'$file_full_name\'";
      &send_mailSelecton($err);
      &print_to_output_and_exit("<b>Cannot upload the file \'$file_full_name\', Please verify that the file exists and contains data.</b>", $err);
   }
   system "cd $WorkingDir; chmod ogu+rx *";
   # to be able to write in the PDB file 
   system "chmod ogu+wxr $full_path";
   print LOG "\nupload_file: system 'chmod ogu+wxr $full_path'\n";
   
    #check file type
    my @type = BIOSEQUENCE_FUNCTIONS::check_file_type($full_path);
    if ($type[0] ne "OK"){
        &sys_error_exit("upload_file: ".$type[1]);
    }
    unless ($type[1] eq "PLAIN_TEXT"){
        &print_to_output_and_exit("The file which you have uploaded, '$file_full_name', could not be read by the server. It seems that its type is: $type[1] and not plain text, as required. Make sure your file is in the <a href=\"$http_path/faq.html#q4\">correct format</a>, and refer to the <a href=\"$http_path/faq.html\">FAQ</a> for further assistance.","upload_file : file type is: $type[1]\n");
    }       
}
###########################################################################################
# extract the SEQRES and ATOM sequences from the PDB file 
sub extract_PDB_info {
   my $html_error = "";
   my $log_error = "";
   
   $cmd = "$extractPDBinfo $PdbFileNameUnc $FORM{chain} $FORM{pdb_ID} none $run_name $OutHtmlFile $pdb_data none none $pdb_to_fasta_error none none none none $pdb_fasta $title_file $QuickHelp selecton";

   ### run the script 
   print LOG "\nextract_PDB_info: going to run $cmd\n";    
   #system 'echo "(cd '.$WorkingDir.'; '.$cmd.')" | /bin/tcsh';
   `cd $WorkingDir;$cmd`;

   ## check if the script has created the file error
   if (-e $WorkingDir.$pdb_to_fasta_error and !(-z $WorkingDir.$pdb_to_fasta_error)){
      print LOG "extract_PDB_info: Found error while running $extractPDBinfo. Read the error from the file: $pdb_to_fasta_error\n";
      unless (open ERROR, $WorkingDir.$pdb_to_fasta_error){
         &print_to_output_and_exit($SystemError, "extract_PDB_info: unfortunately could not open the error file. abort.\n");
         &send_mail();
         &stop_reload;
         exit;
      }
      #read the error file
      while (<ERROR>){
         if (/HTML: (.+)/){
            $html_error = $1;
         }
         elsif(/LOG: (.+)/){
            $log_error = $1;
         }
      }
      close ERROR;
      if ($html_error eq "sys" or $html_error eq "") {$html_error = $SystemError;}
      &print_to_output_and_exit($html_error, "extract_PDB_info: $log_error\n");
   }
   # no error was found, make sure that all the files were created
   if (!(-e $WorkingDir.$pdb_data) or !(-e $WorkingDir.$pdb_fasta) or !(-e $WorkingDir.$title_file)){
      &print_to_output_and_exit($SystemError, "extract_PDB_info: The script $extractPDBinfo did not create one of its outputs. \n");
   }

   print LOG "extract_PDB_info: Extracting SEQRES and pdb sequences from $PdbFileNameUnc and chain $FORM{chain} !\n";
}
######################################################################
sub convertNewline{
   my ($dnaFileUnixPath)=@_;
   my $flip_comm="cd $WorkingDir;dos2unix -q $dnaFileUnixPath";
   print LOG "\nconvertNewLine: running dos2unix -q $dnaFileUnixPath\n";
   system "$flip_comm";   
}
###########################################################################################
# remove bootstrap values
sub removeBPvalues {
   my $treeFile=shift;
   my $oldTreeFile = $WorkingDir .'oldUserTreeFile.txt'; 
   TREE_parser::removeBPvalues($treeFile, $oldTreeFile);
}
###########################################################################################
# check the validity of the newick format of the uploaded tree 
sub check_validity_tree_file {
   my $treeFile=shift;
   my $lineCounter=0;
   my $rightBrackets=0;
   my $leftBrackets=0;
   my @lineArr;
   my $line;
   my $errorBool = 0;
   my $noRegularFormatChar;
   my $treeFileOneLine;
   my $read_right_bracket = "no";
   my $tempTreeFile = $WorkingDir .'TempTreeFile.txt'; 
   open(TREEFILE,"$treeFile");
   while (<TREEFILE>) {
      $line = $_;
      chomp($line);
      $treeFileOneLine .= $line;
      $lineCounter++;
   }
   close TREEFILE;   
   $line =  $treeFileOneLine;
   open OUTPUT, ">>$OutHtmlFile";	       
   if ( $lineCounter>1) {
      open TEMPTREEFILE, ">$tempTreeFile";
      print TEMPTREEFILE  $line; 
      unlink ($treeFile);
      system 'echo "(cp -r '.$tempTreeFile.' '.$treeFile.'; chmod -R ogu+xr '.$treeFile.')" | /bin/tcsh';
      unlink ($tempTreeFile); 
   }   
   @lineArr=split(//,$line); 
   foreach my $chain(@lineArr) {
      if ($chain eq '(')  {
         $leftBrackets++;
         $read_right_bracket = "no";
      }
      elsif ($chain eq ')') {
         $rightBrackets++;
         $read_right_bracket = "yes";
      }	
      elsif ($chain =~ /([\!|\@|\#|\$|\^|\&|\*|\~|\`|\{|\}|\'|\?|\\|\/|\<|\>])/){
         $noRegularFormatChar .= " \"$1\", " if $noRegularFormatChar !~ /\Q$1\E/;
         $read_right_bracket = "no";
      }
      # if right after a right Bracket we read a character which is not legal (ie: , : ;) we output a message to the user, since we don't handle bootstrap values or internal node names
      else{
         if($read_right_bracket eq "yes"){
            if($chain =~ /\d/){
               &print_to_output_and_exit("<b>The <A HREF=$userTree TARGET=MSA_window>TREE file</A> you uploaded includes bootstrap values. Please remove them and resubmit your query.</b>\n", "check_validity_tree_file : found bootstrap values. Abort\n");
            }
            elsif($chain !~ /[,|:|;]/){
               &print_to_output_and_exit("<b>The <A HREF=$userTree TARGET=MSA_window>TREE file</A> you uploaded includes internal nodes names. Please remove them and resubmit your query.</b>\n", "check_validity_tree_file : found internal nodes name in the tree. Abort\n");
            }
         }
        $read_right_bracket = "no";
      }
   }
   if ($leftBrackets ne $rightBrackets) {
      print OUTPUT "\n<p>$ErrorDef<br><b>The <A HREF=$userTree TARGET=MSA_window>TREE file</A> which appears to be in Phylip format is missing parentheses</p>\n";
      print OUTPUT $ContactDef;
      $errorBool++;            
   } 
   if ($noRegularFormatChar =~ /.+/)  {
      $noRegularFormatChar =~ s/\,\s$//;
      print OUTPUT "\n<p>$ErrorDef<br><b>The <A HREF=$userTree TARGET=MSA_window>TREE file</A> which appears to be in Phylip format, contains the following non-standard characters: ". qq($noRegularFormatChar) . ".</p>\n";
      print OUTPUT $ContactDef;
      $errorBool++;            
   }
   close OUTPUT;
   if ($errorBool ne '0') {
      &send_mail();
      &stop_reload;
      exit;
   } 
}
####################################################################################
sub codonAlign{
   my @ans;
   print LOG "\ncodonAlign ";
   
   my $ref_seq_name = shift;
   #my %original_seq_name; # we use this hash to hold the original names for
   
   my %codonTable = (0=>'1', 1=>'15', 2=>'6', 3=>'10', 4=>'2',             # the codon Table input represents a number for each table, as was decided in the html.
                       5=>'5', 6=>'3', 7=>'13', 8=>'9', 9=>'14', 10=>'4',);  # since the bioPerl modoule uses different numbers, we match each "input" number to the bioPerl table number.
       
      # user supplied a codon aligned file -
   if ($upload_unaligned_file_dna eq ""){
      print LOG "calling codonAlign::DNA_checkLegal_and_crate_AAFile($copied_dna_unaligned, $fileDna_aligned, $WorkingDir, $codonTable{$FORM{GENCODE}}, $ref_seq_name, $fileName_amino_aligned, $OutHtmlFile, $WWWdir)\n";
      @ans = codonAlign::DNA_checkLegal_and_crate_AAFile($copied_dna_unaligned, $fileDna_aligned, $WorkingDir, $codonTable{$FORM{GENCODE}}, $ref_seq_name, $fileName_amino_aligned, $OutHtmlFile, $WWWdir);
   }
   else{
      print LOG "calling codonAlign::DNA_align($copied_dna_unaligned, $fileName_amino_aligned, $WorkingDir, $fileDna_aligned, $OutHtmlFile, $muscle, $WWWdir, $codonTable{$FORM{GENCODE}}, $ref_seq_name)\n";
      @ans = codonAlign::DNA_align($copied_dna_unaligned, $fileName_amino_aligned, $WorkingDir, $fileDna_aligned, $OutHtmlFile, $muscle, $WWWdir, $codonTable{$FORM{GENCODE}}, $ref_seq_name);
   }
   unless($ans[0] eq "ok"){
      my $err = $ans[1];
      # in case it is a user error, we let him know via the html file
      if ($ans[0] eq "user"){
         &print_to_output_and_exit("An Error was found in your input file: $err", "codonAlign : error in codonAlign.pm : $err");
      }
      # a system error
      else{
         &sys_error_exit($err);
      }
   }
   print LOG "codonAlign.pm returned OK\n";      
}
######################################################################################
sub change_tree_file_names{
   
   my $ref_sequences_names = shift;
   my $new_tree_mode = shift;    # change_to_numbers OR change_to_sequences
   my $input_tree = shift;
   my $output_tree = shift;
   my ($tree, $err);
   
   unless (open TREE, $input_tree){
      if ($new_tree_mode eq "change_to_numbers") {
         &sys_error_exit("change_tree_file_names : can not open file $input_tree for reading.") ;}
      else{
         print LOG "could not open the file $input_tree for reading. the tree file will be presented with numbers\n";
         return;}
   }
   #check validity of input tree file
   if ($new_tree_mode eq "change_to_numbers"){
      print LOG "\nchange_tree_file_names : reading tree file\n";
      while (<TREE>){
         if ((/.+\r.+;/)||(/.+\n.+;/)){
            $err = "Return or newLine charachters were found inside the tree file.";
            &print_to_output_and_exit("An Error was found in your input tree file: $err<br>Please note: The input tree file should be written in one line. See <a href=\"$tree_faq\">Selecton accepted format</a> for more info.<br>\nPlease correct your input tree file and re-submit your query.<br>\n", "change_tree_file_names :: $err");
         }
         
         # chooping last ^M or \n
         elsif((/^\(.+:.+\).*;\r$/)||(/^\(.+:.+\).*;\n$/)){
            chop;
         }
         # if there is more than one tree - extracting the first tree
         if(/(.+:.+;)(.+)/){
            print OUTPUT "\n<p><ul><li><font color='red'><b>Warning:</b></font> There is more than one tree in your <a href=$userTree_copy>input tree file</a>. The first tree will be used for calculations.</li></ul></p>\n";
            $tree = $1;
            while ($tree =~ m/(.+\;)(.+)/){
               $tree = $1;
            }        
         }
         # the minimum requirments from a tree:
         #elsif(/^\(.+:.+\).*;$/){
         #  $tree=$_;
         #}
         elsif(/^\(.+\).*;$/){
            $tree=$_;
         }
         else{
            $err = "The input tree file is not in a legal format.\n";
            &print_to_output_and_exit("An Error was found in your input tree file:<br>$err See <a href=\"$tree_faq\">Selecton accepted format</a> for more info.<br>\nPlease correct your input tree file and re-submit your query.<br>", "change_tree_file_names :: $err");
         }         
      }
      print LOG "change_tree_file_names : tree input is legal\n";
   }
   else{
      $tree = <TREE>;
   }
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
         if ($new_tree_mode eq "change_to_numbers"){
            for (my $in=1; $in<=($#$ref_sequences_names+1); $in++){               
               $seq_found = "no";
               if ($ref_sequences_names->[$in] eq $2){
                  $seq_found = "yes";
                  $final_tree.= $1.$in.$rest_of_exp;
                  last;
               }               
            }                        
            if ($seq_found eq "no") {
            # in case a sequence was found in the tree and not in the DNA file            
               &print_to_output_and_exit("The sequence name: \"$2\" was found in your tree file, but was not found in your DNA input file.<br>When submitting an input tree file, sequence names of both inputs must be identical. See <a href=\"$tree_faq\">Selecton accepted format</a> for more info.<br>\nPlease correct your input files and re-submit your query.<br>\n","change_tree_file_names : the sequence name $2 appears in tree file, does not appear in DNA input file");
            }
         }
         else{
            $final_tree.= $1.$ref_sequences_names->[$2].$rest_of_exp;
         }
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
######################################################################################
# extract the query sequence from the user-provided MSA file
sub msa_format_extract {
 	
   my $ref_seq_name = shift; 
   
   print LOG "\nentered msa_format_extract : \n";
   
   my $original_seq_name = $FORM{msa_SEQNAME};  # we change the FORM name variable later, in order to search it in the DNA MSA that was built using numbers instead of original names.
   
   if ($upload_PDB_file ne "" or $FORM{pdb_ID} ne "") {  #if user ran selecton with a PDB struct, we change the name of the file ".pdbfasta" which was created by the extract_PDB_info routine and give writing permissions to all

      my $str1 = "mv $pdb_fasta $pdb_msa";
      my $str2 = "chmod ogu+wrx $pdb_msa";
      print LOG "\nmsa_format_extract: mv $pdb_fasta $pdb_msa\n chmod ogu+wrx $pdb_msa\n";
      system 'echo "(cd '.$WorkingDir.'; '.$str1.'; '.$str2.')" | /bin/tcsh';      
   }    
   # Change the sequence name, so it will fit the number that represents the sequence in the DNA file    
   for (my $i=1; $i<=($#$ref_seq_name+1); $i++){
      if ($FORM{msa_SEQNAME} eq $ref_seq_name->[$i]){
         $FORM{msa_SEQNAME} = $i;
         last;
      }
   }

   print LOG "\nmsa_format_extract: SeqName = \'$FORM{msa_SEQNAME}\'\n";

   ### verify that there at least 3 sequences, else - kill the script #changed from 5 to 3 on 1/2/06
   my $msa = $WorkingDir.$fileName_amino_aligned;
   unless (open FILE, "<$msa"){
      &sys_error_exit("msa_format_extract: can not open file $msa for reading\n");
   }
    
   my $counter = 0;
   my $TargetFound = 0;
   my $querySeq="";
   my $space = 0;
   my $firstSequence;
   my $firstSequenceName;
   my $seq_length;
   my $first_seq_length;
    
   # replace '()' with '_'  - that's what clustalw does! 
   my $inFile  = Bio::SeqIO->new('-file' => "$msa" , '-format' => 'Fasta');
   while ( my $seqObj = $inFile->next_seq() ) {
      #obtaining first sequence - if query not found, 1st seq. used as query
      if ($counter == 0) {
         $firstSequence = $seqObj->seq();
         $first_seq_length = length($firstSequence);
         $firstSequenceName = $seqObj->display_id();
         if ($seqObj->desc() ne ""){
            $firstSequenceName .= " ".$seqObj->desc();
         }
      }      
      $counter++;
      my $seq = $seqObj->seq();
      my $name = $seqObj->display_id();
      $seq_length = length($seq);
      if ($seqObj->desc() ne ""){
         $name .= " ".$seqObj->desc();
      }
      #if ($name =~ /^(\Q$FORM{msa_SEQNAME}\E\s*)/){ #find query sequence name
      if ($name =~ /^($FORM{msa_SEQNAME}\s*)$/){ #find query sequence name
         $querySeqFoundinMSA = "yes";
         $query= $1;
         $TargetFound = 1;
         $querySeq = $seq;
         ($querySeq) =~ s/\W+//g;
      }
   }
   $inFile->close();
   unless ($querySeqFoundinMSA eq "yes") { # if query sequence not found - use first sequence
      $querySeq=$firstSequence;
      $query=$firstSequenceName;
   }	
   $querySeq =~ s/-//g; #clear all gaps in the sequence
   close FILE;
   # The query sequnece name is not found in the MSA - warning
   if ($TargetFound == 0){
      open OUTPUT, ">>$OutHtmlFile";
      print OUTPUT "\n<p><ul><li><font color='red'><b>Warning:</b></font> The query sequence name \'$original_seq_name\' is not found in the <A HREF=COPY_$fileUploadName_dna_unaligned TARGET=MSA_window>input file</A>.<br>The calculation continues. The first sequence in MSA is used as query.</li></ul></p>\n";
      close OUTPUT;
      print LOG "\nmsa_format_extract: The query sequence name is not found in the MSA file. Calculation continues with 1st seq. in file\n";
   }
   if ($upload_PDB_file ne "" or $FORM{pdb_ID} ne "") {  #if user ran selecton with a PDB struct.
   ### add the target sequence to the file _PDB_MSA.pdbfasta
      unless (open PDBFASTA, ">>".$WorkingDir.$pdb_msa) {
         &sys_error_exit("msa_format_extract: Can\'t open the file $pdb_msa");
      }
      print PDBFASTA "\n>$query\n$querySeq\n";
      close PDBFASTA;
      
      ### Call 'alignPDB2refseq' to run clustalw and check the results
      &alignPDB2refseq($pdb_msa, $clustal_outFile);
   }  
    # there are less than 3 sequences - write error message and exit
   if ($counter < 3){ #changed from 5 to 3 on 1/2/06
      open OUTPUT, ">>$OutHtmlFile";
      if ($counter == 1){
         &print_to_output_and_exit("The <A HREF=COPY_$fileUploadName_dna_unaligned TARGET=MSA_window>input file</A> contains only one sequence. The minimal number of homologues required for the calculation is 5.<br>(Make sure that the file is saved as plain text and does not contain special characters).", "msa_format_extract: The MSA file contains only $counter sequences");
      }
      else {
         &print_to_output_and_exit("The <A HREF=COPY_$fileUploadName_dna_unaligned TARGET=MSA_window>input file</A> contains only $counter sequences. The minimal number of homologues required for the calculation is 5.", "msa_format_extract: The MSA file contains only $counter sequences");
      }
   }

   # there are less than 10 sequences - write a warning
   elsif ($counter < 10){	
      open OUTPUT, ">>$OutHtmlFile";
      print OUTPUT "\n<p><ul><li><font color='red'><b>Warning:</font></b> The MSA file contains only $counter sequences. It is recommended to use an MSA file with at least 10 homologues. The calculation continues nevertheless.</li></ul></p>\n";
      close OUTPUT;
   }	
 # more than 100 seqeunces - not supported by the server
   elsif ($counter > 100) {
      &print_to_output_and_exit("The <A HREF=COPY_$fileUploadName_dna_unaligned TARGET=MSA_window>MSA file</A> contains $counter sequences. The Selecton server supports only < 100 sequences. For longer runs, please download the source code under SOURCE and install the program locally .", "msa_format_extract: The MSA file contains $counter (>100) sequences, Selecton aborted.");
   }
   else {	
      open OUTPUT, ">>$OutHtmlFile";
      print OUTPUT "\n<p><ul><li>The calculation is performed on the $counter sequences obtained from the MSA file.</li></ul></p>\n";
      close OUTPUT;
   }
   
   #calculate the estimated 
   if ($model eq "M8" or $model eq "MEC"){
      $estimated_run_time = BIOSEQUENCE_FUNCTIONS::selecton_estimated_run_time( $first_seq_length, $counter, $model);
      print LOG "\nmsa_format_extract : sending selecton_estimated_run_time( $first_seq_length, $counter, $model), time is: $estimated_run_time\n";
   }
   open STATISTICS, ">>".$statistics_file;
   flock STATISTICS, 2;
   print STATISTICS "$run_name $model $counter $first_seq_length\n";
   flock STATISTICS, 8;
   close STATISTICS;
   
   return $counter;
}
###########################################################################
# The function runs 'clustalw' with two sequences,
# and checks their pairwise alignment.
# The arguments: 
# 1. the parameters to run clustalw (input-file and additional parameters)
# 2. an output file
###########################################################################
sub alignPDB2refseq{

   my $clustalw_parms = shift;
   my $outfile = shift;

   my %message = (SEQRES => "sequence extracted from the SEQRES field of the PDB file",
		   PDB => "sequence extracted from the ATOM field of the PDB file ",
		   MSA => "query sequence extracted from the MSA file"); 
   my %SeqName = (SEQRES => "SEQRES sequence",
		   PDB => "ATOM sequence",
		   MSA => "query sequence");

   # run clustalw with the given parameters
   my $command = $clustalw . " ". $WorkingDir.$clustalw_parms . " > " . $WorkingDir.$outfile;
   print LOG "\nalignPDB2refseq: running $command\n";
   system 'echo "(cd '.$WorkingDir.'; '.$command.')" | /bin/tcsh';

   my $full_outfile = $WorkingDir . $outfile;
   unless (open OUT, "<$full_outfile"){
      &sys_error_exit("alignPDB2refseq: Can\'t open the file $full_outfile");
   }

   my $msa_seq;
   my $pdb_length ;
   my $msa_seq_length;
   my $score;
   # search the outfile for the name of the sequences, their length
   # and the score of the alignment
   while (<OUT>){
      if ($_ =~ /Sequence\s+\d:\s+PDB_\S?\s+(\d+)\s+aa/i){           
         $pdb_length = $1;
      }
      elsif ($_ =~ /Sequence\s+\d:\s+(.+)\s+(\d+)\s+aa/i){
         $msa_seq = "MSA";
         $msa_seq_length = $2;
      }
      elsif ($_ =~ /Sequences.+Aligned.+Score:\s+(\d+)/i){           
         $score = $1;
      }
   }   
   close OUT;

   my $pdb = $PdbPrefix . ".ent";
   system 'echo "(cd '.$WorkingDir.' ; chmod -R og+rx * )" | /bin/tcsh';    
    ### check if the atoms sequence is longer than the SEQRES sequence
   if ($msa_seq_length < $pdb_length){
      # significant difference - stop the script!
      if ($msa_seq_length < ($pdb_length * 0.9)){	 
         open OUTPUT, ">>$OutHtmlFile";
         print OUTPUT "\n<p><ul><li><font color='red'><b>Warning:</b></font> The $message{MSA} is shorter than the $message{PDB}. The $SeqName{MSA} has $msa_seq_length residues and the $SeqName{PDB} has $pdb_length residues. The calculation continues nevertheless.</li></ul></p>\n";
         close OUTPUT;	 
      }
      # just write a warning
      else {
	 open OUTPUT, ">>$OutHtmlFile"; 
	 if ($FORM{pdb_ID} eq "FILE"){	    
	    print OUTPUT "\n<p><ul><li><font color='red'><b>Warning:</b></font> The $message{MSA} is shorter than the $message{PDB}. The $SeqName{MSA} has $msa_seq_length residues and the $SeqName{PDB} has $pdb_length residues. The calculation continues nevertheless.</li></ul></p>\n";
	 }
	 else {
	    print OUTPUT "\n<p><ul><li>The $message{MSA} is shorter than the $message{PDB}. The $SeqName{MSA} has $msa_seq_length residues and the $SeqName{PDB} has $pdb_length residues. The calculation continues nevertheless.</li></ul></p>\n";
	 }
	 close OUTPUT;
      }
   }

   ### check if the SEQRES is longer than the atoms sequence
   elsif ($pdb_length < $msa_seq_length){
      # significant difference - stop the script!
      if ($pdb_length < ($msa_seq_length * 0.2)){
         open OUTPUT, ">>$OutHtmlFile"; 
         print OUTPUT "\n<p>$ErrorDef<br><b>The $message{PDB} is significantly shorter than the $message{MSA}. The $SeqName{MSA} has $msa_seq_length residues and the $SeqName{PDB} has $pdb_length residues.</b></p>\n";
         close OUTPUT;
         print LOG "\nalignPDB2refseq: The $message{PDB} is significantly shorter than the $message{MSA}\n";
      }
      # write a warning
      else {
         open OUTPUT, ">>$OutHtmlFile";
         if ($FORM{pdb_ID} eq "FILE"){		
            print OUTPUT "\n<p><ul><li><font color='red'><b>Warning:</b></font> The $message{PDB} is shorter than the $message{MSA}. The $SeqName{MSA} has $msa_seq_length residues and the $SeqName{PDB} has $pdb_length residues. The calculation continues nevertheless.</li></ul></p>\n";
         }   
         else {
            print OUTPUT "\n<p><ul><li>The $message{PDB} is shorter than the $message{MSA}. The $SeqName{MSA} has $msa_seq_length residues and the $SeqName{PDB} has $pdb_length residues. The calculation continues nevertheless.</li></ul></p>\n";
         }
         close OUTPUT;
      }
   }

   ### check if the score is not 100
   if ($score < 100){
      # if the alinment score < 60, write a message and stop the script. 
      #if ($score < 60){
      #   &print_to_output_and_exit("The Score of the alignment between the $message{MSA} and the $message{PDB} is ONLY ID% = $score .<br>(<A HREF=$clustal_aligned_file TARGET=Alignment_window>Pairwise Alignment</A>)","alignPDB2refseq: The Score of the alignment between the $message{MSA} and the $message{PDB} is ONLY ID% = $score.");
      #}
      #else {
	 open OUTPUT, ">>$OutHtmlFile";
	 if ($FORM{pdb_ID} eq "FILE"){		
	    print OUTPUT "\n<p><ul><li><font color='red'><b>Warning:</b></font> The Score of the alignment between the $message{MSA} and the $message{PDB} is ID% = $score . The calculation continues nevertheless. (<A HREF=$clustal_aligned_file TARGET=Alignment_window>Pairwise Alignment</A>)</li></ul></p>\n";
	 }
	 else {
	    print OUTPUT "\n<p><ul><li><font color='red'><b>Warning:</b></font>The Score of the alignment between the $message{MSA} and the $message{PDB} is ID% = $score . The calculation continues nevertheless. (<A HREF=$clustal_aligned_file TARGET=Alignment_window>Pairwise Alignment</A>)</li></ul></p>\n";
	 }
	 close OUTPUT;
      }
   #}
}
######################################################################
sub print_data_files{
   # creating a file to hold the sequences names, as we can't send a pointer to the run_calc script
# ------ storable on -------
   # the key is the same string as the variable name in the script Selecton_run_calc, so that the
   # retriving process will be easy.
   print LOG "\nEntered print_data_files. going to create hash\n";
   my %run_data = ();   
   $run_data{run_name} = $run_name; $run_data{WorkingDir} = $WorkingDir; $run_data{WWWdir} = $WWWdir;
   $run_data{epsilonPrecision} = $epsilonPrecision; $run_data{query_seq_name_to_run} = $query; $run_data{optimizeBL} = $optimizeBL; $run_data{querySeqFoundinMSA} = $querySeqFoundinMSA; $run_data{tree_faq} = $tree_faq; $run_data{method} = $method; $run_data{SysErrorDef} = $SysErrorDef; $run_data{ErrorDef} = $ErrorDef; $run_data{ContactDef} = $ContactDef; $run_data{fileDna_aligned} = $fileDna_aligned; $run_data{treeUpload} = $treeUpload; $run_data{OutHtmlFile} = $OutHtmlFile; $run_data{cgi_log_file} = $OutLogFile; $run_data{OutLogFile} = $QsubLogFile; $run_data{fileName_amino_aligned} = $fileName_amino_aligned; $run_data{OutputURL} = $OutputURL; $run_data{estimated_run_time} = $estimated_run_time; $run_data{sequences_names_file} = $sequences_names_file; $run_data{selecton_log_dir} =$Logs_dir;
   if ($upload_MSA_file_dna eq "") {$run_data{upload_MSA_file_dna} = "no";}
      else {$run_data{upload_MSA_file_dna} = "USER_DNA_ALIGNED_FILE";} 
   if ($upload_unaligned_file_dna eq "") {$run_data{upload_unaligned_file_dna} = "no";}
      else {$run_data{upload_unaligned_file_dna} = "USER_DNA_UNALIGNED_FILE";}
   if ($upload_TREE_file ne "") {$run_data{upload_TREE_file} = "GIVEN";}
      else {$run_data{upload_TREE_file} = "NOT_GIVEN";}
   if ($upload_PDB_file ne "" or $FORM{pdb_ID} ne "") {$run_data{was_Pdb_uploaded} = "yes";}
    else {$run_data{was_Pdb_uploaded} = "no";}
   if ($FORM{pdb_ID} ne "") {
      $run_data{PdbFileNameUnc} = $PdbFileNameUnc;
      $run_data{PdbPrefix} = $PdbPrefix;
      $run_data{pdb_data} = $pdb_data;
      $run_data{clustal_aligned_file} = $clustal_aligned_file
   }
   else{
      $run_data{PdbFileNameUnc}= "NOT_GIVEN";
      $run_data{PdbPrefix} = "NOT_GIVEN";
      $run_data{pdb_data} = "NOT_GIVEN";
      $run_data{clustal_aligned_file} = "NOT_GIVEN";
   }
   print LOG "print_data_files : hash contains ".(keys %run_data)." keys. Going to Store in ".$WorkingDir.$runCalcInput."\n";
   unless (open INPUT, ">".$WorkingDir.$runCalcInput){
      print LOG "cannot open ".$WorkingDir.$runCalcInput." $!\n";
      exit;
   }
   print INPUT "\n";
   close INPUT;
   chmod 0755, $WorkingDir.$runCalcInput;
   $! = "";
   unless(store \%run_data, $WorkingDir.$runCalcInput)
   { print LOG "\ncannot store $!\n";
    exit;
   }
   else{
      print LOG "print_data_files : managed to store! $! \n";
   }
   store \%FORM, $WorkingDir.$FormInput; # when storable off, no need for this file

# ------ storable on -------  
   unless (open SEQ_NAMES, ">".$WorkingDir.$sequences_names_file){
      &sys_error_exit("cannot open the file ".$WorkingDir.$sequences_names_file." for writing $!\n");
   }
   my $index=0;
   foreach(@sequences_names){
      print SEQ_NAMES $index." $_\n";
      $index++;
   }
   close SEQ_NAMES;
   chmod 0755, $WorkingDir.$sequences_names_file;
   chmod 0600, $WorkingDir.$runCalcInput;
# ------ storable off -------

   #unless (open RUN_CALC, ">".$WorkingDir.$runCalcInput){
   #   &sys_error_exit("cannot open the file ".$WorkingDir.$runCalcInput." for writing $!\n");}
   #print RUN_CALC "RUN NAME: $run_name\n";
   #print RUN_CALC "WORKING DIR: $WorkingDir\n";   
   #print RUN_CALC "WWW DIR: $WWWdir\n";
   #print RUN_CALC "PRECISION LEVEL: $epsilonPrecision\n";
   #print RUN_CALC "EVOLUTONARY MODEL: $FORM{MODEL}\n";
   #($FORM{EMPIRICAL_MATRIX} ne "") ? print RUN_CALC "EMPIRICAL MATRIX: $FORM{EMPIRICAL_MATRIX}\n" : print RUN_CALC "EMPIRICAL MATRIX: NOT_GIVEN\n"; #ONLY IF IT IS MEC MODEL
   #print RUN_CALC "QUERY NAME TO RUN: $query\n";
   #print RUN_CALC "DISTRIBUTE CATEGORIES: $FORM{CATEGORIES}\n";
   #($upload_TREE_file ne "") ? print RUN_CALC "TREE_WAS_UPLOADED?: $upload_TREE_file\n" : print RUN_CALC "TREE_WAS_UPLOADED?: NOT_GIVEN\n";  #REMARK: CHANGE THIS VAR'S CONTENT IN THE REST OF THE SCRIPT TO TRUE/FALSE.
   #print RUN_CALC "OPTIMIZE BRANCH LENGTH? $optimizeBL\n";
   #print RUN_CALC "GENETIC CODE: $FORM{GENCODE}\n";
   #print RUN_CALC "GIVEN QUERY NAME: $FORM{msa_SEQNAME}\n";
   #($FORM{pdb_ID} ne "") ? print RUN_CALC "PDB ID: $FORM{pdb_ID}\nPDB NAME: $PdbFileNameUnc\nPDB PREFIX: $PdbPrefix\nPDB DATA FILE: $pdb_data\nCLUSTAL ALN: $clustal_aligned_file\n" : print RUN_CALC "PDB ID: NOT_GIVEN\nPDB NAME: NOT_GIVEN\nPDB PREFIX: NOT_GIVEN\nPDB DATA FILE: NOT_GIVEN\nCLUSTAL ALN: NOT_GIVEN\n";
   #($FORM{chain} ne "") ? print RUN_CALC "PDB CHAIN: $FORM{chain}\n" : print RUN_CALC "PDB CHAIN: NOT_GIVEN\n";
   #print RUN_CALC "FOUND QUERY IN MSA?: $querySeqFoundinMSA\n";
   #print RUN_CALC "TREE FAQ: $tree_faq\n";   
   #print RUN_CALC "METHOD: $method\n";
   #($recipient ne "") ? print RUN_CALC "USER EMAIL: $recipient\n" : print RUN_CALC "USER EMAIL: NOT_GIVEN\n";
   #print RUN_CALC "SYS ERROR: $SysErrorDef\n";
   #print RUN_CALC "ERROR DEF: $ErrorDef\n";
   #print RUN_CALC "CONTACT DEFINITION: $ContactDef\n";
   #print RUN_CALC "WAS PDB UPLOADED?: ";
   #($upload_PDB_file ne "" or $FORM{pdb_ID} ne "")? print RUN_CALC "yes\n" : print RUN_CALC "no\n";
   #print RUN_CALC "DNA FILE NAME: $fileDna_aligned\n";
   #print RUN_CALC "UPLOADED TREE PATH: $treeUpload\n";
   #print RUN_CALC "OUTPUT HTML PATH: $OutHtmlFile\n";
   #print RUN_CALC "LOG PATH: $OutLogFile\n";
   #print RUN_CALC "QSUB LOG: $QsubLogFile\n";
   #print RUN_CALC "LOG DIR: $Logs_dir\n";
   #print RUN_CALC "SEQ NAMES FILE: $sequences_names_file\n";
   #print RUN_CALC "AMINO FILE NAME: $fileName_amino_aligned\n";
   #print RUN_CALC "WAS A DNA UNALIGNED FILE UPLOADED?: ";
   #($upload_unaligned_file_dna eq "") ? print RUN_CALC "no\n" : print RUN_CALC $upload_unaligned_file_dna."\n";
   #print RUN_CALC "WAS A DNA ALIGNED FILE UPLOADED?: ";
   #($upload_MSA_file_dna eq "") ? print RUN_CALC "no\n" : print RUN_CALC $upload_MSA_file_dna."\n";   
   #print RUN_CALC "URL OUTPUT: $OutputURL\n";
   
   #close RUN_CALC;
   #chmod 0755, $WorkingDir.$runCalcInput;
# ------ storable off -------   
}
######################################################################
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
###########################################################################################
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
   sleep ($reload_interval);
   open OUTPUT, "<$OutHtmlFile";
   my @output = <OUTPUT>;
   close OUTPUT;   
   open OUTPUT, ">$OutHtmlFile";
   foreach my $line (@output){  # we remove the refresh lines and the button which codes for Selecton cancelled job
      unless ($line =~ /REFRESH/ or $line =~ /NO-CACHE/ or $line =~ /ACTION=\"$kill_job_script/ or
      $line =~ /VALUE=\"Cancel Selecton Job\"/ or $line =~ /TYPE=hidden NAME=/ or $line =~ /<!--job_/){
         print OUTPUT $line;
      }
   }
   close OUTPUT;
   print LOG "\n\nEnd time: ";
   &printTime();
   close LOG;
   
    # remove the job from the running jobs list
    open LIST, "+>>".GENERAL_CONSTANTS::SELECTON_RUNNING_JOBS;
    flock LIST, 2;
    seek LIST, 0, 0; #rewind the pointer to the beginning
    my @all_lines_in_list = <LIST>; # read the contents into the array
    truncate LIST, 0; # remove all the information, The 0 represents the size of the file that we want
    foreach (@all_lines_in_list){
        chomp;
        unless(/$run_name/){
            print LIST $_."\n";
        }
    }
    flock LIST, 8;
    close LIST;
   chmod 0600, $WorkingDir. "user_email.txt";
   
#    if (-e $WorkingDir."core"){
#	print LOG "remove core file from working directory\n";
#	unlink $WorkingDir."core";
#    }
}

#########################################################################################
# Sends an automatic mail when there are errors
sub send_mail { # to user
   
   $email_subject = "Error in Selecton running";
   $email_message = "Hello!\n\nUnfortunately there was an error while running Selecton.\nPlease click on the following link to see more details\n We apologize for the inconvenience\n\n$OutputURL\n";   
   print LOG "send_mail: sending system error to user\n";
   chdir $send_email_dir;
   $email_system_return = system ('./sendEmail -f '.GENERAL_CONSTANTS::ADMIN_EMAIL.' -t $recipient -u '.$email_subject.' -xu '.$userName.' -xp '.$userPass.' -s '.$smtp_server.' -m '.$email_message);
   unless ($email_system_return =~ /successfully/)    {
      print LOG "send_mail: The message was not sent successfully. system returned: $email_system_return\n";
   }   
}
#########################################################################################
sub send_mailSelecton{ # to selecton administrator
   my $email_message = shift;
   $email_subject = "Error in Selecton running $run_name";
   print LOG "send_mailSelecton: send error message to admin\n";
   chdir $send_email_dir;
   $email_system_return = system ('./sendEmail -f '.GENERAL_CONSTANTS::ADMIN_EMAIL.' -t '.GENERAL_CONSTANTS::ADMIN_EMAIL.' -u '.$email_subject.' -xu '.$userName.' -xp '.$userPass.' -s '.$smtp_server.' -m '.$email_message."\n User's email is: $recipient\n");
   unless ($email_system_return =~ /successfully/)    {
      print LOG "The message was not sent successfully. system returned: $email_system_return\n";
   }
}
#########################################################################################
# this function prints the time to the LOG file.
# if used with return arguments: returns time and date in a different format than printed to log (only numbers).
sub printTime {
   my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
   my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
   my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
   my $year = 1900 + $yearOffset;
   my $theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
   print LOG $theTime;   
   
   $second = &convertNum($second);
   $minute = &convertNum($minute);
   $hour = &convertNum($hour);
   $month = &convertNum($month+1);
   $dayOfMonth = &convertNum($dayOfMonth);
   
   return "$hour:$minute:$second $dayOfMonth-".$month."-$year";

}
#########################################################################################
# converts a number from one digit to 2 digits
sub convertNum 
{
    my $input_num = shift;
    if ($input_num < 10)
        {return "0".$input_num;}
    else
        {return $input_num;}
}