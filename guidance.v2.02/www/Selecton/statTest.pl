#!/usr/bin/perl -w


use strict;
use Storable;
use lib "/bioseq/bioSequence_scripts_and_constants";  
use GENERAL_CONSTANTS;
use BIOSEQUENCE_FUNCTIONS;
use Bio::SeqIO;


# print "Content-type:text/html\n\n";

#my $query = new CGI;
#my $run_number = $query->param("run_number");
#my $FORM{MODEL} = $query->param("model");
#my $FORM{email_address} = $query->param("email");
#my $runCalcInput = $query->param("input_file");


my $run_number = shift;
my $runCalcInput = shift;

my %run_data = ();
my %FORM = ();

my ($WWWdir, $linkToOldResults, $outputScoreFile, $params, $log, $oldCgiLog, $oldHtml, $ErrorDef, $SysErrorDef, $ContactDef, $LogsDir, $formInput, $workingDir,$nullDir, $OutHtmlFile);
my $qsub_ans_file = "qsub_ans.txt"; # The file which is read by the daemon process

&read_input;

$WWWdir.="nullModel/";
my $newCgiLog = $LogsDir . $run_number."null_.log";
my $finish_flag = $nullDir."END_OK";



### Sending e-mail from ibis
my $send_email_dir = GENERAL_CONSTANTS::SEND_EMAIL_DIR;
my $smtp_server = GENERAL_CONSTANTS::SMTP_SERVER;
my $userName = GENERAL_CONSTANTS::ADMIN_USER_NAME;
my $userPass = GENERAL_CONSTANTS::ADMIN_PASSWORD;
my $mail = "mailto:".GENERAL_CONSTANTS::ADMIN_EMAIL."?subject=Selecton%20Run%20No:%20$run_number";
my $email_subject;
my $email_message;
my $email_system_return;

my $reload_interval = 30;

my $oldRunCommand = ""; # to be filled with the running command line for selecton

my ($H1_likelihood, $H0_likelihood, $AIC_MEC, $AIC_M8a, $sigLevel, $time);
my $seqLen = 0; # used for AIC calculations for MEC

#******************************

open LOG, ">$newCgiLog";
$time = &BIOSEQUENCE_FUNCTIONS::printTime;
print LOG "BEGINNING OF LOG FILE ".$time."\n\n";


&parseCGI();
&runNullModel();
&comparelikelihoods();

open OUTPUT, ">>" .$OutHtmlFile;  
print OUTPUT "\n<H1><center></a></center></H1>\n";
print OUTPUT "<h3><center>Statistical Significance Results for run number <a href=\"$linkToOldResults\">$run_number</a></center></h3>\n\n";	
print OUTPUT  "\nLikelihood for $FORM{MODEL} model allowing for positive selection: $H1_likelihood<br> \n";
my $null_model_name = "M8a";
$null_model_name = "M7" if ($FORM{MODEL} eq "M8b");
print OUTPUT  "\nLikelihood for null model $null_model_name which does not allow positive selection: $H0_likelihood<br><br> \n";
if ($FORM{MODEL} eq "MEC") {
	print OUTPUT "\nAIC_c score for MEC is: $AIC_MEC<br>\n";
	print OUTPUT "\nAIC_c score for M8a is: $AIC_M8a<br><br>\n"
}
print OUTPUT  "\n<b>Likelihood ratio test</b> between the two models shows a significance level of: <b>$sigLevel</b><br><br><br> \n";
print OUTPUT  "\n<b>Ka/Ks results under the null model $null_model_name:</font></b>\n";

open OLD_OUTPUT, "<$oldHtml"; 
my @contents = <OLD_OUTPUT>;
chomp @contents;
close OLD_OUTPUT;
open OLD_OUTPUT, ">$oldHtml"; 
for  (my $line =0; $line < scalar @contents; $line++){
if ($contents[$line] =~ m/statTest/) { # These are the lines which create a "Statistcial Significance" button: advance by 5 lines, and instead print the LRT results
	$line = $line+5;
	if ($sigLevel =~ m/NON-SIGNIFICANT/) {
		print OLD_OUTPUT "\nPositive selection in the protein is NON-significant.<br> \n";
	}
	else {
		print OLD_OUTPUT "\nPositive selection in the protein is significant.<br> \n";
	}
	print OLD_OUTPUT  "\n<a href=nullModel/output.html>Null model results</a><br>\n";

}
print OLD_OUTPUT $contents[$line]."\n";
}
close OLD_OUTPUT;

### print the links    
print OUTPUT "<p><A HREF= $WWWdir".$outputScoreFile."> Codon Ka/Ks scores (numerical values) </A></p>\n";
print OUTPUT "<p><A HREF = $WWWdir".$params.">Likelihood and parameters of Selecton run</A>\n";
print OUTPUT "<p><A HREF = $WWWdir".$log.">Log-file of Selecton run</A>\n";
print OUTPUT "\n<br><br><p><center>Please <a href= \"$mail\">report any problem</a> in case of need.</center></p>\n";

close OUTPUT;

### write to the output that the job has finished
open OUTPUT, "<$OutHtmlFile";
my @output = <OUTPUT>;
close OUTPUT;   
open OUTPUT, ">$OutHtmlFile";
foreach my $line (@output){
	if ($line =~ /Selecton job status page/i){ #finds the phrase "Selecton" job status page, case-insensitive
		print OUTPUT "<H1 align=center>Selecton Job Status Page - <font color='red'>FINISHED</font></h1>\n";
#	    print OUTPUT "<a href=#finish><H2 align=center>Go to the results</font></H2></a>\n";
	}
	else {
		print OUTPUT $line; 
	}
}
close OUTPUT;   
### stop the automatic reload
&stop_reload;
$email_subject = "Your Selecton Statistical Test Results are ready";
$email_message = "Selecton Statistical Test Results are finished. Please click on the following link to see the results:\n".$WWWdir."output.html";
print LOG "Calculation finished. Sending mail to user\n";
chdir $send_email_dir;
$email_system_return = system ('perl sendEmail.pl -f \''.GENERAL_CONSTANTS::ADMIN_EMAIL.'\' -t '.$FORM{email_address}.' -u \''.$email_subject.'\' -xu '.$userName.' -xp '.$userPass.' -s '.$smtp_server.' -m \''.$email_message.'\'');
unless ($email_system_return =~ /successfully/)    {
	print LOG "The message was not sent successfully. system returned: $email_system_return\n";
}

$time = BIOSEQUENCE_FUNCTIONS::printTime();
print LOG "End time: $time\n";


print LOG "\nSelecton Null model run completed successfully!";
print LOG "\n************** END OF LOG FILE *****************\n";
close LOG;

system 'echo "(cd '.$nullDir.' ; chmod -R og+rx * )" | /bin/tcsh';

# by removing this file, the update daemon will know this run has ended successfully
unlink $nullDir.$qsub_ans_file;

exit;
##### MAIN ENDS HERE ######

############################################################
sub read_input{
# ------ storable on -------
    my $input_data = retrieve($runCalcInput);
    %run_data = %$input_data;    
	
	$WWWdir = $run_data{WWWdir}; $linkToOldResults = $run_data{OutputURL}; $outputScoreFile = $run_data{outputScoreFile}; $params = $run_data{params}; $log = $run_data{log}; $oldCgiLog = $run_data{OutLogFile}; $oldHtml = $run_data{OutHtmlFile}; $ErrorDef = $run_data{ErrorDef}; $SysErrorDef = $run_data{SysErrorDef}; $ContactDef = $run_data{ContactDef}; $LogsDir = $run_data{selecton_log_dir}; $formInput = $run_data{FORMInput}; $workingDir = $run_data{WorkingDir}; $nullDir = $run_data{nullDir}; $OutHtmlFile = $run_data{null_html_out};

	my $FORM_data = retrieve($workingDir."formInput.txt");
	%FORM = %$FORM_data;
# ------ storable on -------
	
# ------ storable off -------
#    unless (open INPUT, $workingDir.$runCalcInput) {&error_and_exit("cannot open $runCalcInput $!");}
#    while(<INPUT>){
#        chomp;
#		if(/WWW DIR: (.+)/){$WWWdir = $1;}		
#		elsif(/URL OUTPUT: (.+)/){$linkToOldResults = $1;}
#		elsif(/SCORE RESULTS FILE: (.+)/){$outputScoreFile = $1;}
#		elsif(/GLOBAL RESULTS FILE: (.+)/){$params = $1;}
#		elsif(/KAKS LOG FILE: (.+)/){$log = $1;}
#		elsif(/QSUB LOG: (.+)/){$oldCgiLog = $1;}
#		elsif(/OUTPUT HTML PATH: (.+)/){$oldHtml = $1;}
#		elsif(/ERROR DEF: (.+)/){$ErrorDef = $1;}
#		elsif(/SYS ERROR: (.+)/){$SysErrorDef = $1;}
#		elsif(/CONTACT DEFINITION: (.+)/){$ContactDef = $1;}
#		elsif(/LOG DIR: (.+)/){$LogsDir = $1;}
#	}
#	close INPUT;
# ------ storable off -------
	open ANS, ">".$nullDir.$qsub_ans_file;
    print ANS "OK";
    close ANS;
	chmod 0755, $qsub_ans_file;

}
############################################################

sub parseCGI {
    open OLDCGI, "<$oldCgiLog";

    while ( <OLDCGI> ) {
		if ($_ =~ m/run_kaks4site\: running (.*)$/) { # this is what the running command for Selecton looks like in cgi.log
			$oldRunCommand = $1;
		}    
    }
	close OLDCGI;
    if ($oldRunCommand eq "") {
		&error_and_exit("Error in parseCGI, cannot find the run command for selecton in cgi.log");
    }
    print LOG "parseCGI: The old run command for Selecton was: $oldRunCommand\n";
}
############################################################
sub parseResultsForLikelihood ($) {
    my $globalResultsFile = $_[0];
    my $likelihood;
    unless  (open RESULTS, "<$globalResultsFile") {
	&error_and_exit("Error in parseResultsForLikelihood: cannot open file $globalResultsFile\n");
    }
    while ( <RESULTS> ) {
		if ($_ =~ m/likelihood\: (-\d+\.?\d+)$/) { # the log-likelihood
			$likelihood = $1;
		}
	}
	close RESULTS;
    if ($likelihood eq "") {
		&error_and_exit("Error in parseResults, cannot find likelihood value\n");
    }
    return $likelihood;
}

############################################################
sub runNullModel {
    my $newRunCommand;
    if ($FORM{MODEL} eq "M8") {
		$newRunCommand = $oldRunCommand . " -w1 -Fw"; # command for M8a, setting w=1
    }
    elsif ($FORM{MODEL} eq "MEC") { # run M8a, AIC
		$newRunCommand = $oldRunCommand;
		my $msaFile;
		unless ($oldRunCommand =~ m/-i(.*?)-/) { # -i msaCodonFile -z <matrix>
			&error_and_exit("Error in runNullModel, cannot find input file name from MEC command line $oldRunCommand\n");
		}
		$msaFile = $1;
		my $in  = Bio::SeqIO->new('-file' => "$msaFile" , '-format' => 'Fasta'); ## opening file to get seqLength for later use in AIC
		$seqLen = ($in->next_seq())->length();
		print LOG "\nSeqeunce length = $seqLen\n";
		unless ($seqLen > 0) {
			&error_and_exit("Error in runNullModel: sequnece length = 0");
		}
		$newRunCommand =~ s/mec.exe/selecton.exe/; # change executable to the standard one (not mec)
		$newRunCommand =~ s/-z \d//; ## remove option for empricial matrix
		$newRunCommand .= " -w1 -Fw"; # M8a model
    }
    else { # not supposed to be here
		&error_and_exit("Error in runNullModel, model is not M8 or M8b");
    }
    print LOG "\nrunNullModel: running ".$newRunCommand."\n";

# run Selecton with null model
	chdir $nullDir;
	
	#close LOG;
	`$newRunCommand`;
    #check for user errors in kaks4ite.log
	my $selectonLogFile = $nullDir."kaks4site.log";
	
	if (!(-e $selectonLogFile) or -z $selectonLogFile) {
       &error_and_exit("Unexpected error. Error in runNullModel, $selectonLogFile doesn't exist.");
	}
	print LOG "\nrunNullModel: command finished successfully\n";
}

############################################################

sub comparelikelihoods{
    $H1_likelihood = &parseResultsForLikelihood($workingDir . $params); 
    $H0_likelihood = &parseResultsForLikelihood($nullDir . "globalResult.txt");
    print LOG "sub comparelikelihoods: H1 L= $H1_likelihood, H0 L= $H0_likelihood\n";

    if ($FORM{MODEL} eq "MEC") {
	my $MEC_params = 5;
	my $M8a_params = 4;
	$AIC_MEC = -2 * $H1_likelihood + 2 *  $MEC_params * $seqLen / ($seqLen - $MEC_params - 1);
	$AIC_M8a = -2 * $H0_likelihood + 2 *  $M8a_params * $seqLen / ($seqLen - $M8a_params - 1);
	if ($AIC_MEC < $AIC_M8a) {
	    $sigLevel = "AIC score of MEC is lower than M8a, significance test passed.";
	}
	else {
	    $sigLevel = "AIC score of MEC is higher than M8a, NON-SIGNIFICANT.";
	}
	return;
    }

    my $twice_diff = 2 * ($H1_likelihood - $H0_likelihood) ;
    my %cutoff;

    if ($FORM{MODEL} eq "M8") { # DF = 1: chi-sqaure cutoff values
	$cutoff{0.05} = 3.84;
	$cutoff{0.01} = 6.64;
	$cutoff{0.001} = 10.83;
    }
    elsif ($FORM{MODEL} eq "M8b") { # DF = 2: chi-sqaure cutoff values
	$cutoff{0.05} = 5.99;
	$cutoff{0.01} = 9.21;
	$cutoff{0.001} = 13.82;
    }
    $sigLevel = "NON-SIGNIFICANT";
    for my $key ( reverse sort {$a<=>$b} keys %cutoff ) {
	if ($twice_diff > $cutoff{$key}) {
	    $sigLevel = $key;
	}
    }


}
######################################################################
# Send an error (unexpected error) and exit
sub error_and_exit {
    my $err = $_[0];
    print LOG $err;
    open OUTPUT, ">>$OutHtmlFile";
    print OUTPUT $SysErrorDef;
    print OUTPUT $ContactDef;
    &send_mail();
    &send_mailSelecton($err);
    &stop_reload;
	close LOG;
    exit;
}

######################################################################
# Stops the reload of the output page
sub stop_reload {
    sleep ($reload_interval);
    open OUTPUT, "<$OutHtmlFile";
    my @output = <OUTPUT>;
    close OUTPUT;   
    open OUTPUT, ">$OutHtmlFile";
    foreach my $line (@output){
		unless ($line =~ /REFRESH/ or $line =~ /NO-CACHE/){
			print OUTPUT $line;
		}
    }
    close OUTPUT;
    
    open FINISH, ">".$finish_flag;
    close FINISH;
    # remove the job from the running jobs list
	&BIOSEQUENCE_FUNCTIONS::remove_job_from_running_log("Selecton", "null_$run_number");
}

#####################################################################
# Sends an automatic mail when there are errors
sub send_mail {
	$email_subject = "Subject: Error in Selecton Statistical Testing";
	$email_message = "Unfortunately there was an error while running Selecton Statistical Testing.\n Please click on the following link to see more details\n We apologize for the inconvenience\n\n".$WWWdir."output.html";
	print LOG "send_mail: sending system error to user\n";
	chdir $send_email_dir;
	$email_system_return = system ('perl sendEmail.pl -f \''.GENERAL_CONSTANTS::ADMIN_EMAIL.'\' -t '.$FORM{email_address}.' -u \''.$email_subject.'\' -xu '.$userName.' -xp '.$userPass.' -s '.$smtp_server.' -m \''.$email_message.'\'');
	unless ($email_system_return =~ /successfully/)    {
		print LOG "send_mail: The message was not sent successfully. system returned: $email_system_return\n";
   }
}
#####################################################################
sub send_mailSelecton{
	my $email_message = shift;
	$email_subject = "Subject: Error in Selecton running";
	print LOG "send_mailSelecton: send error message to admin\n";
	chdir $send_email_dir;
	$email_system_return = system ('perl sendEmail.pl -f \''.GENERAL_CONSTANTS::ADMIN_EMAIL.'\' -t '.GENERAL_CONSTANTS::ADMIN_EMAIL.' -u \''.$email_subject.'\' -xu '.$userName.' -xp '.$userPass.' -s '.$smtp_server.' -m \''.$email_message."\n User's email is: $FORM{email_address}\n\'");
	unless ($email_system_return =~ /successfully/) {
		print LOG "The message was not sent successfully. system returned: $email_system_return\n";
   }   
}

#####################################################################
