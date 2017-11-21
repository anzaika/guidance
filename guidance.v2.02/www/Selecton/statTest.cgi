#!/usr/bin/perl

# create HTML file in the nullModel directory
# call the calculation script and run it via the pbs queuing mechanism

use CGI;
use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);

use strict;
use Storable;
use lib "/bioseq/bioSequence_scripts_and_constants";  
use GENERAL_CONSTANTS;
use BIOSEQUENCE_FUNCTIONS;

my $run_number = param('run_number');
my $runCalcInput = param('input_file');
my $Q_log = param('Q_log');
my $user_mail = param('receipent');

my %run_data = ();


my ($WorkingDir, $outHtml, $WWWdir, $walltime);
&read_input;
my $selecton_url = GENERAL_CONSTANTS::SELECTON_URL;

my $user_ip = $ENV{'REMOTE_ADDR'};
&BIOSEQUENCE_FUNCTIONS::check_if_user_is_allowed("selecton",$user_ip, $user_mail);

#print header;
#print start_html("Get Form");
#print "selecton_url is: $selecton_url";
#print end_html;
#exit;

my $nullDir = $WorkingDir . "nullModel/";
my $null_html_out = $nullDir . "output.html"; #new OUTPUT to the user regarding the statistical test run

my $qstat_script = "/bioseq/Selecton/statTest.pl";
my $qsub_script = $WorkingDir."qsub_statTest.sh";
my $Qstat_No_file = "QSTAT_NO_NULL_MODEL";


#------------------------------------------------------------------------------------------------
#-------------------------------------- S C R I P T ---------------------------------------------
#------------------------------------------------------------------------------------------------

my $curr_time = &BIOSEQUENCE_FUNCTIONS::printTime();

# add this run to the log
open LIST, ">>".GENERAL_CONSTANTS::SELECTON_LOG;
flock LIST, 2;
print LIST $curr_time." null_$run_number $user_ip $user_mail\n";
flock LIST, 8;
close LIST;

# print the user's ip and run name in the running jobs list
open RUN_LIST, ">>".GENERAL_CONSTANTS::SELECTON_RUNNING_JOBS;
flock RUN_LIST, 2;
print RUN_LIST "null_$run_number $user_ip $user_mail\n";
flock RUN_LIST, 8;
close RUN_LIST;      

open LOG, ">>".$Q_log;    
print LOG "\n\n-----------------------------------\n\n**** Running statistical Test ****\n\n";

mkdir $nullDir;
&start_output_html;

&store_new_data;
   

# Move directly to the output file
print "Location: ".$WWWdir."output.html\n\n";
close(STDIN); 
close(STDOUT); 
close(STDERR);
##----------------------
#
##print end_html;

# create shell file with the perl script command
# run this on the Q - get the job id
# if something goes wrong: output a message saying the statistical test has an error

unless (open QSUB_SH, ">$qsub_script"){
    &exit_on_system_error("cannot open the file $qsub_script for writing $!\n", $user_mail);
}
print QSUB_SH '#!/bin/sh';
print QSUB_SH "\nperl $qstat_script $run_number $runCalcInput > $WorkingDir"."statTest_std.out";
close QSUB_SH;

chmod 0755, $qsub_script;
print LOG "submitting qsub job \"SELECTON_NULL_$run_number\"\n";

my $qsub_command = "ssh bioseq\@biocluster qsub -q bioseq -e $WorkingDir -o $WorkingDir -N SELECTON_NULL_$run_number $qsub_script";
if ($walltime ne "none"){
    $qsub_command.=" -l walltime=$walltime";
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

unless (open LIST, ">>".GENERAL_CONSTANTS::QUEUING_JOBS){
    print LOG "Could not open file ".GENERAL_CONSTANTS::QUEUING_JOBS.". Reason: $!\nThe job was not listed in the queuing_jobs list.\n";
    &printTime();
}
else{
    flock LIST, 2; # locks the list, so no other process will write to it. On the same time - if the list is currently locked by another process - it waits until the list file is realeased. The "2" and "8" are the operation symbols for "lock" and "unlock".
    print LIST "$qsub_job_no Selecton null_$run_number ".&BIOSEQUENCE_FUNCTIONS::printTime()."\n";
    flock LIST, 8;
    close LIST;
}

print LOG "\n\nCGI ended successfully.";
close LOG;   


exit;

#######################################################################################
sub read_input{
    # --- storable on ---
    my $data = retrieve($runCalcInput);
    %run_data = %$data;
    $WorkingDir = $run_data{WorkingDir};
    $outHtml = $run_data{OutHtmlFile};
    $WWWdir = $run_data{WWWdir}."nullModel/";
    $walltime = $run_data{estimated_run_time};
    # --- storable on ---
    # --- storable off --- 
    #$WorkingDir = GENERAL_CONSTANTS::SERVERS_RESULTS_DIR."Selecton/".$run_number."/";
    #$outHtml = $WorkingDir."output.html";
    #$WWWdir = GENERAL_CONSTANTS::SELECTON_RESULTS_URL.$run_number."/nullModel/";
    # --- storable off ---
}
#######################################################################################
sub store_new_data{
    unlink $runCalcInput;
    $run_data{nullDir} = $nullDir;
    $run_data{null_html_out} = $null_html_out;
    store \%run_data, $runCalcInput;
    if (!(-e $runCalcInput) or (-z $runCalcInput)){
        &exit_on_system_error("\nstore_new_data : Seems that the store failed. the file $runCalcInput doesn't exists or is of size zero\n", $user_mail);
    }
}
#######################################################################################
# Start writing the output web page of Selecton
sub start_output_html {

    system 'echo "(touch '.$null_html_out.'; chmod oug+rxw '.$null_html_out.')" | /bin/tcsh';
    unless (open OUTPUT, ">$null_html_out"){
        &exit_on_system_error("start_output_html : cannot open $null_html_out $!\n", "problems in statTest.cgi: cannot open the file $null_html_out for writing $!");
    }
    print OUTPUT <<EndOfHTML;

<HTML>
<HEAD> <META HTTP-EQUIV="REFRESH" CONTENT=30> </HEAD>
<HEAD> <META HTTP-EQUIV="PRAGMA" CONTENT="NO-CACHE"> </HEAD>


<TITLE>Selecton Statistical Test Results - Run $run_number</TITLE>

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
<tr >
<td align=center><a href="$selecton_url/index.html" id=menu target=_top>HOME</a></td>
<td align=center><a href="$selecton_url/overview.html"  id=menu target=_top>OVERVIEW</a></td>
<td align=center><a href="$selecton_url/gallery.html"  id=menu target=_top>GALLERY</a></td>
<td align=center><a href="$selecton_url/credits.html"  id=menu target=_top>CREDITS</a></td>
</tr>
</table>

</td><tr></table>


<H1 align=center>Selecton Job Status Page</H1>
<H3 align=center>Statistical Test Results</H3>
<blockquote>

<p><font face=Verdana size=2>
Selecton is now calculating the statistical significance of your results.<br>
<b>Please note that this process will take approximately the same time as the original running time.</b><br>
An email will be sent to the address you supplied once the calculation is finished.<br>
This page will be automatically updated every 30 seconds. You can also reload it manually.<br> 
<br><br>

</font></p>


EndOfHTML


    close OUTPUT;

}
#######################################################################################
sub exit_on_system_error{
    my $log_err = shift;
    my $err = shift;
    
    print LOG $log_err;
    &GENERAL_CONSTANTS::send_mail("Selecton", GENERAL_CONSTANTS::ADMIN_EMAIL, $run_number."_nullModel", 'error', $err);
    &BIOSEQUENCE_FUNCTIONS::remove_job_from_running_log("Selecton", "null_$run_number");
    close LOG;
    exit;
}