#!/usr/bin/perl -w

use lib "/bioseq/bioSequence_scripts_and_constants"; #"/bioseq/scripts_for_servers";  
use GENERAL_CONSTANTS;
use BIOSEQUENCE_FUNCTIONS;

use CGI qw(:standard);
use strict;

my $queryForm = new CGI;
my ($line, $Q_no);

# 3 options: (1)the query was not submitted to the Q, (2)The query is in status Q, (3)The query is in status R
# if (1) : kill the cgi process and all its children
# if (2) or (3): delete the job from the Q
# delete all the files from the running dir
# in both cases: create a new file of the name: deleted_by_user, this will be a flag to the daemon and the update scripts that they should remove this run from their lists.
# load the cancel html page


my $Qstat_no_file = $queryForm->param("Qstat_file");
my $selecton_hp = $queryForm->param("selecton_http");
my $working_dir = $queryForm->param("run_no");
my $cgi_pid = $queryForm->param("cgi_pid");

my $run_name;


# first - try to remove cgi and all child processes of it.
# than - try to remove from Q (might also be in Q)

my $cgi_cmd = `ps $cgi_pid`;
chomp($cgi_cmd);
# if the process is still alive on ibis - it will be written to the ps list
if ($cgi_cmd =~ $cgi_pid){
    my $pid_to_kill .= $cgi_pid." ";
    my $pid_found = "yes";
    my $new_pid;
    while ($pid_found eq "yes"){
        $pid_found = "no";
        $new_pid = `ps U bioseq -o pid,ppid | awk \'\{if \(\$2\~\"$cgi_pid\"\) print \$1\}\'`;
        chomp($new_pid);
        if ($new_pid =~ /\d+/){
            $pid_found = "yes";
            $pid_to_kill .= $new_pid." ";
            $cgi_pid = $new_pid;
        }
    }
    chomp($pid_to_kill);
    `kill $pid_to_kill`;
}

if (-e $working_dir.$Qstat_no_file){
    if (open QSTAT, $working_dir.$Qstat_no_file){
        $line = <QSTAT>;
        if ($line =~ m/^(\d+).+/){
            $Q_no = $1;
        }
        my $cmd = "ssh bioseq\@biocluster qdel $Q_no";
        `$cmd`; 
    }
}

if ($working_dir =~ /\d+\/$/){
    system("rm -f $working_dir".'*');
    open DEL, ">".$working_dir."deleted_by_user";
    close DEL;
}

# remove the line from the running jobs list
if ($working_dir =~ m/\/(\d+)\//){
    $run_name = $1;
    # remove the job from the running jobs list
    &BIOSEQUENCE_FUNCTIONS::remove_job_from_running_log("Selecton", $run_name);
}

print "Location: $selecton_hp\n\n";
exit;


