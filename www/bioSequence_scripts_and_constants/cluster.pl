
#!/usr/bin/perl

#  --------
#  PROGRAM:  1lock.pl  (created by Developer's Daily, http://www.DevDaily.com)
#  --------
use lib "/bioseq/bioSequence_scripts_and_constants";  
use GENERAL_CONSTANTS;
use Fcntl qw(:DEFAULT :flock);
#$running_jobs_list = "list.txt";

$script_order = shift;

$LOCK_EXCLUSIVE = 2;
$UNLOCK         = 8;

# -----------------------
# What's about to happen:
# -----------------------
# open the file, lock the file, sleep, then write, then unlock the
# file, then close the file.

open (FILE, ">>".GENERAL_CONSTANTS::RUNNING_JOBS) || die "problem opening ".GENERAL_CONSTANTS::RUNNING_JOBS."\n";

#open (FILE, ">> $running_jobs_list") || die "problem opening $running_jobs_list\n";
#flock FILE, $LOCK_EXCLUSIVE or die "cluster return: $!";
print FILE "cluster $script_order entered\n";
sleep 7;
print FILE "this line printed by cluster.pl $script_order\n";
#flock FILE, $UNLOCK or die "cluster return: $!";
close(FILE);
