#!/usr/local/bin/perl

use lib "/bioseq/bioSequence_scripts_and_constants";
use GENERAL_CONSTANTS;

use strict;

my $results_dir = GENERAL_CONSTANTS::SERVERS_RESULTS_DIR;
my $logs_dir = GENERAL_CONSTANTS::SERVERS_LOGS_DIR;
my $bioseq_temp = GENERAL_CONSTANTS::BIOSEQ_TEMP;
my @servers_results_to_remove = ("Selecton", "ConSurf", "ConSeq", "pepitope","Selecton/older_versions", "fastml");

my $report = "/bioseq/data/remove_results_report.txt";

my %months = (Jan => {number => 1},
	       Feb => {number => 2},
	       Mar => {number => 3},
	       Apr => {number => 4},
	       May => {number => 5},
	       Jun => {number => 6},
	       Jul => {number => 7},
	       Aug => {number => 8},
	       Sep => {number => 9},
	       Oct => {number => 10},
	       Nov => {number => 11},
	       Dec => {number => 12});

my $date = `date`;
$date =~ /^\S+\s+(\S+)\s+(\d+)\s+/;
my $month = $1;
my $day = $2;
open REPORT, ">".$report;

foreach (@servers_results_to_remove){    
    remove_server_dirs($_);
}

close REPORT;
if (!-z $report){
    GENERAL_CONSTANTS::send_mail("bioseq", GENERAL_CONSTANTS::ADMIN_EMAIL, "remove", "Remove dirs report", "there was an error while running the script to remove old results. see: $report");	
}
# remove files from bioseq temp directory (mainly created by consurfdb server)
system("rm -rf $bioseq_temp".'*');

sub remove_server_dirs{
    my $server= shift;
    my ($dir_month, $dir_day, $dir_name, $dir_to_remove);
	my $results_list = "/bioseq/data/";
	if ($server =~ /older_versions/){
		$results_list .= "Selecton_older_versions"."_results.txt";}
	else{
		$results_list .= $server."_results.txt";}
	
    chdir $results_dir.$server;
    
    my $cmd = 'ls -ltr | awk \'{print $6, $7, $9}\'> '.$results_list;
    `$cmd`;
    if (-z $results_list){
        print REPORT "no info in $results_list\n";
        return;
    }
    unless (open RES, $results_list){
        print REPORT "cannot open $results_list $!\n";
        `rm -rf $results_list`;
        return;}
    while (<RES>){
        my $rm_flag = 0;
        if (/(\w+) (\d+) (\d{10})/){
            $dir_month = $1;
            $dir_day = $2;
            $dir_name = $3;
			
			chmod 0600, $results_dir.$server."/".$dir_name."/user_email.txt" if (-e $results_dir.$server."/".$dir_name."/user_email.txt");
            remove_system_ER($results_dir.$server."/".$dir_name."/") if ($server eq "pepitope");
            my $month_dif = ($months{$month}{number} - $months{$dir_month}{number}) % 12;
            if ($month_dif==3){
                # should remove all days which are equal or smaller than today.
                if ($day > $dir_day){
                    $rm_flag = 1;
                }
            }
            elsif($month_dif>3){
                $rm_flag = 1;
            }

            if ($rm_flag == 1){
                $dir_to_remove = $results_dir.$server."/".$dir_name;
                system ("rm -rf $dir_to_remove");
                system ("rm -rf ".$logs_dir.$server."/".$dir_name.'*');				
                if (-d $dir_to_remove){
                    print REPORT "the directory: $dir_to_remove was not removed\n";
                }
            }            
        }
    }
    `rm -rf $results_list`;
}

sub remove_system_ER{
	my $work_dir = shift;
	open JOB_NUM ,$work_dir."QSTAT_NO" if (-e $work_dir."QSTAT_NO");
	my $job_id = <JOB_NUM>;
	close JOB_NUM;
	$job_id =~ s/\s*//;	
	my $system_error_file = $work_dir.$job_id.".bioc.ER";
	if ($job_id =~ /^\d+$/ and -e $system_error_file){
		unlink $system_error_file;
	}
}