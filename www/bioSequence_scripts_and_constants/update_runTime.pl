#!/usr/local/bin/perl -w

# read the Currently running list and update the output html file of the run with the time which passed since the beginning of this run.
# in order to know whether this run is still running, check if the file "qsub_ans.txt" exists. This file is removed when a run finishes (by the script which is running on the queue)
use lib "/bioseq/bioSequence_scripts_and_constants";
use GENERAL_CONSTANTS;
use BIOSEQUENCE_FUNCTIONS;
use HandleQueue;
use strict;

my ($qsub_job_no, $server_name, $run_number, $time, $date, $dir_path, @all_lines_in_list, $user_email, $qsub_ans,$visited_index, $ans, $null_dir);
my %que_number_status = ();
my @er_ans = ();


my $running_jobs_list = GENERAL_CONSTANTS::RUNNING_JOBS;
my $log = GENERAL_CONSTANTS::UPDATE_RUN_TIME_LOG_FILE;
my $i = 1;
while ($i<2){ # a process that runs forever    
    open LOG, ">>".$log or die "cannot open $log for writing $!";
    open LIST, "+>>".$running_jobs_list or die "cannot open ".$running_jobs_list." for writing: $!";
    flock LIST, 2;
    seek LIST, 0, 0; #rewind the pointer to the beginning
    @all_lines_in_list = <LIST>; # read the contents into the array
    truncate LIST, 0; # remove all the information, The 0 represents the size of the file that we want
    foreach (@all_lines_in_list){
        chomp;
        if(/^(\d+) (\w+) (null_)?(\d+) (.+:\d+) (.+\-\d+)/){            
            $qsub_job_no = $1;
            $server_name = $2;
            if ($3) {$null_dir = "null_";}
            else {$null_dir = "";}
            $run_number = $4;
            $time = $5;
            $date = $6;
            #print $server_name." $qsub_job_no\n";
            
            $dir_path = GENERAL_CONSTANTS::SERVERS_RESULTS_DIR.$server_name."/$run_number"."/";
            $dir_path .= "nullModel/" if ($null_dir ne "");
            # if the run was deleted by the user - remove it from the list
            if (-e $dir_path."deleted_by_user"){
                $_ = "";
            }
            elsif (-e $dir_path."END_OK"){
                # this means that the run was finished through a normal procedure of the calculation methods
                $_ = ""; # removing this line from the lines array
            }
            elsif (-e $dir_path."qsub_ans.txt"){
                # check if the job is indeed running
                &HandleQueue::find_place_in_Q(\%que_number_status);
                my $continue = "yes";
                if (!exists $que_number_status{$qsub_job_no}){
                    # if the job was not found, we take 10 seconds to make sure it did not ended properly.
                    for(my $time = 0; $time<5; $time++){
                        sleep 2;
                        if (-e $dir_path."END_OK"){
                            $_ = "";
                            $continue = "no";
                            $time=5;
                        }
                    }
                    if ($continue eq "yes"){
                        &HandleQueue::check_job($qsub_job_no, $null_dir.$run_number, $server_name, \@er_ans);    
                        print LOG $er_ans[0];
                        if (exists $er_ans[1]  && $er_ans[1] eq "error"){
                            $_ = "";  # removing this line from the lines array
                            #print "going to report error\n";
                            &report_error($qsub_job_no);
                        }
                    }
                }
                #substructing the time
                else{
                    #print "$qsub_job_no exists in the Q\n";
                    # currently the time report mechnism works only with Selecton
                    if ($server_name eq "Selecton"){
                        $ans = &BIOSEQUENCE_FUNCTIONS::subtract_time_from_now($time, $date);
                        $ans = &GENERAL_CONSTANTS::print_Q_status_in_html($dir_path."output.html", "Running", $ans, "none");
                        print LOG $ans if ($ans ne "OK");
                    }
                }
            }
            else{
                $_ = ""; # removing this line from the lines array
            }
        }        
        $ans = "";
    }
    foreach (@all_lines_in_list){
        if (/.+/){
            print LIST $_."\n";
            #print $_."\n";
        }
    }
    flock LIST, 8;
    close LIST;
    close LOG;
    sleep 5;
    @er_ans = ();
    %que_number_status = ();    
}

sub report_error{
    my $qsub_job_no = shift;
    print LOG "Terminating $run_number. ";    
    $user_email = &HandleQueue::report_error_to_user($dir_path, $server_name, $run_number, $run_number.$null_dir."_Q.log", $run_number, "update_runTime.pl", $dir_path."output.html", $qsub_job_no);
    print LOG "output error message to: ".$dir_path."output.html to email: $user_email\n";
} 
