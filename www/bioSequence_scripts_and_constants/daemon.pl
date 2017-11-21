#!/usr/local/bin/perl -w



#---------------------------------------------------------------------------------------------
#                                **************************   
#                                **  SCRIPT DESCRIPTION  **
#                                **************************
#
#  Daemon responsibility is to make sure that bioseq jobs were not lost; When running a query
# in bioseq server, the first stage is taking place on the ibis server, then the cgi is
# sending the job to the biocluster queue. Daemon is in charge of controlling this transfer.
# Daemon runs forever and performs its checks every 5 seconds.
# 1. Read the "queuing jobs" list, check the status for each job in the biocluster bioseq Q
# 2. Decides which job should be removed from the queuing list and move to the running list:
#    This is done by a flag which should be found in the running directory, in case the job
#    began its run
# 3. If a job is in Q status, output to the screen its position in the Q and the time that
#    passed since the beginning of the run.
# $. If something went wrong: stop the run and output a system error to the user.
#--------------------------------------------------------------------------------------------


# REMARK: $server_name must be spelled:
# Selecton, ConSurf, ConSeq - as it reads info from these directories

# first: reads the bioseq Q. build a hash. key: the job number. value: stats. R or Q_<no>
use lib "/bioseq/bioSequence_scripts_and_constants";
use GENERAL_CONSTANTS;
use BIOSEQUENCE_FUNCTIONS;
use HandleQueue;
use strict;

my $queuing_jobs_list = GENERAL_CONSTANTS::QUEUING_JOBS;
my $running_jobs_list = GENERAL_CONSTANTS::RUNNING_JOBS;
my ($qsub_job_no, $server_name, $run_number, $time, $date, $dir_path, $dir_path_null, @all_lines_in_list, $user_email, $qsub_ans,$visited_index, $ans, $null_dir, $add_null, $cmd);
$user_email = "";
$add_null = "";
my %runs_with_no_qsub = ();
my %que_number_status = ();
my @er_ans = ();
my $ans_out_file = "checkjob.out";
my $ans_err_file = "checkjob.err";


my $i = 1;
while ($i<2){ # a process that runs forever
    # read the queing list
    open LOG, ">>".GENERAL_CONSTANTS::DAEMON_LOG_FILE or die "could not open log GENERAL_CONSTANTS::DAEMON_LOG_FILE $!";
    open LIST, "+>>".$queuing_jobs_list or die "cannot open ".$queuing_jobs_list." for writing: $!";
    flock LIST, 2;
    seek LIST, 0, 0; #rewind the pointer to the beginning
    @all_lines_in_list = <LIST>; # read the contents into the array
    truncate LIST, 0; # remove all the information, The 0 represents the size of the file that we want
    foreach (@all_lines_in_list){
        chomp;
        #------------------------------------
        # getting the info for the queing job
        #------------------------------------
        if(/^(\d+) (\w+) (null_)?(\d+) (.+:\d+) (.+\-\d+)( visit\d+)?/){
            $qsub_job_no = $1;
            $server_name = $2;
            if ($3) {$null_dir = "null_";}
            else {$null_dir = "";}
            $run_number = $4;
            $time = $5;
            $date = $6;
                    
            $dir_path = GENERAL_CONSTANTS::SERVERS_RESULTS_DIR.$server_name."/$run_number"."/";
            $dir_path_null = $dir_path;
            $dir_path_null = $dir_path ."nullModel/" if ($null_dir ne "");
            $add_null = "/nullModel" if ($null_dir ne "");
            #---------------------------------------
            # read the status of the Q in biocluster
            #---------------------------------------
            &HandleQueue::find_place_in_Q(\%que_number_status);
            if (exists $que_number_status{$qsub_job_no}){
                #---------------------------------
                # check if the job started running
                #---------------------------------                
                if (-e $dir_path_null."qsub_ans.txt"){
                    open QSUB_ANS, $dir_path_null."qsub_ans.txt";             
                    $qsub_ans = <QSUB_ANS>;
                    close QSUB_ANS;
                    if ($qsub_ans eq "OK"){
                        #print "removing from list: $qsub_job_no $run_number and move it to run list\n";
                        $_ = "";  # removing this line from the lines array
                        unless (open RUN_LIST, ">>".$running_jobs_list)
                            {print LOG "Could not open the file $running_jobs_list. reason: $!. report run $run_number is running\n";}
                        flock RUN_LIST, 2;
                        print RUN_LIST "$qsub_job_no $server_name ".$null_dir."$run_number $time $date\n";            
                        flock RUN_LIST, 8;
                        close RUN_LIST;
                    }
                    else{
                        $_ = "";  # removing this line from the lines array
                        &report_error($qsub_job_no);
                    }            
                }
                #----------------------------------------------------------------------
                # jobs is in the Q, in Q status
                # count the number of the daemon visits to it (not active at the moment)
                #----------------------------------------------------------------------
                else{
                    # if the job was not accepted yet, we output the Q status to the user
                    if (exists $runs_with_no_qsub{$run_number}){                    
                        $runs_with_no_qsub{$run_number}++;
                        $visited_index = $runs_with_no_qsub{$run_number};
                        #$_ =~ s/visit\d+/visit$visited_index/;
                        #$_.= "\n";
                    }
                    else{
                        $runs_with_no_qsub{$run_number} = 1;
                        #$_.= " visit1\n";
                    }
                    #print LOG "going to substruct time from $time $date\n";
                    $ans = &BIOSEQUENCE_FUNCTIONS::subtract_time_from_now($time, $date);
                    #print LOG $ans if ($ans ne "OK");
                    #print LOG "reply from BIOSEQUENCE_FUNCTIONS::subtract_time_from_now  :  $ans\n";
                    $ans = &GENERAL_CONSTANTS::print_Q_status_in_html($dir_path."output.html", $que_number_status{$qsub_job_no}, $ans, "none"); 
                    print LOG $ans if ($ans ne "OK");
                    #print LOG "reply from GENERAL_CONSTANTS::print_Q_status_in_html  :  $ans\n";
                }
            }                                    
            #---------------------
            # jobs is not in the Q
            #---------------------            
            elsif(!(exists $que_number_status{$qsub_job_no})){
                sleep 7; # maybe the run was just finished and the out files were not created yet.
                #-------------------------------------------------------------
                # if the run was deleted by the user - remove it from the list
                #-------------------------------------------------------------
                if (-e $dir_path."deleted_by_user"){
                    $_ = "";
                }
                #-----------------------------------             
                # check if it was not ended properly
                #-----------------------------------
                elsif (-e $dir_path_null."END_OK"){
                    # this means that the run was finished through a normal procedure of the calculation methods
                    $_ = "";  # removing this line from the lines array
                }
                #---------------------------------------------------------------------------------------------
                #if the job is not running, there will be an error when trying to check for it using "checkjob"
                #---------------------------------------------------------------------------------------------
                else{
                    &HandleQueue::check_job($qsub_job_no, $null_dir.$run_number, $server_name, \@er_ans);    
                    print LOG $er_ans[0];
                    if (exists $er_ans[1] && $er_ans[1] eq "error"){
                        $_ = "";  # removing this line from the lines array
                        &report_error($qsub_job_no);
                    }
                }
            }            
        }
    }
    # after reading the file, reprint the content to the list, only for lines which were not removed by the daemon
    foreach (@all_lines_in_list){
        if (/.+/){
            print LIST $_."\n";
        }
    }    
    flock LIST, 8;
    close LIST;
    close LOG;
    sleep 5;
    $add_null = "";
    %que_number_status = ();
    @er_ans = ();
}

#----------------------------------
## lock the HTML file, updating its status in the Q
#sub print_Q_status_in_html{
#    my $html_file = shift;
#    my $_status = shift;
#    $html_file.="output.html";  
#    
#    unless (open HTML, "+>>".$html_file) {
#        print LOG "Could not open file $html_file to update the status. Status is: $_status\n";}
#    else{
#        flock HTML, 2;
#        seek HTML, 0, 0; #rewind the pointer to the beginning
#        my @html_lines = <HTML>; # read the contents into the array
#        truncate HTML, 0; # remove all the information, The 0 represents the size of the file that we want
#        foreach (@html_lines){    
#            if(/Your job status is: (.+)<\/font><br>/){
#                s/$1/$_status/;
#                last;
#            }    
#        }
#        print HTML $_ foreach (@html_lines);
#        flock HTML, 8;
#        close HTML;
#    }
#}
#----------------------------------
sub report_error{
    my $qsub_job_no = shift;
    print LOG "Terminating $run_number. ";
    $user_email = &HandleQueue::report_error_to_user($dir_path, $server_name, $run_number.$add_null, $run_number.$null_dir.".log", $run_number, "daemon.pl", $dir_path_null."output.html", $qsub_job_no);
    print LOG "output error message to: ".$dir_path_null."output.html to email: $user_email\n";    
}
