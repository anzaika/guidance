#!/usr/bin/perl

package HandleQueue;

#use GENERAL_CONSTANTS;
use BIOSEQUENCE_FUNCTIONS;
use strict;
#----------------------------------
# read the Q and assign each job with a status
sub find_place_in_Q{
    my $hash_ref = shift;
    my $cmd = 'qstat bioseq';    
    my $num_in_Q = 0;
    my $num_of_R = 0;
    my $status = "";
    my (@line, $pbs_no);
    
    my @ans = `$cmd`;
    foreach(@ans){
        chomp;
        $num_in_Q++;    
        @line = split(/\s+/, $_);
        if($line[0] =~ /^(\d+)\./){
            $pbs_no = $1 ;
            if($line[4] =~ /R/) {$num_of_R++; $status="Running";}
            elsif($line[4] =~ /E/) {$num_of_R++; $status="Running";}
            elsif($line[4] =~ /Q/) {$status="Queued. Your number in the queue is: ".($num_in_Q-$num_of_R);}
            $hash_ref->{$pbs_no} = $status;
        }
        $pbs_no = "";
        $status="";
    }    
}
#------------------------------------------------------------------
sub check_job{
    my $qsub_job_no= shift;
    my $run_number = shift;
    my $server_name = shift;
    my $ans_ref = shift;
    my $ans_out_file = $qsub_job_no."_checkjob.out";
    my $ans_err_file = $qsub_job_no."_checkjob.err";
    
    my $cmd = "sh -c 'checkjob $qsub_job_no 1>$ans_out_file 2>$ans_err_file'";
    `$cmd`;
    if (-z $ans_err_file and !(-z $ans_out_file)){
        # if the file exists, it means that this job really lives.
        $ans_ref->[0] = "The job $qsub_job_no of the run $server_name $run_number was not found under \"qstat bioseq\", but created a file with \"checkjob\"\n";
        # should be written to the Q next round
    }
    #if the error file is not of size zero - then this job is dead
    elsif(!(-z $ans_err_file) and -z $ans_out_file){
        $ans_ref->[0] = "The job $server_name $run_number was not ended properly. Was not find in the biocluster Q.  ";
        $ans_ref->[1] = "error";
    }
    else{
        $ans_ref->[0] = "The job $qsub_job_no of the run $server_name $run_number was not found under \"qstat bioseq\", and did not create any file with \"checkjob\"\n";
        $ans_ref->[1] = "error";
    }
    unlink $ans_err_file;
    unlink $ans_out_file;
}
#------------------------------------------------------------------
sub report_error_to_user{
    
    my $dir_path = shift;
    my $server_name = shift;
    my $run_dir = shift;
    my $log_run = shift;
    my $run_number = shift;    
    my $sending_script = shift;
    my $output_path = shift;
    my $qsub_job_no = shift;
    
    my $user_email = "";
    my ($err_message, $err_subject);

    # not ok: we write it to HTML, stop reload and send an e-mail to user (if supplied)
    if (-e $dir_path."user_email.txt"){
        if (open MAIL, $dir_path."user_email.txt"){
            $user_email = <MAIL>;
            chomp($user_email);
            close MAIL;                        
        }                    
    }
    if ($user_email eq "NOT_GIVEN" or $user_email eq ""){$user_email = "NO";} 
    # print to output that there was a failoure and send mail to user
    #GENERAL_CONSTANTS::print_to_output($output_path, $server_name, $run_dir, $user_email);
    # remove the job from the running list of the server
    &BIOSEQUENCE_FUNCTIONS::remove_job_from_running_log($server_name, $run_number);
    
    $err_message = "*** MESSAGE FROM $sending_script : an error occured while trying to run the job in biocluster Q. The run was stopped.***";
    if (open RUN_LOG, ">>".GENERAL_CONSTANTS::SERVERS_LOGS_DIR.$server_name."/$log_run"){        
        print RUN_LOG "\n".$err_message."\n";
        close RUN_LOG;
    }
    $err_subject = "Error in $server_name run $run_number";
    if (defined $qsub_job_no){
        $err_subject .= " job $qsub_job_no";
        $err_message .= "\njob $qsub_job_no";    
        # check if the job died for no reason
        my $err_file = $dir_path.$qsub_job_no.".bioc.ER";
        my $out_file = $dir_path.$qsub_job_no.".bioc.OU";
        if (-e $err_file and -z $err_file and -e $out_file and -z $out_file){
            my $host = find_host($qsub_job_no);
            $err_subject .= " died on host $host";
            $err_message .= " died on host $host";
        }
        else{
            $err_message .= "\n\nDid not file $err_file";
        }
    }
    # report about the error to the administrator
    GENERAL_CONSTANTS::send_mail($server_name, GENERAL_CONSTANTS::ADMIN_EMAIL, $run_number, $err_subject, $err_message."\nRun: $run_number User: $user_email\n");
    return $user_email;
}

sub find_host{
    my $job_no = shift;
    my $n = 1;
    my $found = 0;
    my $ret = "";
    
    my $cmd = "tracejob $job_no | grep exec_host";
    while ($found != 1 and $n<10){                        
        my $ans = `$cmd`;
        if ($ans =~ /exec_host=(bioc\d+)\.tau\.ac\.il/){
            $ret = $1;
            #print "$job_no died on $1\n";
            $found = 1;
        }
        else{
            $n++;
            $cmd = "tracejob -n $n $job_no | grep exec_host";
            #print "ans was: $ans\n";
        }
    }
    return $ret;
}
#------------------------------------------------------------------
# return the number of running jobs on specific node on the cluster
sub node_status
{
	my $node=shift;
	my $command="ssh bioseq\@biocluster ssh $node  ps -ef | awk \'\$4!=0\' | grep \"bioseq\" | wc -l|";
	unless (open (NUM_OF_RUNNING,$command)) {return "HandleQueue::node_status Can't Execute \'$command\': $!";}
	my $running_jobs=<NUM_OF_RUNNING>;
	close (NUM_OF_RUNNING);
	chomp $running_jobs;
	#print "running: $running_jobs\n";
	return ("ok",$running_jobs);
}
#------------------------------------------------------------------
# return the number of running jobs on the entire cluster
sub queue_status
{
	unless (open (NUM_OF_RUNNING,"ssh bioseq\@biocluster qstat | grep -v bioseq | grep -c R|")) {return  "HandleQueue::queue_status Cna't Execute $!";}
	my $Num_Of_Jobs=<NUM_OF_RUNNING>;
	chomp ($Num_Of_Jobs);
	close (NUM_OF_RUNNING);
	return ("ok",$Num_Of_Jobs);
}
1;
