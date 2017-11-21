use strict;

##############################################
# MONITOR THE BIOSEQ NODE - RUN ALL THE TIME #
##############################################

# first: reads the bioseq Q. build a hash. key: the job number. value: stats. R or Q_<no>
use lib "/bioseq/bioSequence_scripts_and_constants";
use GENERAL_CONSTANTS;
my $sending_script="Monitor_Bioseq_Node.pl";
my $Bioseq_Node=GENERAL_CONSTANTS::BIOSEQ_NODE;
my $Running_Jobs_Log=GENERAL_CONSTANTS::JOBS_ON_BIOSEQ_NODE;
my $Get_Runnings_Jobs_Command="ssh bioseq\@biocluster.tau.ac.il ssh $Bioseq_Node ps -ef | grep \"bioseq\" |";
my $Log_File="/groups/bioseq.home/test_Monitor.log";


my $i = 1;
while ($i<2){ # a process that runs forever
#open (LOG,">>$Log_File");
#print "$Get_Runnings_Jobs_Command\n";
open (RUNNING_JOBS,$Get_Runnings_Jobs_Command) || die "Can't exec '$Get_Runnings_Jobs_Command' $!";
my %Running_Jobs=();
my $Server_Name;
my $Run_Num;
my $dir_path;
while (my $line=<RUNNING_JOBS>)
{
#bioseq   21151 21150  0 15:20 ?        00:00:00 perl /bioseq/ConSurf_old/consurf_run_calc.pl /bioseq/data/results/ConSurf/1246882804/ input.data
	if ($line=~/bioseq\s+([0-9]+) ([0-9]+).*([0-9]+)(.*)(Con.*)\/([0-9]+).*/)
	{
		my $Job_PID=$1;
		my $Job_Name=$6;
		my $Server_Name=$5;
		$Running_Jobs{$Job_Name}{'SERVER'}=$Server_Name;
		$Running_Jobs{$Job_Name}{'PID'}=$Job_PID;
#		print LOG "JOB_PID:*$Job_PID*\tNUM:*$Job_Name*\tSERVER:*$Server_Name*\n";
	}
#bioseq   21684 21673  0 16:09 ?        00:00:00 /bin/sh /bioseq/data/results//patchfinder/1246972158/patchfinder	
	elsif ($line=~/bioseq\s+([0-9]+) ([0-9]+).*([0-9]+)(.*)([0-9]+)\/(patchfinder).*/)
	{
		my $Job_PID=$1;
		my $Job_Name=$5;
		my $Server_Name=$6;
		$Running_Jobs{$Job_Name}{'SERVER'}="patchfinder";
		$Running_Jobs{$Job_Name}{'PID'}=$Job_PID;
	}
}
close (RUNNING_JOBS);
open (LIST,"+>>".$Running_Jobs_Log) || die "Can't Open the Runnings_Jobs Log '$Running_Jobs_Log' $!";

flock LIST, 2;
seek LIST, 0, 0; #rewind the pointer to the beginning
my @all_lines_in_list = <LIST>; # read the contents into the array
my @new_list=();
truncate LIST, 0; # remove all the information, The 0 represents the size of the file that we want
	foreach my $line (@all_lines_in_list){
       		chomp;
		#ConSurf 1246803545 17:19:13 05-07-2009
		if ($line=~/([A-Za-z]+) ([0-9]+) (.*)/)
		{
			$Server_Name=$1;
			$Run_Num=$2;
#			print LOG "SERVER:$Server_Name\tRUN:$Run_Num\t";
			$dir_path = GENERAL_CONSTANTS::SERVERS_RESULTS_DIR.$Server_Name."/$Run_Num"."/";
			if ($Running_Jobs{$Run_Num}{'SERVER'} eq $Server_Name) #Job Is Running
			{
				push (@new_list,$line);
				#print LIST "$line\n";
#				print LOG "FOUND $Run_Num - $Server_Name\n";
			}
			elsif (!-e $dir_path."END_OK") #Not Endded OK
			{
#				print LOG "NOT FOUND $Run_Num - $Server_Name\n";
				sleep 20; # maybe the run was just finished and the out files were not created yet.
				if (!-e $dir_path."END_OK")
				{
					&report_error($Server_Name,$Run_Num);
					#push (@new_list,$line);
				}
			}
		}
		else
		{
			if ($line =~/[A-Z]+/)
			{
				push (@new_list,$line);
				#print LIST "$line\n";
			}
		}
	}
	foreach my $line (@new_list)
	{
		print LIST "$line\n";
	}
flock LIST, 8;
close LIST;
sleep 5;
#close LOG;

}

sub report_error{
    my $Server_Name=shift;
    my $Run_Num=shift;
    my $Logs_Dir=GENERAL_CONSTANTS::SERVERS_LOGS_DIR;
    my $err_subject = "Error in '$Server_Name' run '$Run_Num' on $Bioseq_Node";
    my $err_message = "*** MESSAGE FROM $sending_script : an error occured while trying to run the job \'$Run_Num\' in BioSeq Node $Bioseq_Node. The run was stopped. Have A look at: $Logs_Dir"."$Server_Name"."/"."$Run_Num".".log ***\n";
    GENERAL_CONSTANTS::send_mail($Server_Name, GENERAL_CONSTANTS::ADMIN_EMAIL, $Run_Num, $err_subject, $err_message);
}
