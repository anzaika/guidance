#!/usr/bin/perl -w
#use lib "/bioseq/Guidance/";

use strict;
use Storable;
use FindBin qw($Bin);
use lib $Bin;

use Guidance;
use Getopt::Long;
use File::Basename;
if (scalar @ARGV<2) {die "USAGE: perl $0 --MSA <Base MSA> --Scores <Scores_File> --FilterdMSA <Out_MSA_File> --Cutoff <Cutoff> --RemovedPos <Out_File_With_Removed_Pos>\n";}
my ($stored_data_file,$Cutoff,$vars_ref,$isServer);
my %VARS=();
if ($ARGV[0] =~ m/^-/) { # cmd_line mode
	if (scalar @ARGV<5) {die "USAGE: perl $0 --MSA <Base MSA> --Scores <Scores_File> --FilterdMSA <Out_MSA_File> --Cutoff <Cutoff> --RemovedPos <Out_File_With_Removed_Pos>\n";}
	# Commandline input mode
	my $getoptResult = GetOptions ("MSA=s"=>\$VARS{Alignment_File},  # = means that this parameter is required, s means string
								   "Scores=s"=>\$VARS{Col_Scores_File},
								   "FilterdMSA=s"=>\$VARS{Alignment_File_without_low_SP_Col},
								   "Cutoff=f"=>\$Cutoff,
                                   "RemovedPos=s"=>\$VARS{removed_low_SP_SITE},
		);								   
	$isServer="NO";
}
else # server mode
{
	$stored_data_file=shift;
	$Cutoff=shift;
	$vars_ref = retrieve($stored_data_file);
	%VARS = %$vars_ref;
	$VARS{Alignment_File_without_low_SP_Col}=$VARS{WorkingDir}."/".$VARS{Alignment_File_without_low_SP_Col}.".$Cutoff";
	$VARS{Col_Scores_File}="$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_col.scr";
	$VARS{Alignment_File}="$VARS{WorkingDir}$VARS{Alignment_File}";           # add path
	$VARS{removed_low_SP_SITE}="$VARS{WorkingDir}$VARS{removed_low_SP_SITE}.$Cutoff"; # add path
	$isServer="YES";
}
open (LOG,">>$VARS{OutLogFile}") || exit_on_error('sys_error', "Can't open Log File: $VARS{OutLogFile} $!") if ($isServer eq "YES");

#remove sites with SP-score < Col sp_cutoff
############################################
print LOG "Guidance::removeLowSPsites_NoBioPerl (\"$VARS{Alignment_File}\",\"$VARS{Col_Scores_File}\",\"$VARS{Alignment_File_without_low_SP_Col}\",$Cutoff,\"$VARS{removed_low_SP_SITE}\");\n" if ($isServer eq "YES");
my @ans=Guidance::removeLowSPsites_NoBioPerl ($VARS{Alignment_File},$VARS{Col_Scores_File},$VARS{Alignment_File_without_low_SP_Col},$Cutoff,$VARS{removed_low_SP_SITE});
if ($ans[0]eq "OK")
{
    $VARS{REMOVED_SITES}=$ans[1];
    $VARS{MSA_LENGTH}=$ans[2];
}

if ($isServer eq "NO")
{
	print "REMOVED_SITES:$VARS{REMOVED_SITES}\n";
	print "MSA_LENGTH:$VARS{MSA_LENGTH}\n";
}
else # server
{
	print LOG "REMOVED_SITES:$VARS{REMOVED_SITES}\n";
	print LOG "MSA_LENGTH:$VARS{MSA_LENGTH}\n";
	$VARS{Alignment_File_without_low_SP_Col_with_Names}=$VARS{Alignment_File_without_low_SP_Col}.".With_Names";
	my $Alignment_File_without_low_SP_Col_with_Names_NO_PATH=basename($VARS{Alignment_File_without_low_SP_Col_with_Names});
	my $removed_low_SP_SITE_NO_PATH=basename($VARS{removed_low_SP_SITE});
	if (-s $VARS{Alignment_File_without_low_SP_Col} > 0) # Not EMPTY
	{
		my @ans=Guidance::codes2nameFastaFrom1($VARS{Alignment_File_without_low_SP_Col},"$VARS{WorkingDir}$VARS{code_fileName}",$VARS{Alignment_File_without_low_SP_Col_with_Names});
		if ($ans[0] ne "OK") {exit_on_error("sys_error","Guidance::codes2nameFastaFrom1: Guidance::codes2nameFastaFrom1(\"$VARS{Alignment_File_without_low_SP_Col}\",\"$VARS{WorkingDir}$VARS{code_fileName}\",\"$VARS{Alignment_File_without_low_SP_Col_with_Names}\") failed:",join("",@ans),"\n");}
	}
    # Update the output page
    #######################################
	open (OUTPUT,"$VARS{WorkingDir}$VARS{output_page}");
	my @out=<OUTPUT>;
	close (OUTPUT);
	open (OUTPUT,">$VARS{WorkingDir}$VARS{output_page}");
	my $Remove_Pos_Section=0;
	foreach my $line (@out)
	{
		if ($line=~/Remove unreliable columns below confidence score/)
		{
			$Remove_Pos_Section=1;
			print OUTPUT $line;
		}
		elsif (($line=~/form/) and $Remove_Pos_Section==1)
		{
			print OUTPUT $line;
			print_message_to_output("<A HREF='$Alignment_File_without_low_SP_Col_with_Names_NO_PATH' TARGET=_blank>The MSA after the removal of unreliable columns (below $Cutoff)</A><font size=-1> (see list of removed columns <A HREF='$removed_low_SP_SITE_NO_PATH' TARGET=_blank>here</A>)</font><br>"); 
			$Remove_Pos_Section=0;
		}
		else
		{
			print OUTPUT $line;
		}
	}
	close (OUTPUT);
}

#---------------------------------------------
sub print_message_to_output{
#---------------------------------------------
    my $msg = shift;
    print OUTPUT "\n<ul><li>$msg</li></ul>\n";
}
