#!/usr/bin/perl -w
#use lib "/bioseq/Guidance/";

use strict;
use Storable;
use FindBin qw($Bin);
use lib $Bin;

use Guidance;
use Getopt::Long;
use File::Basename;

my ($Cutoff,$isServer,$Type);
my %VARS=();
my %FORM=();
if (scalar @ARGV<3) {die "USAGE: perl $0 --MSA <Base MSA> --Scores <Scores_File> --FilterdSeq <Out_Seq_File> --Cutoff <Cutoff> --RemovedSeq <Out_File_With_Removed_Seq>\n";}
if ($ARGV[0] =~ m/^-/) { # cmd_line mode
	if (scalar @ARGV<5) {die "USAGE: perl $0 --MSA <Base MSA> --Scores <Scores_File> --FilterdSeq <Out_Seq_File> --Cutoff <Cutoff> --RemovedSeq <Out_File_With_Removed_Seq> --Type <BySeqName|ByRowNum [optioal]>\n";}
	# Commandline input mode
	my $getoptResult = GetOptions ("MSA=s"=>\$VARS{Alignment_File},  # = means that this parameter is required, s means string
								   "Scores=s"=>\$VARS{Seq_Scores_File},
								   "FilterdSeq=s"=>\$VARS{Seq_File_without_low_SP_SEQ},
								   "Cutoff=f"=>\$Cutoff,
                                   "RemovedSeq=s"=>\$VARS{removed_low_SP_SEQ},
								   "Type:s"=>\$Type, # BySeqName | ByRowNum (optional)
		);								   
	$isServer="NO";
	if ($Type eq ""){$Type="BySeqName";}
}
else
{
	my $stored_data_file=shift;
	my $stored_form_file=shift;
	$Cutoff=shift;

	my $vars_ref = retrieve($stored_data_file);
	%VARS = %$vars_ref;
	
	my $form_ref = retrieve($stored_form_file);
	%FORM = %$form_ref;


	$VARS{Seq_Scores_File}="$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_seq.scr";

	$VARS{Seq_File_without_low_SP_SEQ}=$VARS{WorkingDir}."/".$VARS{Seq_File_without_low_SP_SEQ}.".$Cutoff"; # add path and cutoof suffix
	$VARS{removed_low_SP_SEQ}=$VARS{WorkingDir}."/".$VARS{removed_low_SP_SEQ}.".$Cutoff"; # add path and cutoof suffix
	
	$isServer="YES";
	$Type="ByRowNum";
}
open (LOG,">>$VARS{OutLogFile}") || exit_on_error('sys_error', "Can't open Log File: $VARS{OutLogFile} $!") if ($isServer eq "YES");

#remove sites with SP-score < Col sp_cutoff
############################################

print LOG "Guidance::removeLowSPseq ($VARS{Alignment_File},$VARS{Seq_Scores_File},$VARS{Seq_File_without_low_SP_SEQ},$Cutoff,$VARS{removed_low_SP_SEQ},$Type);\n" if ($isServer eq "YES");
my @ans=Guidance::removeLowSPseq ($VARS{Alignment_File},$VARS{Seq_Scores_File},$VARS{Seq_File_without_low_SP_SEQ},$Cutoff,$VARS{removed_low_SP_SEQ},$Type);
print "ANS:",join("",@ans),"\n";

if ($isServer eq "YES")
{
	$VARS{Seq_File_without_low_SP_SEQ_with_Names}=$VARS{Seq_File_without_low_SP_SEQ}.".With_Names";
	$VARS{removed_low_SP_SEQ_With_Names}=$VARS{removed_low_SP_SEQ}.".With_Names";
	if (-s  $VARS{Seq_File_without_low_SP_SEQ} > 0) # NOT EMPTY
	{
		my @ans=Guidance::codes2nameFastaFrom1($VARS{Seq_File_without_low_SP_SEQ},"$VARS{WorkingDir}$VARS{code_fileName}",$VARS{Seq_File_without_low_SP_SEQ_with_Names});
		if ($ans[0] ne "OK") {exit_on_error("sys_error","Guidance::codes2nameFastaFrom1: Guidance::codes2nameFastaFrom1($VARS{Seq_File_without_low_SP_SEQ},\"$VARS{WorkingDir}$VARS{code_fileName}\",$VARS{Seq_File_without_low_SP_SEQ_with_Names}) failed:",join("",@ans),"\n");}
	}
	
	if (-s $VARS{removed_low_SP_SEQ} > 0) # Seq were removed
	{
		my  @ans=Guidance::codes2nameFastaFrom1($VARS{removed_low_SP_SEQ},"$VARS{WorkingDir}$VARS{code_fileName}",$VARS{removed_low_SP_SEQ_With_Names});
		if ($ans[0] ne "OK") {exit_on_error("sys_error","Guidance::codes2nameFastaFrom1: Guidance::codes2nameFastaFrom1(Guidance::codes2nameFastaFrom1($VARS{removed_low_SP_SEQ},\"$VARS{WorkingDir}$VARS{code_fileName}\",$VARS{removed_low_SP_SEQ_With_Names}) failed:".join("",@ans)."\n");}
	}
	
	
    # Update the output page
    #######################################
	open (OUTPUT,"$VARS{WorkingDir}$VARS{output_page}");
	my @out=<OUTPUT>;
	close (OUTPUT);
	open (OUTPUT,">$VARS{WorkingDir}$VARS{output_page}");
	my $Seq_File_without_low_SP_SEQ_with_Names_NO_PATH=basename($VARS{Seq_File_without_low_SP_SEQ_with_Names});
	my $removed_low_SP_SEQ_With_Names_NO_PATH=basename($VARS{removed_low_SP_SEQ_With_Names});

	my $Remove_Seq_Section=0;
	foreach my $line (@out)
	{
		if ($line=~/Remove unreliable sequences below confidence score/)
		{
			$Remove_Seq_Section=1;
			print OUTPUT $line;
		}
		elsif (($line=~/form/) and $Remove_Seq_Section==1 and ($line!~/form.data/))
		{
			print OUTPUT $line;
			if ($FORM{'Redirect_From_MAFFT'}==1) {print_message_to_output("<A HREF='$Seq_File_without_low_SP_SEQ_with_Names_NO_PATH' TARGET=_blank>The input sequences after the removal of unreliable sequences (with confidence score below $Cutoff)</A><font size=-1> (see list of removed sequences <A HREF='$removed_low_SP_SEQ_With_Names_NO_PATH' TARGET=_blank>here</A></font>)&nbsp;&nbsp;&nbsp;<INPUT TYPE=\"BUTTON\" VALUE=\"run GUIDANCE on the confidently-aligned sequences only\" ONCLICK=\"var answer = confirm('ATTENTION: Running GUIDANCE on the confidently-aligned sequences only, ignores the parameters used for the original run on MAFFT server. It is therefore recommended to adjust these parameters or aligning the confidently-aligned sequences on MAFFT server and run GUIDANCE again from there');if (answer){window.open('http://guidance.tau.ac.il/index_rerun.php?run=$VARS{run_number}&file=$Seq_File_without_low_SP_SEQ_with_Names_NO_PATH')}\"><br>");}
			else {print_message_to_output("<A HREF='$Seq_File_without_low_SP_SEQ_with_Names_NO_PATH' TARGET=_blank>The input sequences after the removal of unreliable sequences (with confidence score below $Cutoff)</A><font size=-1> (see list of removed sequences <A HREF='$removed_low_SP_SEQ_With_Names_NO_PATH' TARGET=_blank>here</A></font>)&nbsp;&nbsp;&nbsp;<INPUT TYPE=\"BUTTON\" VALUE=\"run GUIDANCE on the confidently-aligned sequences only\" ONCLICK=\"window.open('http://guidance.tau.ac.il/index_rerun.php?run=$VARS{run_number}&file=$Seq_File_without_low_SP_SEQ_with_Names_NO_PATH')\"><br>");}
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
