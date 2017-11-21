#!/usr/bin/perl -w

use strict;
use Storable;

#use lib "/bioseq/bioSequence_scripts_and_constants";
#use lib "/bioseq/pupkoSVN/trunk/www/bioSequence_scripts_and_constants/";
#use lib "/bioseq/Guidance/"; # == /bioseq/pupkoSVN/trunk/www/Guidance
#use lib "/bioseq/Selecton/external_scripts";
use FindBin qw($Bin);

use lib $Bin; # == pupkoSVN/www/Guidance
use lib "$Bin/../bioSequence_scripts_and_constants/";
use lib "$Bin/../Selecton/";

use File::Basename;
use codonAlign;

use MSA_parser;
use Guidance;
use GENERAL_CONSTANTS;
use BIOSEQUENCE_FUNCTIONS;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Tools::CodonTable;
use Getopt::Long;
use Bio::Seq;
use File::Copy;
use List::Util qw (shuffle);
# Determine input mode - commandline or in server
my $isServer;
my (%VARS, %FORM);
my $proc_num;
my $outDir;
my $stored_data_file;
my $stored_form_data;
my ($semphy_prog,$mafft_prog, $prank_prog,$clustalw_prog,$ruby_prog,$muscle_prog,$msa_set_score_prog,$pagan_prog,$remove_taxa_prog);
my @ORIG_ARGV=@ARGV;
my $overview_URL="";
my $gallery_URL="";
my $home_URL="";
my $status_file=""; # will follow the status of alternative MSA creation on the web-server
die "USAGE: --seqFile <seqFile> --msaProgram <MAFFT|PRANK|CLUSTALW|MUSCLE|PAGAN> --seqType <aa|nuc|codon> --outDir <full path outDir>
Optional parameters:
  --program <GUIDANCE|HoT|GUIDANCE2> default=GUIDANCE2
  --bootstraps <number of bootstrap iterations> default=100
  --genCode <option value> default=1
                <option value=1>  Nuclear Standard
                <option value=15> Nuclear Blepharisma
           	    <option value=6>  Nuclear Ciliate
                <option value=10> Nuclear Euplotid
                <option value=2>  Mitochondria Vertebrate
                <option value=5>  Mitochondria Invertebrate
                <option value=3>  Mitochondria Yeast
                <option value=13> Mitochondria Ascidian
                <option value=9>  Mitochondria Echinoderm
                <option value=14> Mitochondria Flatworm
                <option value=4>  Mitochondria Protozoan
  --outOrder <aligned|as_input> default=aligned
  --msaFile <msaFile> - not recommended, see documentation online guidance.tau.ac.il
  --seqCutoff <confidence cutoff between 0 to 1> default=0.6
  --colCutoff <confidence cutoff between 0 to 1> default=0.93
  --Z_Seq_Cutoff <Z score as additional criteria to filter sequences> EXPERIMENTAL, default=NA (not active)
  --Z_Col_Cutoff <Z score as additional criteria to filter position> EXPERIMENTAL, default=NA (not active)
  --mafft <path to mafft executable> default=mafft
  --prank <path to prank executable> default=prank
  --clustalw <path to clustalw executable> default=clustalw
  --muscle <path to muscle executable> default=muscle
  --pagan <path to pagan executable> default=pagan
  --ruby <path to ruby executable> default=ruby
  --dataset Unique name for the Dataset - will be used as prefix to outputs (default=MSA)
  --MSA_Param passing parameters for the alignment program. To pass parameter containning '-' in it, add \\ before each '-' e.g. \\-F for PRANK
  --proc_num <num of processors to use> default=1 
" unless (@ARGV >= 1);

if ($ARGV[0] =~ m/^-/) {
	# Commandline input mode
	$isServer = 0;
	$VARS{proc_num} = 1; 
	$semphy_prog = "$Bin/../../programs/semphy/semphy";
	$msa_set_score_prog="$Bin/../../programs/msa_set_score/msa_set_score";
	$mafft_prog = "mafft";
	$prank_prog = "prank";
	$clustalw_prog="clustalw";
	$muscle_prog="muscle";
	$pagan_prog="pagan";
	$ruby_prog = "ruby";
	$remove_taxa_prog= "$Bin/../../programs/removeTaxa/removeTaxa";
	$VARS{align_param} = "";
	$FORM{userMSA_File} = "";
	$FORM{PROGRAM} = "GUIDANCE2";
	$FORM{Bootstraps} = 100;
	$FORM{CodonTable} = 1;  # Nuclear standard code
	$FORM{Align_Order} = "aligned";
	$VARS{SP_SEQ_CUTOFF} = "0.6";
	$VARS{SP_COL_CUTOFF} = "0.93";
	$VARS{Z_Col_Cutoff} = "NA";
	$VARS{Z_Seq_Cutoff} = "NA";
	$VARS{dataset}="MSA";
	$FORM{usrSeq_File}="";
	$FORM{MSA_Program}="";
	$VARS{rooting_type}="BioPerl";
	$VARS{BBL}="NO";
	$VARS{GapPenDist}="";
	my $getoptResult = GetOptions ("seqFile=s"=>\$FORM{usrSeq_File},  # = means that this parameter is required, s means string
								   "outDir=s"=>\$outDir,
 								   "msaFile:s"=>\$FORM{userMSA_File},  # : means that this parameter is optional
 								   "program:s"=>\$FORM{PROGRAM},       # <GUIDANCE|HoT|GUIDANCE2|GUIDANCE3|GUIDANCE3_HOT> default=GUIDANCE2 [GUIDANCE3=pertube ep and op without HOT | GUIDANCE3_HOT=pertube ep and op with HOT]
 								   "msaProgram=s"=>\$FORM{MSA_Program},
 								   "bootstraps:i"=>\$FORM{Bootstraps},
 								   "genCode:i"=>\$FORM{CodonTable},
 								   "seqType=s"=>\$FORM{Seq_Type},
 								   "outOrder:s"=>\$FORM{Align_Order}, # will be as_input or aligned
								   "seqCutoff:s"=>\$VARS{SP_SEQ_CUTOFF},
								   "colCutoff:s"=>\$VARS{SP_COL_CUTOFF},
								   "Z_Seq_Cutoff:s"=>\$VARS{Z_Seq_Cutoff}, # <Z score as additional criteria to filter sequences> EXPERIMENTAL, default=NA (not active)
								   "Z_Col_Cutoff:s"=>\$VARS{Z_Col_Cutoff}, # <Z score as additional criteria to filter position> EXPERIMENTAL, default=NA (not active)
								   "mafft:s"=>\$mafft_prog,
								   "prank:s"=>\$prank_prog,
								   "clustalw:s"=>\$clustalw_prog,
								   "muscle:s"=>\$muscle_prog,
								   "pagan:s"=>\$pagan_prog,
								   "ruby:s"=>\$ruby_prog,
								   "dataset:s"=>\$VARS{dataset},
								   "MSA_Param:s"=>\$VARS{align_param},
								   "proc_num:i"=>\$VARS{proc_num},
								   ### EXPERIMENTAL FEATURES.... MOST ACTIVE ONLY LOCAL
								   "RootingType:s"=>\$VARS{rooting_type}, # BioPerl (default) | MidPoint
								   "BBL:s"=>\$VARS{BBL}, # Do branch length optimization (BBL) - YES | NO - ONLY ACTIVE LOCAL
								   "GapPenDist:s"=>\$VARS{GapPenDist}, # UNIF | EMP - sample gap panalties from unif[orm] or default: emp[erical] distribution => RELEVANT ONLY FOR GUIDANCE 2
		);

	if (uc ($VARS{BBL}) eq "YES"){$semphy_prog="/groups/pupko/haim/Programs/semphy_test_clean_log_BBL/programs/semphy/semphy";}
	
	die "ERROR: Bootstraps parameter must be a number\n" if ($FORM{Bootstraps} !~/^\d+$/);
	die "ERROR: No path for output\n" if ($outDir eq "");
	die "ERROR: seqFile is requiered\n" if ($FORM{usrSeq_File} eq "");
	die "ERROR: number of processors must be >= 1\n" if ($VARS{proc_num} < 1);  
	print "WARNING: PAGAN require that MAFFT is also installed on your system. Otherwise GUIDANCE will fail\n" if ($FORM{MSA_Program} eq "PAGAN");
	$FORM{MSA_Program}=uc($FORM{MSA_Program});
	if (($FORM{MSA_Program} ne "MAFFT") and ($FORM{MSA_Program} ne "PRANK") and ($FORM{MSA_Program} ne "CLUSTALW") and ($FORM{MSA_Program} ne "MUSCLE") and ($FORM{MSA_Program} ne "PAGAN"))
	{
		die "ERROR: msaProgram should be MAFFT or PRANK or CLUSTALW or MUSCLE or PAGAN (case sensative)\n";
	}
	unless ($outDir =~ m/\/$/) {
		$outDir .= "/";
	}
	print "outDir: $outDir\n";
	$VARS{WorkingDir} = $outDir;
	unless (-e $outDir) {
		system("mkdir $outDir");
	}
	if (($VARS{GapPenDist} eq "") and ($FORM{PROGRAM} eq "GUIDANCE3"))
	{
		$VARS{GapPenDist}="UNIF";
	}
	elsif (($VARS{GapPenDist} eq "") and ($FORM{PROGRAM} eq "GUIDANCE3_HOT"))
	{
		$VARS{GapPenDist}="UNIF";
	}
	elsif (($VARS{GapPenDist} eq "") and ($FORM{PROGRAM} eq "GUIDANCE2"))
	{
		# $VARS{GapPenDist}="EMP";
		$VARS{GapPenDist}="UNIF";
	}
	
	if ($FORM{PROGRAM} eq "GUIDANCE2")
	{
		$overview_URL="http://guidance.tau.ac.il/ver2/overview.php";
		$gallery_URL="http://guidance.tau.ac.il/ver2/Gallery.php";
		$home_URL="http://guidance.tau.ac.il/ver2/";
	}
	else
	{
		$overview_URL="http://guidance.tau.ac.il/overview.html";
		$gallery_URL="http://guidance.tau.ac.il/Gallery.htm";
		$home_URL="http://guidance.tau.ac.il/";
	}
	$VARS{OutLogFile} = "$outDir/log";
	open (LOG,">>$VARS{OutLogFile}") || exit_on_error('sys_error', "Can't open Log File: $VARS{OutLogFile} $!");
	print LOG "\n\n========================================= NEW GUIDANCE RUN STARTED ===========================================\n";
	print LOG "GUIDANCE COMMAND: perl $0 ".join (" ",@ORIG_ARGV)."\n";
	# INPUT FILES
	$VARS{SeqsFile} = "Seqs.Orig.fas"; #generic fixed name for the sequence file. the user file will be copied to this file.
	$VARS{SeqsFile_Codons} = "Seqs.Orig_DNA.fas"; #generic fixed name for the DNA CODONS sequence file. the user file will be copied to this file.
    if ($FORM{usrSeq_File} ne ""){
		system ("cp $FORM{usrSeq_File} $VARS{WorkingDir}$VARS{SeqsFile}");
		BIOSEQUENCE_FUNCTIONS::convertNewline($VARS{SeqsFile},$VARS{WorkingDir}) if ($isServer==1);
		BIOSEQUENCE_FUNCTIONS::removeEndLineExtraChars($VARS{SeqsFile},$VARS{WorkingDir});
	}

	if ($FORM{Seq_Type} eq "aa")
	{
		$FORM{Seq_Type} = "AminoAcids";
	}
	elsif ($FORM{Seq_Type} eq "nuc")
	{
		$FORM{Seq_Type} = "Nucleotides" if ($FORM{Seq_Type} eq "nuc");
	}
    elsif ($FORM{Seq_Type} eq "codon") 
	{
		$FORM{Seq_Type} = "Codons";
        system ("mv $VARS{WorkingDir}$VARS{SeqsFile} $VARS{WorkingDir}$VARS{SeqsFile_Codons}");
	}
	else
	{
		die "ERROR: --seqType must be:aa or nuc or codon\n";
	}
	if (($VARS{align_param}=~/addfragments (\S+)/) and ($FORM{MSA_Program} eq "MAFFT")) # support --addfragments option of MAFFT, seqFile is the seq for the coreMSA
	{
		my $user_fragments_file=$1;
		$VARS{fragments_file_name}="add_fragments_file";
		$VARS{fragments_file_name_seqName_coded}="$VARS{fragments_file_name}.SeqNameCoded";
		copy($user_fragments_file,"$VARS{WorkingDir}$VARS{fragments_file_name}") or exit_on_error ("sys_error","Can't copy user fragments file: $user_fragments_file to $VARS{WorkingDir}$VARS{fragments_file_name}\n");
		print LOG "Validating fragments: Guidance::validate_Seqs($VARS{WorkingDir},$VARS{fragments_file_name},$FORM{Seq_Type},No): ";
		my @ans=Guidance::validate_Seqs($VARS{WorkingDir},$VARS{fragments_file_name},$FORM{Seq_Type},"No"); # validate fragmets

		print LOG "return: ","@ans","\n";
		if ($ans[0] eq "sys_error") {exit_on_error ('sys_error',$ans[1]);}
		elsif ($ans[0] ne "OK"){exit_on_error('user_error', join ("",@ans));}
		if (($ans[0] eq "OK") and ($ans[1] ne ""))
		{
			print LOG "Warnning: $ans[1]; Nevertheless calculation is continued\n";
			print "Warnning: $ans[1]; Nevertheless calculation is continued\n";
#			$VARS{SeqsFile}=$ans[2];
		}
		close (LOG);
		$VARS{fragments_file_name}=$ans[2];
		$VARS{NumOfFragments}=$ans[3];

		my @tmp=split (/\Q\-\-\E/,$VARS{align_param}); # \Q and \E are quotemeta 
		#print "JOIN:",join("\\\-\\\-",@tmp);<STDIN>;
		#print "JOIN:",join("*",@tmp);<STDIN>;
		my $tmp_size=@tmp;
		if ($tmp_size>=1) # command line type
		{
			foreach my $param (@tmp)
			{
				#print "PARAM BEFORE:$param\n";
				if ($param=~/addfragments/){$param="addfragments $VARS{WorkingDir}$VARS{fragments_file_name_seqName_coded} ";}
				#print "PARAM AFTER:$param\n";
			}
			#print "TMP:@tmp";
			$VARS{align_param}=join ("\\\-\\\-",@tmp);
			#print "ALN_PARAM L189:$VARS{align_param}\n";<STDIN>;
		}
		else # server type
		{
			my @tmp=split (/\Q--\E/,$VARS{align_param});
			foreach my $param (@tmp)
			{
				if ($param=~/addfragments/){$param="addfragments $VARS{WorkingDir}$VARS{fragments_file_name_seqName_coded}";}
			}
			$VARS{align_param}=join ("\\\-\\\-",@tmp);
		}
		#$VARS{align_param}=~s/addfragments $user_fragments_file/addfragments $VARS{WorkingDir}$VARS{fragments_file_name_seqName_coded}/;
		#print "ALN_PARAM L202:$VARS{align_param}\n";<STDIN>;
		# TO DO CHECK IF NEED TO REMOVE REORDER WITH FRAGMENTS
		if ($VARS{align_param}=~/\-\-reorder/)
		{
			$VARS{align_param}=~s/\-\-reorder//; # if seed is provided reorder must be removed so the seeds will be first
			print "WARNNING: --reorder is not allowed if seed alignment is provided, therefore the --reorder argument will be ignored and the output order will be same as input (with seeds first)\n";
			print LOG "WARNNING: --reorder is not allowed if seed alignment is provided, therefore the --reorder argument will be ignored and the output order will be same as input (with seeds first)\n";
		}
	}
	
	if (($VARS{align_param}=~/\-\-reorder/) and ($VARS{align_param}=~/seed/))
	{
		$VARS{align_param}=~s/\-\-reorder//; # if seed is provided reorder must be removed so the seeds will be first
		print "WARNNING: --reorder is not allowed if seed alignment is provided, therefore the --reorder argument will be ignored and the output order will be same as input (with seeds first)\n";
		print LOG "WARNNING: --reorder is not allowed if seed alignment is provided, therefore the --reorder argument will be ignored and the output order will be same as input (with seeds first)\n";
	}
	if ($VARS{align_param}=~/\-\-retree ([0-9]+)/)
	{
		if ($1>1)
		{
			print "WARNNING: --retree $1 is not supported in GUIDANCE, therefore this argument is ignored.\n";
			print LOG "WARNNING: --retree $1 is not supported in GUIDANCE, therefore this argument is ignored.\n";
			$VARS{align_param}=~s/\-\-retree ([0-9]+)//;
		}
	}

	close (LOG);

	#### OutPut Files
#	$VARS{dataset}="MSA";
	if ($FORM{userMSA_File} ne "") { # The user gave the base MSA as input
		$VARS{Alignment_File}="UserMSA";
		if (!-e $VARS{WorkingDir}.$VARS{Alignment_File}) {
			system ("cp $FORM{userMSA_File} $VARS{WorkingDir}$VARS{Alignment_File}");
			BIOSEQUENCE_FUNCTIONS::convertNewline($VARS{Alignment_File},$VARS{WorkingDir}) if ($isServer==1);
			BIOSEQUENCE_FUNCTIONS::removeEndLineExtraChars($VARS{Alignment_File},$VARS{WorkingDir});
		}
	}
	else {
		$VARS{Alignment_File}="$VARS{dataset}.".$FORM{MSA_Program}.".aln";
	}
	if ($VARS{align_param}=~/addfragments/) {
		$VARS{Core_Alignment_File}="$VARS{dataset}.".$FORM{MSA_Program}.".CORE.aln";
	}
	$VARS{output_page} = "log";
#	if ($FORM{Seq_Type} eq "Codons") {
#		$VARS{output_page} = "log";
#	}
	
	# Validate Seqs NEW
	if ($FORM{userMSA_File} eq "") # Seq file provided
	{
		open (LOG,">>$VARS{OutLogFile}") || exit_on_error('sys_error', "Can't open Log File: $VARS{OutLogFile} $!");
		my @ans;
		if ($FORM{Seq_Type} ne "Codons")
		{
			print LOG "Guidance::validate_Seqs($VARS{WorkingDir},$VARS{SeqsFile},$FORM{Seq_Type},No): ";
			@ans=Guidance::validate_Seqs($VARS{WorkingDir},$VARS{SeqsFile},$FORM{Seq_Type},"No");
		}
		else
		{
			print LOG "Guidance::validate_Seqs($VARS{WorkingDir}$VARS{SeqsFile_Codons},$FORM{Seq_Type},No): ";
			@ans=Guidance::validate_Seqs($VARS{WorkingDir},$VARS{SeqsFile_Codons},$FORM{Seq_Type},"No");
		}
		print LOG "return: ","@ans","\n";
		if ($ans[0] eq "sys_error") {exit_on_error ('sys_error',$ans[1]);}
		elsif ($ans[0] ne "OK"){exit_on_error('user_error', join ("",@ans));}
		if (($ans[0] eq "OK") and ($ans[1] ne ""))
		{
			print LOG "Warnning: $ans[1]; Nevertheless calculation is continued\n";
			print "Warnning: $ans[1]; Nevertheless calculation is continued\n";
#			$VARS{SeqsFile}=$ans[2];
		}
		close (LOG);
		$VARS{SeqsFile}=$ans[2];
		$VARS{NumOfSeq}=$ans[3];
	}
	elsif ((-e "$VARS{WorkingDir}$VARS{Alignment_File}") and (-s "$VARS{WorkingDir}$VARS{Alignment_File}" > 0)) # Alignment provided
	{
		open (LOG,">>$VARS{OutLogFile}") || exit_on_error('sys_error', "Can't open Log File: $VARS{OutLogFile} $!");
		my @ans=();
		if ($FORM{Seq_Type} ne "Codons")
		{
			print LOG "Guidance::validate_Seqs($VARS{WorkingDir},$VARS{Alignment_File},$FORM{Seq_Type},Yes): ";
			@ans=Guidance::validate_Seqs($VARS{WorkingDir},$VARS{Alignment_File},$FORM{Seq_Type},"Yes");
		}
		else
		{
			print LOG "Guidance::validate_Seqs($VARS{WorkingDir},$VARS{Alignment_File},$FORM{Seq_Type},Yes,$FORM{CodonTable}) ";
			@ans=Guidance::validate_Seqs($VARS{WorkingDir},$VARS{Alignment_File},$FORM{Seq_Type},"Yes",$FORM{CodonTable});
		}
		print LOG "return: ","@ans","\n";
		if ($ans[0] eq "sys_error") {exit_on_error ('sys_error',$ans[1]);}
		elsif ($ans[0] ne "OK"){exit_on_error('user_error', join ("",@ans));}
		if (($ans[0] eq "OK") and ($ans[1] ne ""))
		{
			print "Warnning: $ans[1]; Nevertheless calculation is continued\n";
			print LOG "Warnning: $ans[1]; Nevertheless calculation is continued\n";
			#$VARS{SeqsFile}=$ans[2];
			#$VARS{Alignment_File}=$ans[2];
		}
		close (LOG);
		$VARS{Alignment_File}=$ans[2];
		$VARS{NumOfSeq}=$ans[3];
	}
	if (($VARS{NumOfSeq}<4) and ($FORM{PROGRAM} eq "GUIDANCE"))
	{
		exit_on_error('user_error', "Only $VARS{NumOfSeq} sequences were provided, however at least 4 sequences are requiered for GUIDANCE<br>You can run HoT algorithm instead.");
	}
	if (($VARS{NumOfSeq}<4) and ($FORM{PROGRAM} eq "GUIDANCE2"))
	{
		exit_on_error('user_error', "Only $VARS{NumOfSeq} sequences were provided, however at least 4 sequences are requiered for GUIDANCE2<br>You can run HoT algorithm instead.");
	}
} else {
	# Server input mode
	$isServer = 1;

	$stored_data_file = $ARGV[0];
	$stored_form_data = $ARGV[1];

	my $vars_ref = retrieve($stored_data_file);
	%VARS = %$vars_ref;
	my $form_ref = retrieve($stored_form_data);
	%FORM = %$form_ref;

	$VARS{dataset}="MSA";

	# GUIDANCE2 server deafults
	$VARS{rooting_type}="BioPerl";
	$VARS{BBL}="NO";
	$VARS{GapPenDist}="UNIF";

	if (uc $VARS{BBL} eq "NO")
	{
		$semphy_prog = GENERAL_CONSTANTS::SEMPHY;
	}
	else
	{
		$semphy_prog = GENERAL_CONSTANTS::SEMPHY_BBL; # TO DO: Change it location to more stable one
	}
	$mafft_prog = GENERAL_CONSTANTS::MAFFT_GUIDANCE;
	$prank_prog = GENERAL_CONSTANTS::PRANK_LECS;
	$clustalw_prog = GENERAL_CONSTANTS::CLUSTALW_LECS;
	$muscle_prog = GENERAL_CONSTANTS::MUSCLE_LECS;
	$pagan_prog= GENERAL_CONSTANTS::PAGAN_LECS;
	$ruby_prog = GENERAL_CONSTANTS::RUBY;
	$msa_set_score_prog=Guidance::MSA_SET_SCORE;

	# deafults , still not supported by the web server implementation, experimental feature
	if (!exists $VARS{Z_Col_Cutoff}){$VARS{Z_Col_Cutoff} = "NA";}
	if (!exists  $VARS{Z_Seq_Cutoff}){$VARS{Z_Seq_Cutoff} = "NA";}

	if ($FORM{CALLING_SERVER} eq "GUIDANCE2")
	{
		$overview_URL="http://guidance.tau.ac.il/ver2/overview.php";
		$gallery_URL="http://guidance.tau.ac.il/ver2/Gallery.php";
		$home_URL="http://guidance.tau.ac.il/ver2/";
	}
	else
	{
		$overview_URL="http://guidance.tau.ac.il/overview.html";
		$gallery_URL="http://guidance.tau.ac.il/Gallery.htm";
		$home_URL="http://guidance.tau.ac.il/";
	}
}

# Codons handaling
$VARS{TranslateErrors}="xCodons.html";
# Output Files
$VARS{code_fileName}="Seqs.Codes";
$VARS{codded_seq_fileName}="Seqs.numberd.fas";
$VARS{codded_seq_fileName_Codons}="Seqs_DNA.numberd.fas";
#$VARS{Alignment_File}="$VARS{dataset}.".$FORM{MSA_Program}.".aln";
$VARS{Alignment_File_PROT}="$VARS{dataset}.".$FORM{MSA_Program}.".PROT.aln";
$VARS{BootStrap_Dir}=$VARS{WorkingDir}."BP/";
$VARS{BootStrap_MSA_Dir}=$VARS{BootStrap_Dir}."BP_MSA/";
$VARS{HoT_MSAs_Dir}="COS_MSA";
$VARS{GUIDANCE2_MSAs_Dir}=$VARS{BootStrap_Dir}."GUIDANCE2_MSA/";
$VARS{Tree_File}="";
$VARS{Semphy_OutFile}="";
$VARS{Semphy_LogFile}="";
$VARS{Semphy_StdFile}="";
$VARS{COL_SCORES_FIGURE}="Col_Scores_Graph.png";
$VARS{Scoring_Alignments_Dir}=""; #Tha dir with the alignment used to create the score ($VARS{BootStrap_Dir} in GUIDANCE and $VARS{HoT_MSA_Dir} in Hot
$VARS{send_email_dir}=GENERAL_CONSTANTS::SEND_EMAIL_DIR_IBIS;

my %DNA_AA=();
if ($isServer == 1) {
	open (OUTPUT,">>$VARS{WorkingDir}$VARS{output_page}");
}
open (LOG,">>$VARS{OutLogFile}") || exit_on_error('sys_error', "Can't open Log File: $VARS{OutLogFile} $!");

#USER Provide an MSA
if ((-e "$VARS{WorkingDir}$VARS{Alignment_File}") and (-s "$VARS{WorkingDir}$VARS{Alignment_File}" > 0))
  {
	  $VARS{code_fileName_aln}=$VARS{code_fileName}.".fromALN";
	  if ($FORM{Seq_Type} ne "Codons")
	  {
		  print LOG "==== Alignment_File:$VARS{Alignment_File}\tSeq_Type: $FORM{Seq_Type}\n";
		  print LOG "extract_seq_from_MSA($VARS{WorkingDir}$VARS{Alignment_File},$VARS{WorkingDir}$VARS{SeqsFile})\n";
		  extract_seq_from_MSA("$VARS{WorkingDir}$VARS{Alignment_File}","$VARS{WorkingDir}$VARS{SeqsFile}");
	  }
	  else
	  {
		  print LOG "==== Alignment_File:$VARS{Alignment_File}\tSeq_Type: $FORM{Seq_Type}\n";
		  print LOG "extract_seq_from_MSA($VARS{WorkingDir}$VARS{Alignment_File},$VARS{WorkingDir}$VARS{SeqsFile_Codons})\n";
		  extract_seq_from_MSA("$VARS{WorkingDir}$VARS{Alignment_File}","$VARS{WorkingDir}$VARS{SeqsFile_Codons}");
	  }
	  if (($FORM{PROGRAM} eq "GUIDANCE2")||($FORM{PROGRAM} eq "GUIDANCE3_HOT"))
	  {
		  print LOG "Guidance::name2codeFastaFrom1($VARS{WorkingDir}$VARS{Alignment_File}, $VARS{WorkingDir}$VARS{code_fileName_aln}, $VARS{WorkingDir}$VARS{Alignment_File}.WithCodesName,0,seqNum);\n";
		  my @ans=Guidance::name2codeFastaFrom1("$VARS{WorkingDir}$VARS{Alignment_File}", "$VARS{WorkingDir}$VARS{code_fileName_aln}", "$VARS{WorkingDir}$VARS{Alignment_File}.WithCodesName",0,"seqNum");
		  unless ($ans[0] eq "ok") {exit_on_error("sys_error","Guidance::name2codeFastaFrom1: ".join (" ",@ans));}
	  }
	  else
	  {
		  print LOG "Guidance::name2codeFastaFrom1($VARS{WorkingDir}$VARS{Alignment_File}, $VARS{WorkingDir}$VARS{code_fileName_aln}, $VARS{WorkingDir}$VARS{Alignment_File}.WithCodesName);\n";
		  my @ans=Guidance::name2codeFastaFrom1("$VARS{WorkingDir}$VARS{Alignment_File}", "$VARS{WorkingDir}$VARS{code_fileName_aln}", "$VARS{WorkingDir}$VARS{Alignment_File}.WithCodesName");
		  unless ($ans[0] eq "ok") {exit_on_error("sys_error","Guidance::name2codeFastaFrom1: ".join (" ",@ans));}
	  }
	  if (($FORM{Seq_Type} eq "Nucleotides")) # ?? ($FORM{Seq_Type} eq "Codons") ??
	  {
		  if ($FORM{MSA_Program} eq "MAFFT")
		  {
			  convert_fs_to_lower_case("$VARS{WorkingDir}$VARS{Alignment_File}.WithCodesName"); # MAFFT ALWAYS OUTPUT NUC MSA IN LOWER CASE
		  }
		  elsif ($FORM{MSA_Program} eq "PRANK")
		  {
			  convert_fs_to_upper_case("$VARS{WorkingDir}$VARS{Alignment_File}.WithCodesName"); # PRANK ALWAYS OUTPUT NUC MSA IN UPPER CASE
		  }
		  elsif ($FORM{MSA_Program} eq "CLUSTALW")
		  {
			  convert_fs_to_upper_case("$VARS{WorkingDir}$VARS{Alignment_File}.WithCodesName"); # CLUSTALW ALWAYS OUTPUT NUC MSA IN UPPER CASE
		  }
		  elsif ($FORM{MSA_Program} eq "MUSCLE")
		  {
			  convert_fs_to_upper_case("$VARS{WorkingDir}$VARS{Alignment_File}.WithCodesName"); # MUSCLE ALWAYS OUTPUT NUC MSA IN UPPER CASE
		  }
		  elsif ($FORM{MSA_Program} eq "PAGAN")
		  {
			  convert_fs_to_upper_case("$VARS{WorkingDir}$VARS{Alignment_File}.WithCodesName"); # PAGAN ALWAYS OUTPUT NUC MSA IN UPPER CASE
		  }
	  }
	  if ((-z "$VARS{WorkingDir}$VARS{Alignment_File}.WithCodesName") or !(-e "$VARS{WorkingDir}$VARS{Alignment_File}.WithCodesName")) {exit_on_error("user_error","Sequences were not found on the <A htef=\"$VARS{Alignment_File}\">uploades file</A><br>Make sure it is a Plain text FASTA Format<br>");}
	  if ($FORM{PROGRAM} eq "HoT")
	  {
		  names_according_CoS("$VARS{WorkingDir}$VARS{Alignment_File}.WithCodesName");
      }
	  system ("cp $VARS{WorkingDir}$VARS{Alignment_File} $VARS{WorkingDir}$VARS{Alignment_File}.ORIG");
	  system ("mv $VARS{WorkingDir}$VARS{Alignment_File}.WithCodesName $VARS{WorkingDir}$VARS{Alignment_File}");
	  
  }
# Convert Fasta names to Codes
####################################
my $seq_counter=0;
if ($VARS{align_param}=~/\-\-seed/)
{
	my @seeds=$VARS{align_param}=~m/\-\-seed ($VARS{WorkingDir}seedFile[0-9]+)/g;
	foreach my $seed_file (@seeds)
	{
		print LOG "Guidance::name2codeFasta_without_codded_out($seed_file, $VARS{WorkingDir}$VARS{code_fileName},$seq_counter)\n";
		my @ans=Guidance::name2codeFasta_without_codded_out($seed_file, "$VARS{WorkingDir}$VARS{code_fileName}", $seq_counter);
		if ($ans[0] ne "ok") {exit_on_error("sys_error","Guidance::name2codeFasta_without_codded_out: ".join (" ",@ans));}
		else {$seq_counter=$ans[1];}
	}
}
if ($FORM{Seq_Type} ne "Codons")
{
	my @ans=();
	if (($FORM{PROGRAM} eq "GUIDANCE2")||($FORM{PROGRAM} eq "GUIDANCE3_HOT"))
	{
		print LOG "Guidance::name2codeFastaFrom1($VARS{WorkingDir}$VARS{SeqsFile}, $VARS{WorkingDir}$VARS{code_fileName}, $VARS{WorkingDir}$VARS{codded_seq_fileName},$seq_counter,seqNum);\n";
		@ans=Guidance::name2codeFastaFrom1("$VARS{WorkingDir}$VARS{SeqsFile}", "$VARS{WorkingDir}$VARS{code_fileName}", "$VARS{WorkingDir}$VARS{codded_seq_fileName}",$seq_counter,"seqNum");
	}
	else
	{
		print LOG "Guidance::name2codeFastaFrom1($VARS{WorkingDir}$VARS{SeqsFile}, $VARS{WorkingDir}$VARS{code_fileName}, $VARS{WorkingDir}$VARS{codded_seq_fileName},$seq_counter);\n";
		@ans=Guidance::name2codeFastaFrom1("$VARS{WorkingDir}$VARS{SeqsFile}", "$VARS{WorkingDir}$VARS{code_fileName}", "$VARS{WorkingDir}$VARS{codded_seq_fileName}",$seq_counter);
	}
	unless ($ans[0] eq "ok") {exit_on_error("sys_error","Guidance::name2codeFastaFrom1: ".join (" ",@ans));}
	if ((-z "$VARS{WorkingDir}$VARS{codded_seq_fileName}") or !(-e "$VARS{WorkingDir}$VARS{codded_seq_fileName}")) {exit_on_error("user_error","Sequences were not found on the <A htef=\"$VARS{SeqsFile}\">uploades file</A><br>Make sure it is a Plain text FASTA Format<br>");}

	if (($VARS{align_param}=~/addfragments/) and ($FORM{MSA_Program} eq "MAFFT")) # support --addfragments option of MAFFT, seqFile is the seq for the coreMSA
	{
		my $offset=$VARS{NumOfSeq}+1;
		if (($FORM{PROGRAM} eq "GUIDANCE2") || ($FORM{PROGRAM} eq "GUIDANCE3_HOT"))
		{
			pprint LOG "Guidance::name2codeFastaFrom1($VARS{WorkingDir}$VARS{fragments_file_name},$VARS{WorkingDir}$VARS{code_fileName},$VARS{WorkingDir}$VARS{fragments_file_name_seqName_coded},$offset,seqNum);\n";
			@ans=Guidance::name2codeFastaFrom1("$VARS{WorkingDir}$VARS{fragments_file_name}", "$VARS{WorkingDir}$VARS{code_fileName}","$VARS{WorkingDir}$VARS{fragments_file_name_seqName_coded}",$offset,"seqNum");
		}
		else
		{
			print LOG "Guidance::name2codeFastaFrom1($VARS{WorkingDir}$VARS{fragments_file_name},$VARS{WorkingDir}$VARS{code_fileName},$VARS{WorkingDir}$VARS{fragments_file_name_seqName_coded},$offset);\n";
			@ans=Guidance::name2codeFastaFrom1("$VARS{WorkingDir}$VARS{fragments_file_name}", "$VARS{WorkingDir}$VARS{code_fileName}","$VARS{WorkingDir}$VARS{fragments_file_name_seqName_coded}",$offset);
		}
		unless ($ans[0] eq "ok") {exit_on_error("sys_error","Guidance::name2codeFastaFrom1: ".join (" ",@ans));}
		if ((-z "$VARS{WorkingDir}$VARS{fragments_file_name_seqName_coded}") or !(-e "$VARS{WorkingDir}$VARS{fragments_file_name_seqName_coded}")) {exit_on_error("user_error","Sequences were not found on the <A htef=\"$VARS{fragments_file_name}\">uploades file</A><br>Make sure it is a Plain text FASTA Format<br>");}
		# create file with the list of fragments name
		$VARS{fragments_codes}="add_fragments.seqCodes";
		open (FRAGMENTS_NAMES,">$VARS{WorkingDir}$VARS{fragments_codes}") || exit_on_error ("sys_error", "Can't open fragments codes file name: '$VARS{WorkingDir}$VARS{fragments_codes}' $!");
		for (my $i=($VARS{NumOfSeq}+1);$i<=($VARS{NumOfSeq}+$VARS{NumOfFragments});$i++)
		{
			print FRAGMENTS_NAMES "$i\n";
		}
		close FRAGMENTS_NAMES;
	}

}
else
{
	my @ans=();
	if (($FORM{PROGRAM} eq "GUIDANCE2") or ($FORM{PROGRAM} eq "GUIDANCE3") or ($FORM{PROGRAM} eq "GUIDANCE3_HOT")) # TO DO: check if conversion is needed in GUIDANCE3 without HoT, did not do it for non codons seq
	{
		print LOG "Guidance::name2codeFastaFrom1($VARS{WorkingDir}$VARS{SeqsFile_Codons}, $VARS{WorkingDir}$VARS{code_fileName}, $VARS{WorkingDir}$VARS{codded_seq_fileName},$seq_counter,seqNum);\n";
		@ans=Guidance::name2codeFastaFrom1("$VARS{WorkingDir}$VARS{SeqsFile_Codons}", "$VARS{WorkingDir}$VARS{code_fileName}", "$VARS{WorkingDir}$VARS{codded_seq_fileName_Codons}",$seq_counter,"seqNum");
	}
	#elsif ($FORM{PROGRAM} ne "GUIDANCE2")
	else
	{
		print LOG "Guidance::name2codeFastaFrom1($VARS{WorkingDir}$VARS{SeqsFile_Codons}, $VARS{WorkingDir}$VARS{code_fileName}, $VARS{WorkingDir}$VARS{codded_seq_fileName},$seq_counter);\n";
		@ans=Guidance::name2codeFastaFrom1("$VARS{WorkingDir}$VARS{SeqsFile_Codons}", "$VARS{WorkingDir}$VARS{code_fileName}", "$VARS{WorkingDir}$VARS{codded_seq_fileName_Codons}",$seq_counter);
	}
	unless ($ans[0] eq "ok") {exit_on_error("sys_error","Guidance::name2codeFastaFrom1: ".join (" ",@ans));}
	if ((-z "$VARS{WorkingDir}$VARS{codded_seq_fileName_Codons}") or !(-e "$VARS{WorkingDir}$VARS{codded_seq_fileName_Codons}")) {exit_on_error("user_error","Sequences were not found on the <A htef=\"$VARS{SeqsFile_Codons}\">uploades file</A><br>Make sure it is a Plain text in FASTA Format<br>");}
	
	### TO DO: Take care of codons and addfragments
}
close (OUTPUT);
# IF CODONS TRANSLATE TO AA
if ($FORM{Seq_Type} eq "Codons")
{
	my @ans=();
	if (($FORM{PROGRAM} eq "GUIDANCE2") or ($FORM{PROGRAM} eq "GUIDANCE3") or ($FORM{PROGRAM} eq "GUIDANCE3_HOT")) # TO DO: check if conversion is needed in GUIDANCE3 without HoT, did not do it for non codons seq
	{
		print LOG "Translating DNA to AA: codonAlign::translate_DNA_to_AA(\"$VARS{WorkingDir}$VARS{codded_seq_fileName_Codons}\",\"$VARS{WorkingDir}$VARS{codded_seq_fileName}\",$FORM{CodonTable},\"$VARS{WorkingDir}$VARS{TranslateErrors}\",\"$VARS{WorkingDir}$VARS{output_page}\",\%DNA_AA,\"\",\"\",seqNum);\n";
		@ans=codonAlign::translate_DNA_to_AA("$VARS{WorkingDir}$VARS{codded_seq_fileName_Codons}","$VARS{WorkingDir}$VARS{codded_seq_fileName}",$FORM{CodonTable},"$VARS{WorkingDir}$VARS{TranslateErrors}","$VARS{WorkingDir}$VARS{output_page}",\%DNA_AA,"","","seqNum");
	}
	else
	{
		print LOG "Translating DNA to AA: codonAlign::translate_DNA_to_AA(\"$VARS{WorkingDir}$VARS{codded_seq_fileName_Codons}\",\"$VARS{WorkingDir}$VARS{codded_seq_fileName}\",$FORM{CodonTable},\"$VARS{WorkingDir}$VARS{TranslateErrors}\",\"$VARS{WorkingDir}$VARS{output_page}\",\%DNA_AA);\n";
		@ans=codonAlign::translate_DNA_to_AA("$VARS{WorkingDir}$VARS{codded_seq_fileName_Codons}","$VARS{WorkingDir}$VARS{codded_seq_fileName}",$FORM{CodonTable},"$VARS{WorkingDir}$VARS{TranslateErrors}","$VARS{WorkingDir}$VARS{output_page}",\%DNA_AA);
	}
	if ($ans[0] ne "ok")
	{
		if ($ans[0] eq "user")
		{
			exit_on_error("user_error",join(" ",$ans[1]));			
		}
		elsif ($ans[0] eq "sys")
		{
			exit_on_error("sys_error",join(" ",@ans));
		}
	}
}
# NOW WE ALWAYS WITH AA SEQ
if ($isServer ==1) {
	open (OUTPUT,">>$VARS{WorkingDir}$VARS{output_page}");
	print OUTPUT "<h4><font face=Verdana><u>Running Messages:</u></h4></font>\n";
}

if (($FORM{PROGRAM} eq "GUIDANCE")||($FORM{PROGRAM} eq "GUIDANCE3"))
{
	run_Guidance()
}
elsif ($FORM{PROGRAM} eq "HoT")
{
	run_HoT();
}
elsif (($FORM{PROGRAM} eq "GUIDANCE2")||($FORM{PROGRAM} eq "GUIDANCE3_HOT"))
{
	run_Guidance2();
}

# run Giddy's program for calculating SP scores
#################################################
if ($isServer ==1) {
	if ($FORM{PROGRAM} eq "GUIDANCE")
	{
		print_message_to_output("Calculating GUIDANCE scores");
	}	
	if ($FORM{PROGRAM} eq "HoT")
	{
		print_message_to_output("Calculating HoT scores");
	}
	if ($FORM{PROGRAM} eq "GUIDANCE2")
	{
		print_message_to_output("Calculating GUIDANCE2 scores");
	}
	if ($FORM{PROGRAM} eq "GUIDANCE3")
	{
		print_message_to_output("Calculating GUIDANCE3 scores");
	}
}

if (($FORM{PROGRAM} eq "GUIDANCE") || ($FORM{PROGRAM} eq "HoT"))
{
    $VARS{Output_Prefix}=$VARS{dataset}.".$FORM{MSA_Program}."."Guidance";
}
elsif ($FORM{PROGRAM} eq "GUIDANCE2")
{
    $VARS{Output_Prefix}=$VARS{dataset}.".$FORM{MSA_Program}."."Guidance2";
}
my $cmd="";
if (($FORM{userMSA_File} ne "") and ($FORM{Seq_Type} eq "Codons")) # user gave codon alignment
{
	# the user gave codon alignment but we are all working with AA sequences so we convert the user codon alignment to AA alignment
	$VARS{Alignment_File_translated_from_user_codon_alignmet}=$VARS{Alignment_File}.".TranslatedProt";
	print LOG "Codon_Aln_to_AA_Aln($VARS{WorkingDir}$VARS{Alignment_File},$VARS{WorkingDir}$VARS{Alignment_File_translated_from_user_codon_alignmet},$FORM{CodonTable},XCodonsFromALN.html)\n";
	my @ans=Codon_Aln_to_AA_Aln("$VARS{WorkingDir}$VARS{Alignment_File}","$VARS{WorkingDir}$VARS{Alignment_File_translated_from_user_codon_alignmet}",$FORM{CodonTable},"XCodonsFromALN.html");
	if ($ans[0] ne "OK") {exit_on_error('user_error',$ans[0]);} # error
	elsif ($ans[1] ne "") #warnning
	{
		if ($isServer == 1) # server
		{
			print OUTPUT "<br><b><font color=\"red\" size4=>Warnning:</b></font><font size=\"4\"> $ans[1]</font>";
			print LOG "Warnning: $ans[1]\n";
		} 
		else
		{
			print "Warnning: $ans[1]\n";
			print LOG "Warnning: $ans[1]\n";
		}
	}
	if ((-z "$VARS{WorkingDir}/$VARS{Alignment_File_translated_from_user_codon_alignmet}") or (!-e "$VARS{WorkingDir}/$VARS{Alignment_File_translated_from_user_codon_alignmet}")){exit_on_error("sys_error","$VARS{WorkingDir}$VARS{Alignment_File_translated_from_user_codon_alignmet} does not exist/empty\n");}
	$cmd=$msa_set_score_prog." $VARS{WorkingDir}/$VARS{Alignment_File_translated_from_user_codon_alignmet} $VARS{WorkingDir}/$VARS{Output_Prefix} -d  $VARS{Scoring_Alignments_Dir} > $VARS{WorkingDir}/$VARS{dataset}.$FORM{MSA_Program}.msa_set_score.std";
}
else
{
	$cmd=$msa_set_score_prog." $VARS{WorkingDir}/$VARS{Alignment_File} $VARS{WorkingDir}/$VARS{Output_Prefix} -d  $VARS{Scoring_Alignments_Dir} > $VARS{WorkingDir}/$VARS{dataset}.$FORM{MSA_Program}.msa_set_score.std";# 2>\&1";
}
print LOG "calculating SP scores: $cmd\n";
system ($cmd);
if (!-e "$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_res.scr" or -z "$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_res.scr")
{
	exit_on_error("sys_error","$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_res.scr does not exist/empty\n");
}
if ($FORM{PROGRAM} eq "HoT") {
	open (ORIG_ALIGN,"$VARS{WorkingDir}$VARS{Alignment_File}") || die "Can't open $VARS{WorkingDir}$VARS{Alignment_File} $!";
	open (NEW_ALIGN,">$VARS{WorkingDir}$VARS{Alignment_File}.NEW") || die "Can't open $VARS{WorkingDir}$VARS{Alignment_File}.NEW $!";
	while (my $line=<ORIG_ALIGN>)
	{
		if ($line=~/>seq/)
		{
			if ($line=~/>seq[0]+([1-9]+[0-9]*)$/){print NEW_ALIGN ">".($1+1)."\n";}
			else {print NEW_ALIGN ">1\n";}
		}
		else
		{
			print NEW_ALIGN $line;
		}	
	}
	close (NEW_ALIGN);
	close (ORIG_ALIGN);
	system ("mv $VARS{WorkingDir}$VARS{Alignment_File} $VARS{WorkingDir}$VARS{Alignment_File}.ORIG");
	system ("mv $VARS{WorkingDir}$VARS{Alignment_File}.NEW $VARS{WorkingDir}$VARS{Alignment_File}");
}
if ($FORM{Seq_Type} eq "Codons")
{
	# We should modify the Scores files to be for CODONS - i.e each col score is repeated 3 times for the coland col+1,col+2
	
	# the alignment file should be back to CODONS alignment
	if ($FORM{userMSA_File} eq "") # ONLY if User not provided CodonAlignment
	{
		system ("cp $VARS{WorkingDir}$VARS{Alignment_File} $VARS{WorkingDir}$VARS{Alignment_File_PROT}");
		print LOG "cp $VARS{WorkingDir}$VARS{Alignment_File} $VARS{WorkingDir}$VARS{Alignment_File_PROT}";
		# copy the hash - AA_to_DNA_aligned change the hash so we copy it for later usage
		my %DNA_AA_cp=();
		foreach my $key (keys %DNA_AA)
		{
			my @Arr;
			$Arr[0]=$DNA_AA{$key}->[0];
			$Arr[1]=$DNA_AA{$key}->[1];
			$DNA_AA_cp{$key}=[@Arr];
		}
		print LOG "codonAlign::AA_to_DNA_aligned($VARS{WorkingDir}$VARS{Alignment_File_PROT},$VARS{WorkingDir}$VARS{Alignment_File},\%DNA_AA_cp);\n";
		my @ans=codonAlign::AA_to_DNA_aligned("$VARS{WorkingDir}$VARS{Alignment_File_PROT}","$VARS{WorkingDir}$VARS{Alignment_File}",\%DNA_AA_cp);
		if ($ans[0] ne "ok")
		{
			if ($ans[1] eq "user")
			{
				exit_on_error("user_error",join(" ",@ans));			
			}
			elsif ($ans[1] eq "sys")
			{
				exit_on_error("sys_error",join(" ",@ans));
			}
		}
	}
	# CP the scores files calculated for PROT
	#############################################
	system ("mv $VARS{WorkingDir}$VARS{Output_Prefix}_col_col.scr $VARS{WorkingDir}$VARS{Output_Prefix}_col_col.PROT.scr");
	system ("mv $VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_col.scr $VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_col.PROT.scr");
#        system ("mv $VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_seq.scr $VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_seq.PROT.scr");
	system ("mv $VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_seq_pair.scr $VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_seq_pair.PROT.scr");
	system ("mv $VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_res.scr $VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_res.PROT.scr");
	system ("mv $VARS{WorkingDir}$VARS{Output_Prefix}_res_pair.scr $VARS{WorkingDir}$VARS{Output_Prefix}_res_pair.PROT.scr");
	# Convert the scores files calculated for PRO to codons
	#########################################################
	print LOG "Guidance::Convert_to_Codons_Numbering($VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_col.PROT.scr,$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_col.scr);\n";
	my @ans=Guidance::Convert_to_Codons_Numbering("$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_col.PROT.scr","$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_col.scr");
	if ($ans[0] ne "OK"){exit_on_error("sys_error",join(" ",@ans));}
	print LOG "Guidance::Convert_to_Codons_Numbering($VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_res.PROT.scr,$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_res.scr);\n";
	@ans=Guidance::Convert_to_Codons_Numbering("$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_res.PROT.scr","$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_res.scr");
	if ($ans[0] ne "OK"){exit_on_error("sys_error",join(" ",@ans));}
	print LOG "Guidance::Convert_to_Codons_Numbering($VARS{WorkingDir}$VARS{Output_Prefix}_res_pair.PROT.scr,$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair.scr);\n";
	@ans=Guidance::Convert_to_Codons_Numbering("$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair.PROT.scr","$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair.scr");
	if ($ans[0] ne "OK"){exit_on_error("sys_error",join(" ",@ans));}
}	

if ($isServer == 1) {
	Round_Scores_File("$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_seq.scr");
	Round_Scores_File("$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_col.scr");
}


# Return names to the Seq score file
############################################
$VARS{Seq_Scores}=$VARS{Output_Prefix}."_res_pair_seq.scr"."_with_Names";

print LOG "Guidance::codes2name_scoresFile_NEW(\"$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_seq.scr\",\"$VARS{WorkingDir}$VARS{code_fileName}\",\"$VARS{WorkingDir}$VARS{Alignment_File}\",\"$VARS{WorkingDir}$VARS{Seq_Scores}\");\n";
my @ans=Guidance::codes2name_scoresFile_NEW("$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_seq.scr","$VARS{WorkingDir}$VARS{code_fileName}","$VARS{WorkingDir}$VARS{Alignment_File}","$VARS{WorkingDir}$VARS{Seq_Scores}");
if ($ans[0] ne "OK")
{
	print LOG "Guidance::codes2name_scoresFile:".join("",@ans);
}



#remove sites with SP-score < Col sp_cutoff
############################################
$VARS{Alignment_File_without_low_SP_Col}="$VARS{dataset}.".$FORM{MSA_Program}.".Without_low_SP_Col";
$VARS{removed_low_SP_SITE}=$VARS{SeqsFile}.".$VARS{dataset}.".$FORM{MSA_Program}.".Removed_Col";

print LOG "Guidance::removeLowSPsites_NoBioPerl (\"$VARS{WorkingDir}$VARS{Alignment_File}\",\"$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_col.scr\",\"$VARS{WorkingDir}$VARS{Alignment_File_without_low_SP_Col}\",$VARS{SP_COL_CUTOFF},\"$VARS{WorkingDir}$VARS{removed_low_SP_SITE}\");\n";
@ans=Guidance::removeLowSPsites_NoBioPerl ("$VARS{WorkingDir}$VARS{Alignment_File}","$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_col.scr","$VARS{WorkingDir}$VARS{Alignment_File_without_low_SP_Col}",$VARS{SP_COL_CUTOFF},"$VARS{WorkingDir}$VARS{removed_low_SP_SITE}");

if ($ans[0]eq "OK")
{
    $VARS{REMOVED_SITES}=$ans[1];
    $VARS{MSA_LENGTH}=$ans[2];
}
print LOG "REMOVED_SITES:$VARS{REMOVED_SITES}\n";
print LOG "MSA_LENGTH:$VARS{MSA_LENGTH}\n";
$VARS{Alignment_File_without_low_SP_Col_with_Names}=$VARS{Alignment_File_without_low_SP_Col}.".With_Names";
if (-s "$VARS{WorkingDir}$VARS{Alignment_File_without_low_SP_Col}" > 0) # Not EMPTY
{
	my @ans=Guidance::codes2nameFastaFrom1("$VARS{WorkingDir}$VARS{Alignment_File_without_low_SP_Col}","$VARS{WorkingDir}$VARS{code_fileName}","$VARS{WorkingDir}$VARS{Alignment_File_without_low_SP_Col_with_Names}");
	if ($ans[0] ne "OK") {exit_on_error("sys_error","Guidance::codes2nameFastaFrom1: Guidance::codes2nameFastaFrom1(\"$VARS{WorkingDir}$VARS{Alignment_File_without_low_SP_Col}\",\"$VARS{WorkingDir}$VARS{code_fileName}\",\"$VARS{WorkingDir}$VARS{Alignment_File_without_low_SP_Col_with_Names}\") failed:",join("",@ans),"\n");}
}

if ($VARS{Z_Col_Cutoff} ne "NA")
{
    #remove sites with SP-score < Col sp_cutoff if Z is below cutoff
    ##################################################################
	$VARS{Alignment_File_without_low_SP_Z_Col}="$VARS{dataset}.".$FORM{MSA_Program}.".Without_low_SP_and_low_Z_Col";
	$VARS{removed_low_SP_Z_SITE}=$VARS{SeqsFile}."$VARS{dataset}.".$FORM{MSA_Program}.".Removed_SP_and_Z_Col";
	
	print LOG "Guidance::removeLowSPsites_NoBioPerl_Consider_Z (\"$VARS{WorkingDir}$VARS{Alignment_File}\",\"$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_col.scr\",\"$VARS{WorkingDir}$VARS{Alignment_File_without_low_SP_Z_Col}\",$VARS{SP_COL_CUTOFF},$VARS{Z_Col_Cutoff},\"$VARS{WorkingDir}$VARS{removed_low_SP_Z_SITE}\");\n";
	@ans=Guidance::removeLowSPsites_NoBioPerl_Consider_Z ("$VARS{WorkingDir}$VARS{Alignment_File}","$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_col.scr","$VARS{WorkingDir}$VARS{Alignment_File_without_low_SP_Z_Col}",$VARS{SP_COL_CUTOFF},$VARS{Z_Col_Cutoff},"$VARS{WorkingDir}$VARS{removed_low_SP_Z_SITE}");
	
	if ($ans[0]eq "OK")
	{
		$VARS{REMOVED_Z_SITES}=$ans[1];
		$VARS{MSA_LENGTH}=$ans[2];
	}
	print LOG "REMOVED_SITES_SP_AND_Z:$VARS{REMOVED_Z_SITES}\n";
	print LOG "MSA_LENGTH:$VARS{MSA_LENGTH}\n";
	$VARS{Alignment_File_without_low_SP_Z_Col_with_Names}=$VARS{Alignment_File_without_low_SP_Z_Col}.".With_Names";
	if (-s "$VARS{WorkingDir}$VARS{Alignment_File_without_low_SP_Z_Col}" > 0) # Not EMPTY
	{
		my @ans=Guidance::codes2nameFastaFrom1("$VARS{WorkingDir}$VARS{Alignment_File_without_low_SP_Z_Col}","$VARS{WorkingDir}$VARS{code_fileName}","$VARS{WorkingDir}$VARS{Alignment_File_without_low_SP_Z_Col_with_Names}");
		if ($ans[0] ne "OK") {exit_on_error("sys_error","Guidance::codes2nameFastaFrom1: Guidance::codes2nameFastaFrom1(\"$VARS{WorkingDir}$VARS{Alignment_File_without_low_SP_Z_Col}\",\"$VARS{WorkingDir}$VARS{code_fileName}\",\"$VARS{WorkingDir}$VARS{Alignment_File_without_low_SP_Z_Col_with_Names}\") failed:",join("",@ans),"\n");}
	}
}

# Prepare for the Plot
##########################
print LOG "ans=Convert_to_CSV(\"$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_col.scr\",\"$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_col.scr.csv\");\n";
@ans=Convert_to_CSV("$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_col.scr","$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_col.scr.csv");
print LOG "ANS:",join (" ",@ans),"\n";

$VARS{res_pair_res_html_file}=$VARS{dataset}.".".$FORM{MSA_Program}.".Guidance_res_pair_res.html";

print LOG "Guidance::printColoredAlignment_With_CSS(\"$VARS{WorkingDir}$VARS{Alignment_File}\",\"$VARS{WorkingDir}$VARS{res_pair_res_html_file}\",\"$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_res.scr\",\"$VARS{WorkingDir}$VARS{code_fileName}\",\"$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_col.scr.csv\",$FORM{PROGRAM},\"$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_seq.scr\");\n";
@ans=Guidance::printColoredAlignment_With_CSS("$VARS{WorkingDir}$VARS{Alignment_File}","$VARS{WorkingDir}$VARS{res_pair_res_html_file}","$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_res.scr","$VARS{WorkingDir}$VARS{code_fileName}","$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_col.scr.csv",$FORM{PROGRAM},"$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_seq.scr");
print LOG join("",@ans);

# bring back names to MSA
###########################
$VARS{Alignment_File_With_Names}=$VARS{Alignment_File}.".With_Names";
@ans=Guidance::codes2nameFastaFrom1("$VARS{WorkingDir}$VARS{Alignment_File}","$VARS{WorkingDir}$VARS{code_fileName}","$VARS{WorkingDir}$VARS{Alignment_File_With_Names}");
if ($ans[0] ne "OK") {exit_on_error("sys_error","Guidance::codes2nameFastaFrom1: Guidance::codes2nameFastaFrom1($VARS{WorkingDir}$VARS{Alignment_File},$VARS{WorkingDir}$VARS{code_fileName},$VARS{WorkingDir}$VARS{Alignment_File_With_Names}) failed:",join("",@ans),"\n");}

#remove seq with SP-score < Seq sp_cutoff
############################################
$VARS{Seq_File_without_low_SP_SEQ}=$VARS{SeqsFile}.".Without_low_SP_Seq";
$VARS{removed_low_SP_SEQ}=$VARS{SeqsFile}.".Removed_Seq";
$VARS{Seq_File_without_low_SP_SEQ_with_Names}=$VARS{Seq_File_without_low_SP_SEQ}.".With_Names";
$VARS{removed_low_SP_SEQ_With_Names}=$VARS{removed_low_SP_SEQ}.".With_Names";

print LOG "Guidance::removeLowSPseq (\"$VARS{WorkingDir}$VARS{Alignment_File}\",\"$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_seq.scr\",\"$VARS{WorkingDir}$VARS{Seq_File_without_low_SP_SEQ}\",$VARS{SP_SEQ_CUTOFF},\"$VARS{WorkingDir}$VARS{removed_low_SP_SEQ}\");\n";
Guidance::removeLowSPseq ("$VARS{WorkingDir}$VARS{Alignment_File}","$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_seq.scr","$VARS{WorkingDir}$VARS{Seq_File_without_low_SP_SEQ}",$VARS{SP_SEQ_CUTOFF},"$VARS{WorkingDir}$VARS{removed_low_SP_SEQ}");
if (-s  "$VARS{WorkingDir}$VARS{Seq_File_without_low_SP_SEQ}" > 0) # NOT EMPTY
{
	my @ans=Guidance::codes2nameFastaFrom1("$VARS{WorkingDir}$VARS{Seq_File_without_low_SP_SEQ}","$VARS{WorkingDir}$VARS{code_fileName}","$VARS{WorkingDir}$VARS{Seq_File_without_low_SP_SEQ_with_Names}");
	if ($ans[0] ne "OK") {exit_on_error("sys_error","Guidance::codes2nameFastaFrom1: Guidance::codes2nameFastaFrom1(\"$VARS{WorkingDir}$VARS{Seq_File_without_low_SP_SEQ}\",\"$VARS{WorkingDir}$VARS{code_fileName}\",\"$VARS{WorkingDir}$VARS{Seq_File_without_low_SP_SEQ_with_Names}\") failed:",join("",@ans),"\n");}
}

if (-s "$VARS{WorkingDir}$VARS{removed_low_SP_SEQ}" > 0) # Seq were removed
{
	my  @ans=Guidance::codes2nameFastaFrom1("$VARS{WorkingDir}$VARS{removed_low_SP_SEQ}","$VARS{WorkingDir}$VARS{code_fileName}","$VARS{WorkingDir}$VARS{removed_low_SP_SEQ_With_Names}");
	if ($ans[0] ne "OK") {exit_on_error("sys_error","Guidance::codes2nameFastaFrom1: Guidance::codes2nameFastaFrom1(Guidance::codes2nameFastaFrom1(\"$VARS{WorkingDir}$VARS{removed_low_SP_SEQ}\",\"$VARS{WorkingDir}$VARS{code_fileName}\",\"$VARS{WorkingDir}$VARS{removed_low_SP_SEQ_With_Names}\") failed:".join("",@ans)."\n");}
}

#remove seq with SP-score < Seq sp_cutoff if Z<(-Z_Cutoff)
##############################################################
if ($VARS{Z_Seq_Cutoff} ne "NA")
{
	$VARS{Seq_File_without_low_SP_and_Z_SEQ}=$VARS{SeqsFile}.".Without_low_Z_and_low_SP_Seq";
	$VARS{removed_low_SP_Z_SEQ}=$VARS{SeqsFile}.".Removed_Z_SP_Seq";
	$VARS{Seq_File_without_low_SP_and_Z_SEQ_with_Names}=$VARS{Seq_File_without_low_SP_and_Z_SEQ}.".With_Names";
	$VARS{removed_low_SP_Z_SEQ_With_Names}=$VARS{removed_low_SP_Z_SEQ}.".With_Names";
	
	print LOG "Guidance::removeLowSPseq_Consider_Z (\"$VARS{WorkingDir}$VARS{Alignment_File}\",\"$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_seq.scr\",\"$VARS{WorkingDir}$VARS{Seq_File_without_low_SP_Z_SEQ}\",$VARS{SP_SEQ_CUTOFF},$VARS{Z_Seq_Cutoff},\"$VARS{WorkingDir}$VARS{removed_low_SP_Z_SEQ}\");\n";
	Guidance::removeLowSPseq_Consider_Z ("$VARS{WorkingDir}$VARS{Alignment_File}","$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_seq.scr","$VARS{WorkingDir}$VARS{Seq_File_without_low_SP_Z_SEQ}",$VARS{SP_SEQ_CUTOFF},$VARS{Z_Seq_Cutoff},"$VARS{WorkingDir}$VARS{removed_low_SP_Z_SEQ}");
	if (-s  "$VARS{WorkingDir}$VARS{Seq_File_without_low_SP_Z_SEQ}" > 0) # NOT EMPTY
	{
		my @ans=Guidance::codes2nameFastaFrom1("$VARS{WorkingDir}$VARS{Seq_File_without_low_SP_Z_SEQ}","$VARS{WorkingDir}$VARS{code_fileName}","$VARS{WorkingDir}$VARS{Seq_File_without_low_SP_SEQ_with_Names}");
		if ($ans[0] ne "OK") {exit_on_error("sys_error","Guidance::codes2nameFastaFrom1: Guidance::codes2nameFastaFrom1(\"$VARS{WorkingDir}$VARS{Seq_File_without_low_SP_SEQ}\",\"$VARS{WorkingDir}$VARS{code_fileName}\",\"$VARS{WorkingDir}$VARS{Seq_File_without_low_SP_SEQ_with_Names}\") failed:",join("",@ans),"\n");}
	}
	
	if (-s "$VARS{WorkingDir}$VARS{removed_low_SP_SEQ}" > 0) # Seq were removed
	{
		my  @ans=Guidance::codes2nameFastaFrom1("$VARS{WorkingDir}$VARS{removed_low_SP_Z_SEQ}","$VARS{WorkingDir}$VARS{code_fileName}","$VARS{WorkingDir}$VARS{removed_low_SP_Z_SEQ_With_Names}");
		if ($ans[0] ne "OK") {exit_on_error("sys_error","Guidance::codes2nameFastaFrom1: Guidance::codes2nameFastaFrom1(Guidance::codes2nameFastaFrom1(\"$VARS{WorkingDir}$VARS{removed_low_SP_Z_SEQ}\",\"$VARS{WorkingDir}$VARS{code_fileName}\",\"$VARS{WorkingDir}$VARS{removed_low_SP_Z_SEQ_With_Names}\") failed:".join("",@ans)."\n");}
	}
}
# Make JalView outputs
$VARS{'JalView_Page'}=$VARS{Output_Prefix}."_JalView.html";
$VARS{'JalView_Features'}=$VARS{Output_Prefix}."_JalView_Features";
$VARS{'JalView_Annotations'}=$VARS{Output_Prefix}."_JalView_Annot";
if ($isServer == 1) # JalView Applet output (for web-server)
{
	print LOG "Prepare JalView outputs: Guidance::make_JalView_output($VARS{WorkingDir}$VARS{'JalView_Page'},$VARS{WorkingDir},$VARS{run_url},$VARS{Alignment_File},$VARS{Alignment_File_With_Names},$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_res.scr,$VARS{'JalView_Features'},$VARS{WorkingDir}$VARS{code_fileName},$VARS{'JalView_Annotations'},$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_col.scr.csv,$FORM{PROGRAM} scores);\n";
	Guidance::make_JalView_output($VARS{WorkingDir}.$VARS{'JalView_Page'},$VARS{WorkingDir},$VARS{run_url},$VARS{Alignment_File},$VARS{Alignment_File_With_Names},"$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_res.scr",$VARS{'JalView_Features'},"$VARS{WorkingDir}$VARS{code_fileName}",$VARS{'JalView_Annotations'},"$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_col.scr.csv",$FORM{PROGRAM}." scores");
}
else # creates only the Features and annotation graph without the applet page
{
	print LOG "Guidance::make_Jalview_Color_MSA(\"$VARS{WorkingDir}$VARS{Alignment_File}\",\"$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_res.scr\",\"$VARS{WorkingDir}$VARS{'JalView_Features'}\",\"$VARS{WorkingDir}$VARS{code_fileName}\");\n";
	Guidance::make_Jalview_Color_MSA("$VARS{WorkingDir}$VARS{Alignment_File}","$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_res.scr","$VARS{WorkingDir}$VARS{'JalView_Features'}","$VARS{WorkingDir}$VARS{code_fileName}");
	print LOG "Guidance::make_Jalview_AnnotationGraph(\"$VARS{WorkingDir}$VARS{'JalView_Annotations'}\",\"$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_col.scr.csv\",$FORM{PROGRAM}.\" scores\");\n";
	Guidance::make_Jalview_AnnotationGraph("$VARS{WorkingDir}$VARS{'JalView_Annotations'}","$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_col.scr.csv",$FORM{PROGRAM}." scores");
}
# FLAG THAT RUN FINISHED OK
if ($isServer == 1)
{
	open (ENDS_OK,">$VARS{WorkingDir}GUIDANCE_$VARS{run_number}.END_OK");
	close (ENDS_OK);
}
else
{
	open (ENDS_OK,">$VARS{WorkingDir}ENDS_OK");
	close (ENDS_OK);
}
if ((($FORM{PROGRAM} eq "GUIDANCE") or ($FORM{PROGRAM} eq "GUIDANCE2")) and ($FORM{Seq_Type} eq "Codons"))
{
	if (((-s "$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair.PROT.scr")/1048576)>100)
	{
		system ("zip -j $VARS{WorkingDir}$VARS{Output_Prefix}_res_pair.PROT.scr.zip $VARS{WorkingDir}$VARS{Output_Prefix}_res_pair.PROT.scr; rm -f $VARS{WorkingDir}$VARS{Output_Prefix}_res_pair.PROT.scr");
	}
	if (((-s "$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_res.PROT.scr")/1048576)>100)
	{
		system ("zip -j $VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_res.PROT.scr.zip $VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_res.PROT.scr; rm -f $VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_res.PROT.scr");
	}
}
# creat tar of alternative alignment
if ($FORM{PROGRAM} eq "GUIDANCE2")
{
    system ("mkdir $VARS{GUIDANCE2_MSAs_Dir}/wNames/");
    $VARS{GUIDANCE2_MSAs_Name_Pattern}=$VARS{GUIDANCE2_MSAs_Dir}."*.fasta";
    # print "$VARS{GUIDANCE2_MSAs_Name_Pattern}\n";
    foreach my $AltMSA (<$VARS{GUIDANCE2_MSAs_Name_Pattern}>)
    {
		# print "MSA: $AltMSA\n";
		my $AltMSA_No_Path=basename($AltMSA);
		if ($FORM{Seq_Type} eq "Codons")
		{
			# Copy the hash - AA_to_DNA_aligned change the reference...
			my %DNA_AA_cp=();
			foreach my $key (keys %DNA_AA)
			{
				my @Arr;
				$Arr[0]=$DNA_AA{$key}->[0];
				$Arr[1]=$DNA_AA{$key}->[1];
				$DNA_AA_cp{$key}=[@Arr];
			}
			# print "codonAlign::AA_to_DNA_aligned(\"$AltMSA\",\"$AltMSA.AA\",\%DNA_AA_cp);\n";#<STDIN>;
			codonAlign::AA_to_DNA_aligned("$AltMSA","$AltMSA.AA",\%DNA_AA_cp);
			# print "Guidance::codes2nameFastaFrom1(\"$AltMSA.AA\",\"$VARS{WorkingDir}$VARS{code_fileName}\",\"$VARS{GUIDANCE2_MSAs_Dir}/wNames/$AltMSA_No_Path\");\n";#<STDIN>; 
			Guidance::codes2nameFastaFrom1("$AltMSA.AA","$VARS{WorkingDir}$VARS{code_fileName}","$VARS{GUIDANCE2_MSAs_Dir}/wNames/$AltMSA_No_Path"); 
			unlink ("$AltMSA.AA");
		}
		else
		{
			# print "Guidance::codes2nameFastaFrom1(\"$AltMSA\",\"$VARS{WorkingDir}$VARS{code_fileName}\",\"$VARS{GUIDANCE2_MSAs_Dir}/wNames/$AltMSA_No_Path\")\n";
			Guidance::codes2nameFastaFrom1("$AltMSA","$VARS{WorkingDir}$VARS{code_fileName}","$VARS{GUIDANCE2_MSAs_Dir}/wNames/$AltMSA_No_Path");
		}
    }
    # print "cd $VARS{GUIDANCE2_MSAs_Dir}/wNames/;tar -czf $VARS{WorkingDir}$VARS{Output_Prefix}_AlternativeMSA.tar.gz *.fasta\n";
    system ("cd $VARS{GUIDANCE2_MSAs_Dir}/wNames/;tar -czf $VARS{WorkingDir}$VARS{Output_Prefix}_AlternativeMSA.tar.gz *.fasta");
}
if (($FORM{PROGRAM} eq "GUIDANCE") or ($FORM{PROGRAM} eq "GUIDANCE2") or ($FORM{PROGRAM} eq "GUIDANCE3"))
{
# tar and remove the BP dir
	$cmd="cd $VARS{WorkingDir};tar -czf $VARS{Output_Prefix}_BP_Dir.tar.gz ./BP";
#print "$cmd\n";
	system ($cmd);
	system ("rm -r -f $VARS{BootStrap_Dir}");
}
if ($FORM{PROGRAM} eq "HoT")
{
# tar and remove the BP dir
	$cmd="cd $VARS{WorkingDir};tar -czf $VARS{dataset}".".$FORM{MSA_Program}_HoT_Dir.tar.gz ./$VARS{HoT_MSAs_Dir} ./$VARS{dataset}_cos_$VARS{HoT_MSA_Program}";
	# print "$cmd\n";
	system ($cmd);
	system ("rm -r -f $VARS{WorkingDir}$VARS{dataset}_cos_$VARS{HoT_MSA_Program}");
	system ("rm -r -f $VARS{WorkingDir}$VARS{HoT_MSAs_Dir}");
}
# Output
###########################
if ($isServer == 1) {
	print OUTPUT "\n<H1><center><a name=finish>GUIDANCE calculation is finished:</a></center></H1>\n";
		
		print OUTPUT "<h3><b>Results:</b></h3>\n";
		print_message_to_output("<A HREF='$VARS{res_pair_res_html_file}' TARGET=results_Colored>MSA Colored according to the confidence score</A></font>&nbsp;&nbsp;&nbsp;<A HREF='$VARS{'JalView_Page'}' target=\"_blank\">(open with Jalview)</A><br>");
	    my $MEAN_RES_PAIR_SCORE=extract_MEAN_RES_PAIR_SCORE("$VARS{WorkingDir}$VARS{Output_Prefix}_msa.scr");
	    print_message_to_output("GUIDANCE alignment score: $MEAN_RES_PAIR_SCORE</font><br>");
		print OUTPUT "<h3><b>Output Files:</b></h3>\n";
		print_message_to_output("<A HREF='$VARS{Alignment_File_With_Names}' TARGET=MSA>MSA file</A></font><br>");
#print_message_to_output("<A HREF='$VARS{Output_Prefix}_col_col.scr' TARGET=results_col_col>GUIDANCE column column scores</A></font><br>");
		print_message_to_output("<A HREF='$VARS{Output_Prefix}_res_pair_col.scr' TARGET=results_res_pair_col>GUIDANCE column scores</A></font><br>");
#print_message_to_output("<A HREF='$VARS{Output_Prefix}_res_pair_seq.scr' TARGET=results_res_pair_seq>GUIDANCE sequence scores</A></font><br>");
		print_message_to_output("<A HREF='$VARS{Seq_Scores}' TARGET=results_res_pair_seq>GUIDANCE sequence scores</A></font><br>");
#print_message_to_output("<A HREF='$VARS{Output_Prefix}_res_pair_seq_pair.scr' TARGET=results_res_pair_seq_pair>GUIDANCE residue pair sequence pair scores</A></font><br>");
	if (((-s "$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_res.scr")/1048576)>100)
	{
		system ("zip -j $VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_res.scr.zip $VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_res.scr; rm -f $VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_res.scr");
		print_message_to_output("<A HREF='$VARS{Output_Prefix}_res_pair_res.scr.zip' TARGET=results_res_pair_res>GUIDANCE residue scores</A></font><br>");
	}
	else
	{
		print_message_to_output("<A HREF='$VARS{Output_Prefix}_res_pair_res.scr' TARGET=results_res_pair_res>GUIDANCE residue scores</A></font><br>");
	}
	if (((-s "$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair.scr")/1048576)>100)
	{
		system ("zip -j $VARS{WorkingDir}$VARS{Output_Prefix}_res_pair.scr.zip $VARS{WorkingDir}$VARS{Output_Prefix}_res_pair.scr; rm -f $VARS{WorkingDir}$VARS{Output_Prefix}_res_pair.scr");
		print_message_to_output("<A HREF='$VARS{Output_Prefix}_res_pair.scr.zip' TARGET=results_res_pair_res>GUIDANCE residue pair scores</A></font><br>");
	}
	else
	{
		print_message_to_output("<A HREF='$VARS{Output_Prefix}_res_pair.scr' TARGET=results_res_pair_res>GUIDANCE residue pair scores</A></font><br>");
	}
		
	my $cmd1;
	if (($FORM{Seq_Type} eq "AminoAcids") or ($FORM{Seq_Type} eq "Codons")) {$cmd1.="aa";}
	elsif ($FORM{Seq_Type} eq "Nucleotides") {$cmd1.="nuc";}
	
	print_message_to_output("<form enctype=\"multipart/form-data\" action=\"/mask.php\" method=\"post\"> <font color=red face=\"Comic Sans MS\">New feature </font> - Mask specific residues below a certain cutoff: <input type=\"hidden\" name=\"run_num\" value=\"$VARS{run_number}\"/> <input type=\"hidden\" name=\"VARS\" value=\"$stored_data_file\"/>  ".print_remove_site_selection_mask("$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_res.scr",$VARS{MSA_LENGTH})."<input type=\"hidden\" name=\"type_a\" value=\"$cmd1\"/> <input type=\"submit\" value=\"Mask Residues\"/> <font size =\"-1\">(see <a href=\"$overview_URL#Mask_RES\" target=\"_blank\">help</a></font>)</form></li>\n</ul>");
	
	print_message_to_output("<form enctype=\"multipart/form-data\" action=\"/remove_pos.php\" method=\"post\">\n Remove unreliable columns below confidence score ".print_remove_site_selection_box("$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_col.scr",$VARS{MSA_LENGTH}));
###print OUTPUT "<ul>";
###print_message_to_output("<A HREF='$VARS{Alignment_File_without_low_SP_Col_with_Names}' TARGET=_blank>The MSA after the removal of unreliable columns (below $VARS{SP_COL_CUTOFF})</A><font size=-1> (see list of removed columns <A HREF='$VARS{removed_low_SP_SITE}' TARGET=_blank>here</A>)</font><br>");
###print OUTPUT "</ul>\n";
	
	if (!-e "$VARS{WorkingDir}$VARS{Alignment_File_without_low_SP_Col_with_Names}")
	{
		print OUTPUT "<ul>";
		print_message_to_output("<font color='red'><B>ATTENTION:</font></B> All positions had score below ".sprintf("%.3f", $VARS{SP_COL_CUTOFF})."<br>");
		print OUTPUT "</ul>\n";
	}
	elsif (!-e "$VARS{WorkingDir}$VARS{removed_low_SP_SITE}")
	{
		print OUTPUT "<ul>";
		print_message_to_output(" All positions had score higher than ".sprintf("%.3f", $VARS{SP_COL_CUTOFF})."<br>");
		print OUTPUT "</ul>\n";
	}
	else
	{
		print OUTPUT "<ul>";
		print_message_to_output("<A HREF='$VARS{Alignment_File_without_low_SP_Col_with_Names}' TARGET=_blank>The MSA after the removal of unreliable columns (below $VARS{SP_COL_CUTOFF})</A><font size=-1> (see list of removed columns <A HREF='$VARS{removed_low_SP_SITE}' TARGET=_blank>here</A>)</font><br>");
###	print_message_to_output("<A HREF='$VARS{Alignment_File_without_low_SP_Col_with_Names}' TARGET=Alignmnet_without_low_SP>The MSA after the removal of unreliable columns (below ".sprintf("%.3f", $VARS{SP_COL_CUTOFF}).")</A><font size=-1> (see list of removed columns <A HREF='$VARS{removed_low_SP_SITE}' TARGET=_blank>here</A>)</font><br>"); 
		print OUTPUT "</ul>\n";
	}
	print LOG "MSA_Depth($VARS{WorkingDir}$VARS{Alignment_File});\n";
	$VARS{MSA_Depth}=MSA_Depth("$VARS{WorkingDir}$VARS{Alignment_File}");
	print LOG "MSA DEPTH:$VARS{MSA_Depth}\n";
	print_message_to_output("<form enctype=\"multipart/form-data\" action=\"/remove_seq.php\" method=\"post\">\n Remove unreliable sequences below confidence score ".print_remove_seq_selection_box("$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_seq.scr",$VARS{MSA_Depth}));
	if (!-e "$VARS{WorkingDir}$VARS{Seq_File_without_low_SP_SEQ_with_Names}")
	{
		print OUTPUT "<ul>";
		print_message_to_output("<font color='red'><B>ATTENTION:</font></B> All sequences had score below  $VARS{SP_SEQ_CUTOFF}<br>");
		print OUTPUT "</ul>";
	}
	elsif (!-e "$VARS{WorkingDir}$VARS{removed_low_SP_SEQ_With_Names}")
	{
		print OUTPUT "<ul>";
		print_message_to_output(" All sequences had score higher than $VARS{SP_SEQ_CUTOFF}<br>");
		print OUTPUT "</ul>";
	}
	else 
	{
		print OUTPUT "<ul>";
		if ($FORM{'Redirect_From_MAFFT'}==1) {print_message_to_output("<A HREF='$VARS{Seq_File_without_low_SP_SEQ_with_Names}' TARGET=_blank>The input sequences after the removal of unreliable sequences (with confidence score below $VARS{SP_SEQ_CUTOFF})</A><font size=-1> (see list of removed sequences <A HREF='$VARS{removed_low_SP_SEQ_With_Names}' TARGET=_blank>here</A></font>)&nbsp;&nbsp;&nbsp;<INPUT TYPE=\"BUTTON\" VALUE=\"run GUIDANCE on the confidently-aligned sequences only\" ONCLICK=\"var answer = confirm('ATTENTION: Running GUIDANCE on the confidently-aligned sequences only, ignores the parameters used for the original run on MAFFT server. It is therefore recommended to adjust these parameters or aligning the confidently-aligned sequences on MAFFT server and run GUIDANCE again from there');if (answer){window.open('$home_URL/index_rerun.php?run=$VARS{run_number}&file=$VARS{Seq_File_without_low_SP_SEQ_with_Names}')}\"><br>");}
		else {print_message_to_output("<A HREF='$VARS{Seq_File_without_low_SP_SEQ_with_Names}' TARGET=_blank>The input sequences after the removal of unreliable sequences (with confidence score below $VARS{SP_SEQ_CUTOFF})</A><font size=-1> (see list of removed sequences <A HREF='$VARS{removed_low_SP_SEQ_With_Names}' TARGET=_blank>here</A></font>)&nbsp;&nbsp;&nbsp;<INPUT TYPE=\"BUTTON\" VALUE=\"run GUIDANCE on the confidently-aligned sequences only\" ONCLICK=\"window.open('$home_URL/index_rerun.php?run=$VARS{run_number}&file=$VARS{Seq_File_without_low_SP_SEQ_with_Names}')\"><br>");}
		print OUTPUT "</ul>";
	}
	print_message_to_output("<A HREF='$VARS{Output_Prefix}_AlternativeMSA.tar.gz'> The set of alternative MSAs<br>");
#finish_output();
	print OUTPUT "<br>\n";
	close OUTPUT;
	sleep (5);
	print LOG "update_output_that_run_finished(\"$VARS{WorkingDir}$VARS{output_page}\");\n";
	&update_output_that_run_finished("$VARS{WorkingDir}$VARS{output_page}");
	&send_finish_email_to_user();
	Prepare_rerun_param();
#stop_reload("$VARS{working_dir}/$VARS{output_page}");
	&store_data();
	close LOG;
}

sub run_Guidance2
{
	$VARS{Scoring_Alignments_Dir}=$VARS{GUIDANCE2_MSAs_Dir};
	
	# INIT
	if (($FORM{MSA_Program} eq "MAFFT") or ($FORM{MSA_Program} eq "MAFFT_LINSI")) 
	{
		$VARS{HoT_MSA_Program}="MFT";
		$VARS{HoT_MSA_Program_path}=$mafft_prog;
		my $check_mafft_profile=`which $VARS{HoT_MSA_Program_path}-profile`;
		if ($check_mafft_profile=~/Command not found/){die "It seems that $VARS{HoT_MSA_Program_path}-profile is not properly installed or found in PATH variable. Please fix that and/or provide GUIDANCE with the full path to mafft installation using the --mafft argument\n";}
	}
	if ($FORM{MSA_Program} eq "MUSCLE") {exit_on_error('user_error', "GUIDANCE2 currently does not support MUSCLE, please run GUIDANCE<br>");}
	if ($FORM{MSA_Program} eq "PAGAN") {exit_on_error('user_error', "GUIDANCE2 currently does not support PAGAN, please run GUIDANCE<br>");}
	elsif ($FORM{MSA_Program} eq "CLUSTALW") {$VARS{HoT_MSA_Program}="CLW";$VARS{HoT_MSA_Program_path}=$clustalw_prog;}
    elsif ($FORM{MSA_Program} eq "PRANK") {$VARS{HoT_MSA_Program}="PRK";$VARS{HoT_MSA_Program_path}=$prank_prog;}
	my $OP_DistFile;
	my $EP_DistFile;
	my $Num_of_Aln_from_HoT_per_Run=4;# number of alignment to sample out of HoT (GUIDANCE2) pertubed alignment (per split,tree and op value, 4 or 8 HoT outputs are created)
	#HoT assumes that all the sequences are upper case
	print LOG "convert_fs_to_upper_case($VARS{WorkingDir}$VARS{codded_seq_fileName})\n";
	convert_fs_to_upper_case("$VARS{WorkingDir}$VARS{codded_seq_fileName}");
	
	# Align
	##############
	if ($FORM{userMSA_File} eq "") #align if user did not supplied alignment
	{
		Align();
		convert_fs_to_upper_case("$VARS{WorkingDir}$VARS{Alignment_File}");# HOT ASSUME THAT THE SEQUENCES ARE ALL UPPER CASE, SO WE CONVERT THE ALN TO UPPER CASE
	}
	# TO DO: handle the adjustdirection

	# BootStrap Trees
	##################
	Bootstrap_Trees();

	# pull out the trees
	######################
	my $numUniqueTrees="";
	my @numRepeats4UniqueTree=();
	if (uc $VARS{BBL} eq "YES")
	{
		print LOG "Guidance::pullOutBPtrees_BBL($VARS{WorkingDir},$VARS{dataset},$FORM{Bootstraps},$FORM{MSA_Program});\n";
		@ans=Guidance::pullOutBPtrees_BBL($VARS{WorkingDir},$VARS{dataset},$FORM{Bootstraps},$FORM{MSA_Program});
		unless ($ans[0] eq "ok") {exit_on_error("sys_error","Guidance::pullOutBPtrees_BBL: ".join (" ",@ans)."\n");}
		$numUniqueTrees=$ans[1] if ($FORM{MSA_Program} ne "MAFFT");
		@numRepeats4UniqueTree=@{$ans[2]} if ($FORM{MSA_Program} ne "MAFFT");
	}
	else
	{
		print LOG "Guidance::pullOutBPtrees($VARS{WorkingDir},$VARS{dataset},$FORM{Bootstraps},$FORM{MSA_Program});\n";
		@ans=Guidance::pullOutBPtrees($VARS{WorkingDir},$VARS{dataset},$FORM{Bootstraps},$FORM{MSA_Program});
		unless ($ans[0] eq "ok") {exit_on_error("sys_error","Guidance::pullOutBPtrees: ".join (" ",@ans)."\n");}
		$numUniqueTrees=$ans[1] if ($FORM{MSA_Program} ne "MAFFT");
		@numRepeats4UniqueTree=@{$ans[2]} if ($FORM{MSA_Program} ne "MAFFT");
	}

	# Convert trees to MAFFT format
	#################################
	if (($FORM{MSA_Program} eq "MAFFT") or ($FORM{MSA_Program} eq "MAFFT_LINSI")) #FOR MAFFT BUILED ALIGNMENT FOR ALL TREES
	{
		# if addfragments is used first prune the trees for the the 'core' alignment and afterwards will use the created 'BP core' and the 'full tree' to create the full alignment
		if ($VARS{align_param}=~/addfragments/)
		{
			$VARS{prune_BootStrap_Dir}=$VARS{BootStrap_Dir}."PRUNE_BP_FOR_CORE_ALN/";
			# prune trees
			unless (-e $VARS{prune_BootStrap_Dir}){system ("mkdir $VARS{prune_BootStrap_Dir}");}
			foreach my $full_tree_dir (<$VARS{BootStrap_Dir}tree*>)
			{
				#print $full_tree_dir;<STDIN>;
				if (-d $full_tree_dir)
				{
					if ($full_tree_dir=~/(tree_[0-9]+)/)
					{
						my $tree_num=$1;
						unless (-e "$VARS{prune_BootStrap_Dir}$tree_num"){system ("mkdir $VARS{prune_BootStrap_Dir}$tree_num");}
						my $FullTree="$full_tree_dir/$VARS{dataset}.$FORM{MSA_Program}.semphy.$tree_num";
						my $PrunedTree="$VARS{prune_BootStrap_Dir}/$tree_num/$VARS{dataset}.$FORM{MSA_Program}.semphy.$tree_num"."CORE";
						system ("$remove_taxa_prog $FullTree $VARS{WorkingDir}$VARS{fragments_codes} $PrunedTree");
					}
				}
			}
			# root pruned trees
			print LOG "Guidance::root_BP_trees($VARS{prune_BootStrap_Dir},$VARS{dataset}, $FORM{MSA_Program}, $FORM{Bootstraps},\"\",$VARS{rooting_type});\n"; # TO DO: Consider to be $numUniqueTrees
			@ans=Guidance::root_BP_trees($VARS{prune_BootStrap_Dir},$VARS{dataset}, $FORM{MSA_Program}, $FORM{Bootstraps},"",$VARS{rooting_type});
			unless ($ans[0] eq "ok"){exit_on_error("sys_error","Guidance::root_BP_trees: ".join(" ",@ans)."\n");}
			if (-z "$VARS{prune_BootStrap_Dir}tree_".($FORM{Bootstraps}-1)."/$VARS{dataset}.$FORM{MSA_Program}.semphy.tree_".($FORM{Bootstraps}-1)."CORE.rooted" or !-e "$VARS{prune_BootStrap_Dir}tree_".($FORM{Bootstraps}-1)."/$VARS{dataset}.$FORM{MSA_Program}.semphy.tree_".($FORM{Bootstraps}-1)."CORE.rooted")
			{
				exit_on_error("sys_error","$VARS{prune_BootStrap_Dir}tree_".($FORM{Bootstraps}-1)."/$VARS{dataset}.$FORM{MSA_Program}.semphy.tree_".($FORM{Bootstraps}-1)."CORE.rooted"." does not exist/empty\n");  # TO DO: Consider to be $numUniqueTrees instead Bootstraps
			}
		}
		# prepare the trees
		print LOG "Guidance::root_BP_trees($VARS{BootStrap_Dir},$VARS{dataset}, $FORM{MSA_Program}, $FORM{Bootstraps},\"\",$VARS{rooting_type});\n"; # TO DO: consider $numUniqueTrees
		@ans=Guidance::root_BP_trees($VARS{BootStrap_Dir},$VARS{dataset}, $FORM{MSA_Program}, $FORM{Bootstraps},"",$VARS{rooting_type});
		unless ($ans[0] eq "ok"){exit_on_error("sys_error","Guidance::root_BP_trees: ".join(" ",@ans)."\n");}
		if (-z "$VARS{BootStrap_Dir}tree_".($FORM{Bootstraps}-1)."/$VARS{dataset}.$FORM{MSA_Program}.semphy.tree_".($FORM{Bootstraps}-1).".rooted" or !-e "$VARS{BootStrap_Dir}tree_".($FORM{Bootstraps}-1)."/$VARS{dataset}.$FORM{MSA_Program}.semphy.tree_".($FORM{Bootstraps}-1).".rooted")
		{
			exit_on_error("sys_error","$VARS{BootStrap_Dir}tree_".($FORM{Bootstraps}-1)."/$VARS{dataset}.$FORM{MSA_Program}.semphy.tree_".($FORM{Bootstraps}-1).".rooted"." does not exist/empty\n");
		}
	}
	
	$VARS{MSA_Depth}=MSA_Depth("$VARS{WorkingDir}$VARS{Alignment_File}");
	# sample OP
	if (($FORM{MSA_Program} eq "MAFFT") and ($FORM{PROGRAM} eq "GUIDANCE2"))
	{
		$OP_DistFile=Guidance::MAFFT_OP_DIST;
	}
	elsif (($FORM{MSA_Program} eq "MAFFT") and ($FORM{PROGRAM} eq "GUIDANCE3_HOT"))
	{
		$OP_DistFile=Guidance::MAFFT_OP_DIST_0_25;
		$EP_DistFile=Guidance::MAFFT_EP_DIST_0_25;
	}
	my $OutEP="$VARS{WorkingDir}/SampledEPVals.log";
	my $OutOP="$VARS{WorkingDir}/SampledOPVals.log";
	my $op_vals_arr_ref="";
	my $ep_vals_arr_ref="";
	if ($FORM{PROGRAM} eq "GUIDANCE2")
	{
		if (uc ($VARS{GapPenDist}) eq "EMP")
		{
			print LOG "Sample op according to empiric distribution: Guidance::SampelFromEmpiricDistribution($OP_DistFile,$OutOP,$FORM{Bootstraps})\n";
			$op_vals_arr_ref=Guidance::SampelFromEmpiricDistribution($OP_DistFile,$OutOP,$FORM{Bootstraps});
		}
		elsif (uc ($VARS{GapPenDist}) eq "UNIF")
		{
			if ($FORM{MSA_Program} eq "MAFFT")
			{
				print LOG "Sample op according to uniform distribution: Guidance::SampleFromUnifomDist(1,3,$OutOP,$FORM{Bootstraps})\n";
				$op_vals_arr_ref=Guidance::SampleFromUnifomDist(1,3,$OutOP,$FORM{Bootstraps}); # according to mafft web-site defaults: http://mafft.cbrc.jp/alignment/server/index.html
			}
			elsif ($FORM{MSA_Program} eq "PRANK")
			{
#			print LOG "Sample op according to uniform distribution: Guidance::SampleFromUnifomDist(0,0.01,$OutOP,$FORM{Bootstraps})\n";
#			$op_vals_arr_ref=Guidance::SampleFromUnifomDist(0,0.01,$OutOP,$FORM{Bootstraps}); # for prank v.140110 the defaults are: dna 0.025 / prot 0.005
				print LOG "Sample op according to uniform distribution: Guidance::SampleFromUnifomDist(0,0.5,$OutOP,$FORM{Bootstraps})\n";
				$op_vals_arr_ref=Guidance::SampleFromUnifomDist(0,0.5,$OutOP,$FORM{Bootstraps}); # for prank v.140110 the defaults are: dna 0.025 / prot 0.005
			}
			elsif ($FORM{MSA_Program} eq "CLUSTALW")
			{
				print LOG "Sample gap opening panelty according to uniform distribution: Guidance::SampleFromUnifomDist(4,16,$OutOP,$FORM{Bootstraps})\n";
				$op_vals_arr_ref=Guidance::SampleFromUnifomDist(4,16,$OutOP,$FORM{Bootstraps});
			}
		}
	}
	if ($FORM{PROGRAM} eq "GUIDANCE3_HOT")
	{
		if (uc ($VARS{GapPenDist}) eq "EMP")
		{
			print LOG "Sample op according to empiric distribution: Guidance::SampelFromEmpiricDistribution($OP_DistFile,$OutOP,$FORM{Bootstraps})\n";
			$op_vals_arr_ref=Guidance::SampelFromEmpiricDistribution($OP_DistFile,$OutOP,$FORM{Bootstraps});
			print LOG "Sample ep according to empiric distribution: Guidance::SampelFromEmpiricDistribution($EP_DistFile,$OutEP,$FORM{Bootstraps})\n";
			$ep_vals_arr_ref=Guidance::SampelFromEmpiricDistribution($EP_DistFile,$OutEP,$FORM{Bootstraps});
		}
		if (uc ($VARS{GapPenDist}) eq "UNIF")
		{
			print LOG "Sample op according to uniform distribution: Guidance::SampleFromUnifomDist(0,6,$OutOP,$FORM{Bootstraps})\n";
			$op_vals_arr_ref=Guidance::SampleFromUnifomDist(0,6,$OutOP,$FORM{Bootstraps});
			print LOG "Sample ep according to uniform distribution: Guidance::SampleFromUnifomDist(0,4,$OutEP,$FORM{Bootstraps})\n";
			$ep_vals_arr_ref=Guidance::SampleFromUnifomDist(0,4,$OutEP,$FORM{Bootstraps});
		}
	}
	if ($isServer == 1) {
		$status_file="$VARS{WorkingDir}MSA_STATUS.txt";
		open (STATUS,">$status_file");
		print STATUS "<ul><li><p><font face=Verdana size=2>Start creating alternative alignments<br></li></ul>\n";
		close (STATUS);
		print OUTPUT "<?php\n\tif (file_exists('MSA_STATUS.txt'))\n\t{\n\t\t\$fil = fopen('MSA_STATUS.txt', r);\n\t\t\$dat = fread(\$fil, filesize('MSA_STATUS.txt'));\n\t\techo \"\$dat\";\n\tfclose(\$fil);\n\t}\n?>\n";
	}
	# Get random branches for HoT
	my $NumOfBrances=((2*$VARS{MSA_Depth})-3);
	my @RandomBranches=();
	for (my $i=0;$i<$FORM{Bootstraps};$i++)
	{
		$RandomBranches[$i]=int(rand($NumOfBrances));
	}
	# CREATE THE PERTUBED ALN DIR
	mkdir ($VARS{Scoring_Alignments_Dir});


		# Running parallel processes using fork, each will run an equal share of the BP alignments 
#	my $proc_num = 8;
	my $bp_per_proc = int($FORM{Bootstraps} / $VARS{proc_num}) + 1;
	my @children;  # array of pid of children
	for (my $proc=0; $proc<$VARS{proc_num}; $proc++) {
		print "Running proc num $proc\n";
		
		my $pid = fork();
		if ($pid) {
			
			# parent
			push(@children, $pid);
			
		} elsif ($pid == 0) {
			
			# child
			for (my $tree_num=0; $tree_num<$bp_per_proc; $tree_num++) {
				my $countTrees = $proc * $bp_per_proc + $tree_num;
				#for (my $i=0;$i<$FORM{Bootstraps};$i++) # NON PARALLEL FOR
				#{
				print "proc num $proc\ttree num $tree_num --> global tree index $countTrees\n";
				last if ($countTrees == $FORM{Bootstraps});
				#my $Branch=int(rand($NumOfBrances));
				my $Branch=$RandomBranches[$countTrees];
				my $tree="";
				if ($FORM{MSA_Program} eq "MAFFT")
				{
					if ($VARS{align_param}=~/addfragments/)
					{
						$tree="$VARS{prune_BootStrap_Dir}tree_".$countTrees."/$VARS{dataset}.$FORM{MSA_Program}.semphy.tree_".$countTrees."CORE.rooted";
					}
					else
					{
						$tree="$VARS{BootStrap_Dir}tree_".$countTrees."/$VARS{dataset}.$FORM{MSA_Program}.semphy.tree_".$countTrees.".rooted";
					}
				}
				else
				{
					$tree="$VARS{BootStrap_Dir}nonUniqueTrees/tree_".$countTrees."/$VARS{dataset}.$FORM{MSA_Program}.semphy.tree_".$countTrees;
				}
					if ($FORM{MSA_Program} eq "PRANK")
					{
						if ($VARS{align_param}=~/iterate/)
						{
							$VARS{align_param}=~s/\\\-iterate\=\d+//;
							print LOG "[WARNNING] -iterate argument is ignored for the pertubed alignments stage\n";
							print "[WARNNING] -iterate argument is ignored when reconstucting the pertubed alignments\n";
						}
						$VARS{align_param}=$VARS{align_param}." -once" if ($VARS{align_param}!~/once/);
					}
				my $tree_good_BranchLength=$tree.".GoodBranchLength";
				Guidance::reformat_trees_branch_length($tree,$tree_good_BranchLength);
				my $HOT_COS_GUIDANCE2_cmd="cd $VARS{WorkingDir}; perl ".Guidance::HOT_GUIDANCE2_PROGRAM." $VARS{dataset}_$countTrees $VARS{HoT_MSA_Program}";
				
				if (($FORM{Seq_Type} eq "AminoAcids") or ($FORM{Seq_Type} eq "Codons")) {$HOT_COS_GUIDANCE2_cmd.=" aa";}
				elsif ($FORM{Seq_Type} eq "Nucleotides") {$HOT_COS_GUIDANCE2_cmd.=" nt";}
					
				# $HOT_COS_GUIDANCE2_cmd.=" $VARS{codded_seq_fileName} . ../MSA_STATUS.txt 0 $VARS{HoT_MSA_Program_path} $tree_good_BranchLength $Branch";
				$HOT_COS_GUIDANCE2_cmd.=" $VARS{codded_seq_fileName} . \"\" 0 $VARS{HoT_MSA_Program_path} $tree_good_BranchLength $Branch";
				if ($FORM{MSA_Program} eq "MAFFT")
				{
					if ($FORM{PROGRAM} eq "GUIDANCE2")
					{
						if ($VARS{align_param} ne "") {$HOT_COS_GUIDANCE2_cmd.=" --- $VARS{align_param} --op $op_vals_arr_ref->[$countTrees]";}
						else {$HOT_COS_GUIDANCE2_cmd.=" --- --op $op_vals_arr_ref->[$countTrees]";}
					}
					elsif ($FORM{PROGRAM} eq "GUIDANCE3_HOT")
					{
						if ($VARS{align_param} ne "") {$HOT_COS_GUIDANCE2_cmd.=" --- $VARS{align_param} --op $op_vals_arr_ref->[$countTrees] --ep $ep_vals_arr_ref->[$countTrees]";}
						else {$HOT_COS_GUIDANCE2_cmd.=" --- --op $op_vals_arr_ref->[$countTrees] --ep $ep_vals_arr_ref->[$countTrees]";}
					}
				}
				elsif ($FORM{MSA_Program} eq "PRANK") # TO DO - add PRANK params
				{
					if ($VARS{align_param} ne "") {$HOT_COS_GUIDANCE2_cmd.=" --- $VARS{align_param} -gaprate=$op_vals_arr_ref->[$countTrees]";}
					else {$HOT_COS_GUIDANCE2_cmd.=" --- -gaprate=$op_vals_arr_ref->[$countTrees]";}
				}
				elsif ($FORM{MSA_Program} eq "CLUSTALW")
				{
					if ($VARS{align_param} ne "") {$HOT_COS_GUIDANCE2_cmd.=" --- $VARS{align_param} -GAPOPEN=$op_vals_arr_ref->[$countTrees]";}
					else {$HOT_COS_GUIDANCE2_cmd.=" --- -GAPOPEN=$op_vals_arr_ref->[$countTrees]";}
				}
				$HOT_COS_GUIDANCE2_cmd.=" >> COS.std";
				
				#print "$HOT_COS_GUIDANCE2_cmd\n";
				print LOG "run_HOT_COS_GUIDANCE2: $HOT_COS_GUIDANCE2_cmd\n";
				system ($HOT_COS_GUIDANCE2_cmd);
				
				#cp relevant MSAs files to $VARS{Scoring_Alignments_Dir}
				my @pertubed_aln=();
				
				if ($VARS{NumOfSeq}>2)
				{
					my $local_dataset="$VARS{dataset}_$countTrees"."_cos_";
					@pertubed_aln=<$VARS{WorkingDir}$local_dataset$VARS{HoT_MSA_Program}/b[01]*.fasta>;
				}
				else
				{
					my $local_dataset="$VARS{dataset}_$countTrees"."_cos_";
					@pertubed_aln=<$VARS{WorkingDir}$local_dataset$VARS{HoT_MSA_Program}/hot*.fasta>;
				}
				# print "ALNS:\n",join("\n",@pertubed_aln),"\n";<STDIN>; # QA
				# print "PERTUBED:@pertubed_aln\n";
				my @shuffled=shuffle(@pertubed_aln); 
				for (my $j=0;$j<$Num_of_Aln_from_HoT_per_Run;$j++) # TO DO: check what happen if 2 MSA only
				{
					my $aln=$shuffled[$j]; # random aln of the pertubed set
					my $base=getBasename($aln);
					# print "BASE:$base\n";<STDIN>;
					if ($FORM{PROGRAM} eq "GUIDANCE2")
					{
						$cmd="cp $VARS{WorkingDir}$VARS{dataset}_$countTrees"."_cos_$VARS{HoT_MSA_Program}/$base".".fasta $VARS{Scoring_Alignments_Dir}$base"."_tree_".$countTrees."_OP_".$op_vals_arr_ref->[$countTrees]."_Split_".$Branch.".fasta";
					}
					elsif ($FORM{PROGRAM} eq "GUIDANCE3_HOT")
					{
						$cmd="cp $VARS{WorkingDir}$VARS{dataset}_$countTrees"."_cos_$VARS{HoT_MSA_Program}/$base".".fasta $VARS{Scoring_Alignments_Dir}$base"."_tree_".$countTrees."_OP_".$op_vals_arr_ref->[$countTrees]."_EP_".$ep_vals_arr_ref->[$countTrees]."_Split_".$Branch.".fasta";
					}
					# print $cmd;
					# <STDIN>;
					system ($cmd);
				}
				#move out dir to BP dir and rename
				my $cmd="";
				if ($FORM{PROGRAM} eq "GUIDANCE2")
				{
					$cmd="mv $VARS{WorkingDir}$VARS{dataset}_$countTrees"."_cos_"."$VARS{HoT_MSA_Program} $VARS{BootStrap_Dir}$VARS{dataset}"."_cos_".$VARS{HoT_MSA_Program}."_tree_".$countTrees."_OP_".$op_vals_arr_ref->[$countTrees]."_Split_".$Branch;
				}
				elsif ($FORM{PROGRAM} eq "GUIDANCE3_HOT")
				{
					$cmd="mv $VARS{WorkingDir}$VARS{dataset}_$countTrees"."_cos_"."$VARS{HoT_MSA_Program} $VARS{BootStrap_Dir}$VARS{dataset}"."_cos_".$VARS{HoT_MSA_Program}."_tree_".$countTrees."_OP_".$op_vals_arr_ref->[$countTrees]."EP_".$ep_vals_arr_ref->[$countTrees]."_Split_".$Branch;
				}
				# print $cmd;<STDIN>;
				## paralleization
				
				
				system ($cmd);
				if ($isServer==1)
				{
					my $i=$tree_num/$VARS{proc_num};
					open (PROGRESS,">$status_file");
					print PROGRESS "\n<ul><li>",($i+1)*4," out of ",$FORM{Bootstraps}*4," alternative alignments were created</li></ul>\n";
					close (PROGRESS);
				}
			}
			exit 0;
		} else {
			die "ERROR: fork failed: $!\n";
		}
		
	}
	
	# Wait for child processes to end
	foreach (@children) {
        my $pid = waitpid($_, 0);
		print "done with pid $pid\n";
	}
	#print "FINISHED ALL ALNS...\n";#<STDIN>;
	# validate all alns were created
	my $aln_count=count_files_on_dir($VARS{Scoring_Alignments_Dir});
	if ($aln_count<$Num_of_Aln_from_HoT_per_Run*$FORM{Bootstraps})
	{
		exit_on_error("sys_error","run_Guidance2: Only $aln_count alignments were created on $VARS{Scoring_Alignments_Dir} while expecting ".$Num_of_Aln_from_HoT_per_Run*$FORM{Bootstraps}."\n");
	}
}
sub run_HoT
{
#	perl -w /groups/pupko/osnatz/projects/alignment/bootstrap/Giddy/cos/cos_code_example/COS.pl caseID msa_method seq_type input_fasta_file . output_dir >& COS.std
#msa_method: MA0 = mafft ; CW2 = clustalW.
#seq_type: aa = amino-acid ; nt = nucleotides
#input_fasta_file = the input sequences file

#The base MSA is: output_dir_cos_msa_method /hot_H.fasta
#So you can copy it:
#cp output_dir_cos_msa_method/hot_H.fasta ./caseID_mafft.fasta

#The MSA sets should be copied to a directory:
#mkdir COS_MSA
#cp ./output_dir_cos_msa_method/b0#*.fasta ./COS_MSA/
	
	
	if (($FORM{MSA_Program} eq "MAFFT") or ($FORM{MSA_Program} eq "MAFFT_LINSI")) 
	{
		$VARS{HoT_MSA_Program}="MFT";
		$VARS{HoT_MSA_Program_path}=$mafft_prog;
		my $check_mafft_profile=`which $VARS{HoT_MSA_Program_path}-profile`;
		if ($check_mafft_profile=~/Command not found/){die "It seems that $VARS{HoT_MSA_Program_path}-profile is not properly installed or found in PATH variable. Please fix that and/or provide GUIDANCE with the full path to mafft installation using the --mafft argument\n";}
	}
	if ($FORM{MSA_Program} eq "MUSCLE") {exit_on_error('user_error', "HoT currently does not support MUSCLE, please run GUIDANCE<br>");}
	if ($FORM{MSA_Program} eq "PAGAN") {exit_on_error('user_error', "HoT currently does not support PAGAN, please run GUIDANCE<br>");}
	elsif ($FORM{MSA_Program} eq "CLUSTALW") {$VARS{HoT_MSA_Program}="CLW";$VARS{HoT_MSA_Program_path}=$clustalw_prog;}
    elsif ($FORM{MSA_Program} eq "PRANK") {$VARS{HoT_MSA_Program}="PRK";$VARS{HoT_MSA_Program_path}=$prank_prog;}

	my $cmd="cd $VARS{WorkingDir}; perl ".Guidance::HOT_PROGRAM." $VARS{dataset} $VARS{HoT_MSA_Program}";
	
	if (($FORM{Seq_Type} eq "AminoAcids") or ($FORM{Seq_Type} eq "Codons")) {$cmd.=" aa";}
	elsif ($FORM{Seq_Type} eq "Nucleotides") {$cmd.=" nt";}
	
#	$cmd.=" $VARS{codded_seq_fileName} $VARS{WorkingDir} $VARS{WorkingDir}MSA_STATUS.txt > $VARS{WorkingDir}COS.std";# 2>&1";
	print LOG "convert_fs_to_upper_case($VARS{WorkingDir}$VARS{codded_seq_fileName})\n";
	convert_fs_to_upper_case("$VARS{WorkingDir}$VARS{codded_seq_fileName}"); #HoT assumes that all the sequences are upper case
	if ($VARS{align_param} eq ""){$cmd.=" $VARS{codded_seq_fileName} . ../MSA_STATUS.txt 0 $VARS{HoT_MSA_Program_path} > COS.std";}
	else {$cmd.=" $VARS{codded_seq_fileName} . ../MSA_STATUS.txt 0 $VARS{HoT_MSA_Program_path} --- $VARS{align_param} > COS.std";}
	print LOG "run_HoT: $cmd\n";
	if ($isServer == 1) {
		my $status_file="$VARS{WorkingDir}MSA_STATUS.txt";
		open (STATUS,">$status_file");
		print STATUS "\n";
		close (STATUS);
		open (STATUS0,">$status_file.0");
		print STATUS0 "\n";
		close (STATUS0);
		print OUTPUT "<?php\n\tif (file_exists('MSA_STATUS.txt.0'))\n\t{\n\t\t\$fil =fopen('MSA_STATUS.txt.0', r);\n\t\t\$dat = fread(\$fil, filesize('MSA_STATUS.txt.0'));\n\t\techo \"\$dat\";\n\tfclose(\$fil);\n\t}\n?>\n";
		
		print OUTPUT "<?php\n\tif (file_exists('MSA_STATUS.txt'))\n\t{\n\t\t\$fil = fopen('MSA_STATUS.txt', r);\n\t\t\$dat = fread(\$fil, filesize('MSA_STATUS.txt'));\n\t\techo \"\$dat\";\n\tfclose(\$fil);\n\t}\n?>\n";
	}
	system ($cmd);
	
	if (!-e "$VARS{WorkingDir}$VARS{Alignment_File}")
	{
		$cmd="cp $VARS{WorkingDir}$VARS{dataset}_cos_$VARS{HoT_MSA_Program}/hot_H.fasta $VARS{WorkingDir}$VARS{Alignment_File}";
		print LOG "run_HoT: $cmd\n";  # Copy Alignment
		system ($cmd);
	}
	if ($FORM{Align_Order} eq "as_input")
	{
		print "MSA_parser::sort_alignment($VARS{WorkingDir}$VARS{Alignment_File},fasta);\n";
		my @ans=MSA_parser::sort_alignment("$VARS{WorkingDir}$VARS{Alignment_File}","fasta");
		print join("",@ans);
		$VARS{Alignment_File_NOT_SORTED}=$VARS{Alignment_File};
		$VARS{Alignment_File}=$VARS{Alignment_File}.".Sorted";
	}
	
	system ("mkdir $VARS{WorkingDir}$VARS{HoT_MSAs_Dir}");
	if ($VARS{NumOfSeq}>2)
	{
		$cmd="cp $VARS{WorkingDir}$VARS{dataset}_cos_$VARS{HoT_MSA_Program}/b[01]*.fasta $VARS{WorkingDir}$VARS{HoT_MSAs_Dir}";
	}
	else
	{
		$cmd="cp $VARS{WorkingDir}$VARS{dataset}_cos_$VARS{HoT_MSA_Program}/hot*.fasta $VARS{WorkingDir}$VARS{HoT_MSAs_Dir}";
	}
	print LOG "run_HoT: $cmd\n";
	system ($cmd);
	convert_fs_to_upper_case("$VARS{WorkingDir}$VARS{Alignment_File}.WithCodesName") if (-e "$VARS{WorkingDir}$VARS{Alignment_File}.WithCodesName"); #HoT assumes that all the sequences are upper case
	$VARS{Scoring_Alignments_Dir}="$VARS{WorkingDir}$VARS{HoT_MSAs_Dir}";
	
}
sub run_Guidance
{	
	# Align
	##############
	Align() if ($FORM{userMSA_File} eq ""); #align if user did not supplied alignment

	# handle the adjustdirection
	if (($VARS{align_param}=~/adjustdirection/) or ($VARS{align_param}=~/adjustdirectionaccurately/)) #support the --adjustdirectionaccurately and --adjustdirection option in MAFFT
	{
		if (($FORM{MSA_Program} eq "MAFFT") and ($FORM{Seq_Type} eq "Nucleotides") and ($FORM{userMSA_File} eq "")) # effective only when MAFFT is used, nucleotides sequence and user provided sequence
		{
			# back the original seq provided
			
			copy ("$VARS{WorkingDir}$VARS{codded_seq_fileName}","$VARS{WorkingDir}$VARS{codded_seq_fileName}.OrigDirection");
			copy("$VARS{WorkingDir}$VARS{Alignment_File}","$VARS{WorkingDir}$VARS{Alignment_File}.OrigMAFFT");
			# After the base alignment take all the sequences out of the alignment and continue with them
			print LOG "adjustdirection type option is in use; going to:\n";
			print LOG "\t\t(1) copy $VARS{WorkingDir}$VARS{codded_seq_fileName} to $VARS{WorkingDir}$VARS{codded_seq_fileName}.OrigDirection\n";
			print LOG "\t\t(2) copy $VARS{WorkingDir}$VARS{Alignment_File} to $VARS{WorkingDir}$VARS{Alignment_File}.OrigMAFFT\n";
			print LOG "\t\t(3) Fix $VARS{WorkingDir}$VARS{codded_seq_fileName} to contain the sequences in the correct direction\n";
			print LOG "\t\t(4) Fix $VARS{WorkingDir}$VARS{Alignment_File} to remove the '_R_' prefix added by MAFFT to the filipped seq names\n";
			
			open (my $MSA,"<","$VARS{WorkingDir}$VARS{Alignment_File}.OrigMAFFT") || exit_on_error ("sys_error","Can't open '$VARS{WorkingDir}$VARS{Alignment_File}' for reading $!");
			open (my $SEQ_FILE,">","$VARS{WorkingDir}$VARS{codded_seq_fileName}") || exit_on_error ("sys_error","Can't open '$VARS{WorkingDir}$VARS{codded_seq_fileName}' for writing $!");
			open (my $MSA_NEW,">","$VARS{WorkingDir}$VARS{Alignment_File}") || exit_on_error ("sys_error","Can't open '$VARS{WorkingDir}$VARS{Alignment_File}' for writing $!");
			my ($FRAGMENT_FILE,$NEW_CORE_ALN);
			if ($VARS{align_param}=~/addfragments/) # alignment with addfragments
			{
				print LOG "\t\t(5) Copy original core alignment $VARS{WorkingDir}$VARS{Core_Alignment_File} to $VARS{WorkingDir}$VARS{Core_Alignment_File}.OrigMAFFT\n";
				copy ("$VARS{WorkingDir}$VARS{Core_Alignment_File}","$VARS{WorkingDir}$VARS{Core_Alignment_File}.OrigMAFFT");
				print LOG "\t\t(6) Copy $VARS{WorkingDir}$VARS{fragments_file_name_seqName_coded} to $VARS{WorkingDir}$VARS{fragments_file_name_seqName_coded}.OrigDirection\n";
				copy ("$VARS{WorkingDir}$VARS{fragments_file_name_seqName_coded}","$VARS{WorkingDir}$VARS{fragments_file_name_seqName_coded}.OrigDirection");
				print LOG "\t\t(7) Fix $VARS{WorkingDir}$VARS{fragments_file_name_seqName_coded} to contain the fragments in the correct direction\n";
				print LOG "\t\t(8) Fix $VARS{WorkingDir}$VARS{Core_Alignment_File} to remove the '_R_' prefix added by MAFFT to the filipped seq names\n";
				open ($FRAGMENT_FILE,">","$VARS{WorkingDir}$VARS{fragments_file_name_seqName_coded}") || exit_on_error ("sys_error","Can't open '$VARS{WorkingDir}$VARS{fragments_file_name_seqName_coded}' for writing $!");
				open (my $CORE_ALN,"<","$VARS{WorkingDir}$VARS{Core_Alignment_File}.OrigMAFFT") || exit_on_error ("sys_error","Can't open '$VARS{WorkingDir}$VARS{Core_Alignment_File}.OrigMAFFT' for writing $!");
				open ($NEW_CORE_ALN,">","$VARS{WorkingDir}$VARS{Core_Alignment_File}") || exit_on_error ("sys_error","Can't open '$VARS{WorkingDir}$VARS{Core_Alignment_File}' for writing $!");
				print LOG "TOTAL SEQ IN CORE ALN: $VARS{NumOfSeq}\n";
				while (my $line=<$CORE_ALN>)
				{
					chomp ($line);
					if ($line=~/^>_R_(.*)/)
					{
						print $NEW_CORE_ALN ">$1\n";
						print LOG "[NOTICE] MAFFT REVERSE SEQUENCE $1 in the CORE ALIGNMENT\n";
					}
					elsif ($line=~/^>/)
					{
						print $NEW_CORE_ALN "$line\n";
					}
					else
					{
						print $NEW_CORE_ALN "$line\n";
					}
				}
				close ($CORE_ALN);
				close ($NEW_CORE_ALN);
			}
				
			my $IsReversedSeq=0;
			my $i=0;
			while (my $line=<$MSA>)
			{
				chomp ($line);
				if ($line=~/^>_R_(.*)/)
				{
					my $reverse_seq_name=$1;
					$i++;
					$IsReversedSeq=1;
					if ($i<=$VARS{NumOfSeq})
					{
						print $SEQ_FILE ">$reverse_seq_name\n";
						print LOG "[NOTICE] MAFFT REVERSE SEQUENCE $reverse_seq_name\n";
					}
					elsif ($VARS{align_param}=~/addfragments/)
					{
						print $FRAGMENT_FILE ">$reverse_seq_name\n";
						print LOG "[NOTICE] MAFFT REVERSE FRAGMENT SEQUENCE $reverse_seq_name\n";
					}
					print $MSA_NEW ">$reverse_seq_name\n";
				}
				elsif ($line=~/^>/)
				{
					$i++;
					if ($i<=$VARS{NumOfSeq})
					{
						print $SEQ_FILE "$line\n";
					}
					elsif ($VARS{align_param}=~/addfragments/) # fragment
					{
						print $FRAGMENT_FILE "$line\n";
					}
					print $MSA_NEW "$line\n";
				}
				else
				{
					print $MSA_NEW "$line\n";
					$line=~s/-//g;
					if ($line ne "")
					{
						if ($i<=$VARS{NumOfSeq})
						{
							print $SEQ_FILE "$line\n";
						}
						elsif ($VARS{align_param}=~/addfragments/) # fragment
						{
							print $FRAGMENT_FILE "$line\n";
						}
					}
				}
			}
			close ($MSA);
			close ($SEQ_FILE);
			close ($MSA_NEW);
			close ($FRAGMENT_FILE) if ($VARS{align_param}=~/addfragments/);
			if ($IsReversedSeq==0)
			{
				copy ("$VARS{WorkingDir}$VARS{Alignment_File}.OrigMAFFT","$VARS{WorkingDir}$VARS{Alignment_File}");
				copy ("$VARS{WorkingDir}$VARS{codded_seq_fileName}.OrigDirection","$VARS{WorkingDir}$VARS{codded_seq_fileName}");
				unlink ("$VARS{WorkingDir}$VARS{codded_seq_fileName}.OrigDirection");
				unlink ("$VARS{WorkingDir}$VARS{Alignment_File}.OrigMAFFT");
				if ($VARS{align_param}=~/\-\-addfragments/)
				{
						copy ("$VARS{WorkingDir}$VARS{Core_Alignment_File}.OrigMAFFT","$VARS{WorkingDir}$VARS{Core_Alignment_File}");
						copy ("$VARS{WorkingDir}$VARS{fragments_file_name_seqName_coded}.OrigDirection","$VARS{WorkingDir}$VARS{fragments_file_name_seqName_coded}");
						unlink ("$VARS{WorkingDir}$VARS{fragments_file_name_seqName_coded}.OrigDirection");
						unlink ("$VARS{WorkingDir}$VARS{Core_Alignment_File}.OrigMAFFT");
				}
				print LOG "No fillped sequences were found... NOTHING OF THE ABOVEMENTIONED WAS NEEDED!!\n";
			}
			
		}

		my @ans=Guidance::ConvertNamesOfAlignWithSeed("$VARS{WorkingDir}$VARS{Alignment_File}","$VARS{WorkingDir}$VARS{Alignment_File}"."_new");
		if ($ans[0] ne "ok"){exit_on_error("sys_error","ConvertNamesOfAlignWithSeed("."$VARS{WorkingDir}$VARS{Alignment_File}","$VARS{WorkingDir}$VARS{Alignment_File}"."_new"."): ".join(" ",@ans)."\n");}
		else {rename("$VARS{WorkingDir}$VARS{Alignment_File}"."_new","$VARS{WorkingDir}$VARS{Alignment_File}");}
	}
	
	# BootStrap Trees
	##################
	Bootstrap_Trees();

	# pull out the trees
	######################
	my $numUniqueTrees="";
	my @numRepeats4UniqueTree=();
	if (uc $VARS{BBL} eq "NO")
	{
		print LOG "Guidance::pullOutBPtrees($VARS{WorkingDir},$VARS{dataset},$FORM{Bootstraps},$FORM{MSA_Program});\n";
		@ans=Guidance::pullOutBPtrees($VARS{WorkingDir},$VARS{dataset},$FORM{Bootstraps},$FORM{MSA_Program});
		unless ($ans[0] eq "ok") {exit_on_error("sys_error","Guidance::pullOutBPtrees: ".join (" ",@ans)."\n");}
		$numUniqueTrees=$ans[1] if ($FORM{MSA_Program} ne "MAFFT");
		@numRepeats4UniqueTree=@{$ans[2]} if ($FORM{MSA_Program} ne "MAFFT");
	}
	else
	{
		print LOG "Guidance::pullOutBPtrees_BBL($VARS{WorkingDir},$VARS{dataset},$FORM{Bootstraps},$FORM{MSA_Program});\n";
		@ans=Guidance::pullOutBPtrees_BBL($VARS{WorkingDir},$VARS{dataset},$FORM{Bootstraps},$FORM{MSA_Program});
		unless ($ans[0] eq "ok") {exit_on_error("sys_error","Guidance::pullOutBPtrees_BBL: ".join (" ",@ans)."\n");}
		$numUniqueTrees=$ans[1] if ($FORM{MSA_Program} ne "MAFFT");
		@numRepeats4UniqueTree=@{$ans[2]} if ($FORM{MSA_Program} ne "MAFFT");
	}
	# Convert trees to MAFFT format
	#################################
	if (($FORM{MSA_Program} eq "MAFFT") or ($FORM{MSA_Program} eq "MAFFT_LINSI")) #FOR MAFFT BUILED ALIGNMENT FO ALL TREES
	{
		# if addfragments is used first prune the trees for the the 'core' alignment and afterwards will use the created 'BP core' and the 'full tree' to create the full alignment
		if ($VARS{align_param}=~/addfragments/)
		{
			$VARS{prune_BootStrap_Dir}=$VARS{BootStrap_Dir}."PRUNE_BP_FOR_CORE_ALN/";
			# prune trees
			unless (-e $VARS{prune_BootStrap_Dir}){system ("mkdir $VARS{prune_BootStrap_Dir}");}
			foreach my $full_tree_dir (<$VARS{BootStrap_Dir}tree*>)
			{
				#print $full_tree_dir;<STDIN>;
				if (-d $full_tree_dir)
				{
					if ($full_tree_dir=~/(tree_[0-9]+)/)
					{
						my $tree_num=$1;
						unless (-e "$VARS{prune_BootStrap_Dir}$tree_num"){system ("mkdir $VARS{prune_BootStrap_Dir}$tree_num");}
						my $FullTree="$full_tree_dir/$VARS{dataset}.$FORM{MSA_Program}.semphy.$tree_num";
						my $PrunedTree="$VARS{prune_BootStrap_Dir}/$tree_num/$VARS{dataset}.$FORM{MSA_Program}.semphy.$tree_num"."CORE";
						system ("$remove_taxa_prog $FullTree $VARS{WorkingDir}$VARS{fragments_codes} $PrunedTree");
					}
				}
			}
			# convert pruned tree to MAFFT format
			print LOG "Guidance::convertBPTrees2MafftFormat($VARS{prune_BootStrap_Dir},$VARS{dataset}, $FORM{MSA_Program}, $FORM{Bootstraps},$ruby_prog,CORE,$VARS{rooting_type});\n";
			@ans=Guidance::convertBPTrees2MafftFormat($VARS{prune_BootStrap_Dir},$VARS{dataset}, $FORM{MSA_Program}, $FORM{Bootstraps},$ruby_prog,"CORE",$VARS{rooting_type});
			unless ($ans[0] eq "ok"){exit_on_error("sys_error","Guidance::convertBPTrees2MafftFormat: ".join(" ",@ans)."\n");}
			if (-z "$VARS{prune_BootStrap_Dir}tree_".($FORM{Bootstraps}-1)."/$VARS{dataset}.$FORM{MSA_Program}.semphy.tree_".($FORM{Bootstraps}-1)."CORE.mafftFormat" or !-e "$VARS{prune_BootStrap_Dir}tree_".($FORM{Bootstraps}-1)."/$VARS{dataset}.$FORM{MSA_Program}.semphy.tree_".($FORM{Bootstraps}-1)."CORE.mafftFormat")
			{
				exit_on_error("sys_error","$VARS{BootStrap_Dir}tree_".($FORM{Bootstraps}-1)."/$VARS{dataset}.$FORM{MSA_Program}.semphy.tree_".($FORM{Bootstraps}-1)."CORE.mafftFormat"." does not exist/empty\n");
			}
		}
		# prepare the trees
		print LOG "Guidance::convertBPTrees2MafftFormat($VARS{BootStrap_Dir},$VARS{dataset}, $FORM{MSA_Program}, $FORM{Bootstraps},$ruby_prog,\"\",$VARS{rooting_type});\n";
		@ans=Guidance::convertBPTrees2MafftFormat($VARS{BootStrap_Dir},$VARS{dataset}, $FORM{MSA_Program}, $FORM{Bootstraps},$ruby_prog,"",$VARS{rooting_type});
		unless ($ans[0] eq "ok"){exit_on_error("sys_error","Guidance::convertBPTrees2MafftFormat: ".join(" ",@ans)."\n");}
		if (-z "$VARS{BootStrap_Dir}tree_".($FORM{Bootstraps}-1)."/$VARS{dataset}.$FORM{MSA_Program}.semphy.tree_".($FORM{Bootstraps}-1).".mafftFormat" or !-e "$VARS{BootStrap_Dir}tree_".($FORM{Bootstraps}-1)."/$VARS{dataset}.$FORM{MSA_Program}.semphy.tree_".($FORM{Bootstraps}-1).".mafftFormat")
		{
			exit_on_error("sys_error","$VARS{BootStrap_Dir}tree_".($FORM{Bootstraps}-1)."/$VARS{dataset}.$FORM{MSA_Program}.semphy.tree_".($FORM{Bootstraps}-1).".mafftFormat"." does not exist/empty\n");
		}	
	}
	# Root trees for MUSCLE and PAGAN
	if (($FORM{MSA_Program} eq "MUSCLE") or ($FORM{MSA_Program} eq "PAGAN"))
	{
		print LOG "Guidance::root_BP_trees($VARS{BootStrap_Dir},$VARS{dataset}, $FORM{MSA_Program}, $numUniqueTrees,\"\",$VARS{rooting_type});\n";
		@ans=Guidance::root_BP_trees($VARS{BootStrap_Dir},$VARS{dataset}, $FORM{MSA_Program}, $numUniqueTrees,"",$VARS{rooting_type});
		unless ($ans[0] eq "ok"){exit_on_error("sys_error","Guidance::root_BP_trees: ".join(" ",@ans)."\n");}
		if (-z "$VARS{BootStrap_Dir}tree_".($numUniqueTrees-1)."/$VARS{dataset}.$FORM{MSA_Program}.semphy.tree_".($numUniqueTrees-1).".rooted" or !-e "$VARS{BootStrap_Dir}tree_".($numUniqueTrees-1)."/$VARS{dataset}.$FORM{MSA_Program}.semphy.tree_".($numUniqueTrees-1).".rooted")
		{
			exit_on_error("sys_error","$VARS{BootStrap_Dir}tree_".($numUniqueTrees-1)."/$VARS{dataset}.$FORM{MSA_Program}.semphy.tree_".($numUniqueTrees-1).".rooted"." does not exist/empty\n");
		}	
	}
	
	# runAlignBPtrees
	##################################
	my $status_file = "";
	if ($isServer == 1) {
		$status_file="$VARS{WorkingDir}MSA_STATUS.txt";
		open (STATUS,">$status_file");
		print STATUS "<ul><li><p><font face=Verdana size=2>Start creating alternative alignments<br></li></ul>\n";
		close (STATUS);
		print OUTPUT "<?php\n\tif (file_exists('MSA_STATUS.txt'))\n\t{\n\t\t\$fil = fopen('MSA_STATUS.txt', r);\n\t\t\$dat = fread(\$fil, filesize('MSA_STATUS.txt'));\n\t\techo \"\$dat\";\n\tfclose(\$fil);\n\t}\n?>\n";
	}
	
	if ($FORM{Seq_Type} eq "Codons")
	{
		if ($FORM{MSA_Program} eq "MAFFT") # FOR MAFFT BUILED ALIGNMENT FOR ALL TREES
		{
			if ($FORM{PROGRAM} eq "GUIDANCE")
			{
				print LOG "Guidance::runAlignBPtrees ($FORM{MSA_Program},$VARS{dataset},$VARS{WorkingDir},$VARS{WorkingDir}$VARS{codded_seq_fileName},AminoAcids,$FORM{MSA_Program},$FORM{Bootstraps},$VARS{align_param},\"\",$status_file,$mafft_prog,$prank_prog,$clustalw_prog,$muscle_prog,$pagan_prog,$VARS{proc_num});\n";
				@ans=Guidance::runAlignBPtrees ($FORM{MSA_Program},$VARS{dataset},$VARS{WorkingDir},"$VARS{WorkingDir}$VARS{codded_seq_fileName}","AminoAcids",$FORM{MSA_Program},$FORM{Bootstraps},$VARS{align_param},"",$status_file,$mafft_prog,$prank_prog,$clustalw_prog,$muscle_prog,$pagan_prog,$VARS{proc_num});
				unless ($ans[0] eq "ok"){exit_on_error("sys_error","Guidance::runAlignBPtrees: ".join(" ",@ans)."\n");}
			}
			elsif ($FORM{PROGRAM} eq "GUIDANCE3")
			{
				print LOG "Guidance::runAlignBPtrees_GUIDANCE3 ($FORM{MSA_Program},$VARS{dataset},$VARS{WorkingDir},$VARS{WorkingDir}$VARS{codded_seq_fileName},AminoAcids,$FORM{MSA_Program},$FORM{Bootstraps},$VARS{align_param},\"\",$status_file,$mafft_prog,$prank_prog,$clustalw_prog,$muscle_prog,$pagan_prog,$VARS{proc_num},$VARS{GapPenDist});\n";
				@ans=Guidance::runAlignBPtrees_GUIDANCE3 ($FORM{MSA_Program},$VARS{dataset},$VARS{WorkingDir},"$VARS{WorkingDir}$VARS{codded_seq_fileName}","AminoAcids",$FORM{MSA_Program},$FORM{Bootstraps},$VARS{align_param},"",$status_file,$mafft_prog,$prank_prog,$clustalw_prog,$muscle_prog,$pagan_prog,$VARS{proc_num},$VARS{GapPenDist});
				unless ($ans[0] eq "ok"){exit_on_error("sys_error","Guidance::runAlignBPtrees: ".join(" ",@ans)."\n");}
			}
		}
		else # BUILED ALIGNMENT ONLY FOR UNIQ TREES
		{
			print LOG "Guidance::runAlignBPtrees ($FORM{MSA_Program},$VARS{dataset},$VARS{WorkingDir},$VARS{WorkingDir}$VARS{codded_seq_fileName},AminoAcids,$FORM{MSA_Program},$numUniqueTrees,$VARS{align_param},\"\",$status_file,$mafft_prog,$prank_prog,$clustalw_prog,$muscle_prog,$pagan_prog,$VARS{proc_num});\n";
			@ans=Guidance::runAlignBPtrees ($FORM{MSA_Program},$VARS{dataset},$VARS{WorkingDir},"$VARS{WorkingDir}$VARS{codded_seq_fileName}","AminoAcids",$FORM{MSA_Program},$numUniqueTrees,$VARS{align_param},"",$status_file,$mafft_prog,$prank_prog,$clustalw_prog,$muscle_prog,$pagan_prog,$VARS{proc_num});
			unless ($ans[0] eq "ok"){exit_on_error("sys_error","Guidance::runAlignBPtrees: ".join(" ",@ans)."\n");}
		}
	}
	else # NOT CODONS
	{
		if ($FORM{MSA_Program} eq "MAFFT") # FOR MAFFT BUILED ALIGNMENT FOR ALL TREES
		{
			if ($FORM{PROGRAM} eq "GUIDANCE")
			{
				print LOG "Guidance::runAlignBPtrees ($FORM{MSA_Program},$VARS{dataset},$VARS{WorkingDir},$VARS{WorkingDir}$VARS{codded_seq_fileName},$FORM{Seq_Type},$FORM{MSA_Program},$FORM{Bootstraps},$VARS{align_param},\"\",$status_file,$mafft_prog,$prank_prog,$clustalw_prog,$muscle_prog,$pagan_prog,$VARS{proc_num});\n";
				@ans=Guidance::runAlignBPtrees 	($FORM{MSA_Program},$VARS{dataset},$VARS{WorkingDir},"$VARS{WorkingDir}$VARS{codded_seq_fileName}",$FORM{Seq_Type},$FORM{MSA_Program},$FORM{Bootstraps},$VARS{align_param},"",$status_file,$mafft_prog,$prank_prog,$clustalw_prog,$muscle_prog,$pagan_prog,$VARS{proc_num});
				unless ($ans[0] eq "ok"){exit_on_error("sys_error","Guidance::runAlignBPtrees: ".join(" ",@ans)."\n");}
			}
			elsif ($FORM{PROGRAM} eq "GUIDANCE3")
			{
				print LOG "Guidance::runAlignBPtrees_GUIDANCE3 ($FORM{MSA_Program},$VARS{dataset},$VARS{WorkingDir},$VARS{WorkingDir}$VARS{codded_seq_fileName},$FORM{Seq_Type},$FORM{MSA_Program},$FORM{Bootstraps},$VARS{align_param},\"\",$status_file,$mafft_prog,$prank_prog,$clustalw_prog,$muscle_prog,$pagan_prog,$VARS{proc_num},$VARS{GapPenDist});\n";
				@ans=Guidance::runAlignBPtrees_GUIDANCE3 ($FORM{MSA_Program},$VARS{dataset},$VARS{WorkingDir},"$VARS{WorkingDir}$VARS{codded_seq_fileName}",$FORM{Seq_Type},$FORM{MSA_Program},$FORM{Bootstraps},$VARS{align_param},"",$status_file,$mafft_prog,$prank_prog,$clustalw_prog,$muscle_prog,$pagan_prog,$VARS{proc_num},$VARS{GapPenDist});
				unless ($ans[0] eq "ok"){exit_on_error("sys_error","Guidance::runAlignBPtrees: ".join(" ",@ans)."\n");}
			}
			
		}
		else # BUILED ALIGNMENT ONLY FOR UNIQ TREES
		{
			print LOG "Guidance::runAlignBPtrees ($FORM{MSA_Program},$VARS{dataset},$VARS{WorkingDir},$VARS{WorkingDir}$VARS{codded_seq_fileName},$FORM{Seq_Type},$FORM{MSA_Program},$numUniqueTrees,$VARS{align_param},\"\",$status_file,$mafft_prog,$prank_prog,$clustalw_prog,$muscle_prog,$pagan_prog,$VARS{proc_num});\n";
			@ans=Guidance::runAlignBPtrees ($FORM{MSA_Program},$VARS{dataset},$VARS{WorkingDir},"$VARS{WorkingDir}$VARS{codded_seq_fileName}",$FORM{Seq_Type},$FORM{MSA_Program},$numUniqueTrees,$VARS{align_param},"",$status_file,$mafft_prog,$prank_prog,$clustalw_prog,$muscle_prog,$pagan_prog,$VARS{proc_num});
			unless ($ans[0] eq "ok"){exit_on_error("sys_error","Guidance::runAlignBPtrees: ".join(" ",@ans)."\n");}
		}
	}
    # copy all BP MSA into one directory
	#######################################
	if ($FORM{MSA_Program} eq "MAFFT") # FOR MAFFT BUILED ALIGNMENT FOR ALL TREES
	{
		mkdir ($VARS{BootStrap_MSA_Dir});
		system ("cp $VARS{BootStrap_Dir}/tree_*/*.$FORM{MSA_Program}.aln $VARS{BootStrap_MSA_Dir}");
	}
	else 
	{	
		print LOG "Guidance::copyBootstrapMSA2oneDir($VARS{BootStrap_MSA_Dir},$VARS{BootStrap_Dir},$FORM{MSA_Program},\@numRepeats4UniqueTree);\n";
		@ans = Guidance::copyBootstrapMSA2oneDir($VARS{BootStrap_MSA_Dir},$VARS{BootStrap_Dir},$FORM{MSA_Program},\@numRepeats4UniqueTree);
		unless ($ans[0] eq "ok"){exit_on_error("sys_error","Guidance::copyBootstrapMSA2oneDir: ".join(" ",@ans)."\n");}
	}
	$VARS{Scoring_Alignments_Dir}=$VARS{BootStrap_MSA_Dir};
	if ($VARS{align_param}=~/\-\-seed/) # alignment with seed, the names should be converted to consider the names starting with _seed_
	{
		
		foreach my $aln_file (<$VARS{BootStrap_MSA_Dir}/*.aln>)
		{
			my @ans=Guidance::ConvertNamesOfAlignWithSeed($aln_file,$aln_file."_new");
			if ($ans[0] ne "ok"){exit_on_error("sys_error","ConvertNamesOfAlignWithSeed($aln_file,$aln_file"."_new"."): ".join(" ",@ans)."\n");}
			else {rename($aln_file."_new", $aln_file);}
		}
	}
	my $aln_count=count_files_on_dir($VARS{BootStrap_MSA_Dir});
	if ($aln_count<$FORM{Bootstraps})
	{
		exit_on_error("sys_error","run_Guidance: Only $aln_count alignment were created on $VARS{BootStrap_MSA_Dir} while expecting $FORM{Bootstraps}\n");
	}

}
#---------------------------------------------
sub send_finish_email_to_user
{
	my $email_subject;
	my $HttpPath=$VARS{run_url}."/".$VARS{output_page};
	if ($FORM{JOB_TITLE} ne "")
	{
		$email_subject = "'Your Guidance results for $FORM{JOB_TITLE} are ready'";
	}
	elsif ($FORM{usrSeq_File} ne "")
	{
		$email_subject = "'Your Guidance results for $FORM{usrSeq_File} are ready'";
	}
	else
	{
		$email_subject = "'Your Guidance results for run number $VARS{run_number} are ready'";
	}
	my $email_message = "'Hello,\\n\\nThe results for your Guidance run are ready at:\\n".$HttpPath."\\n\\nRunning Parameters:\\n";
	$email_message.="Job Title:$FORM{JOB_TITLE}\\n" if ($FORM{JOB_TITLE} ne "");
	$email_message.="Sequences File: $FORM{usrSeq_File}\\nMSA Algorithm: $FORM{MSA_Program}\\nNumber of Bootstraps: $FORM{Bootstraps}\\n";
	$email_message.="Scoreing Method: $FORM{PROGRAM}\\n";
	
	$email_message.="\\nPlease note: the results will be kept on the server for three months.\\n\\nThanks\\nGUIDANCE Team'";

        my $msg = './sendEmail.pl -f \''.GENERAL_CONSTANTS::ADMIN_EMAIL.'\' -t \''.$FORM{user_mail}.'\' -u '.$email_subject.' -xu '.$VARS{userName}.' -xp '.$VARS{userPass}.' -s '.$VARS{smtp_server}.' -m '.$email_message;
#	my $msg = "ssh bioseq\@lecs \"cd $VARS{send_email_dir}; ".'perl sendEmail.pl -f \''.GENERAL_CONSTANTS::ADMIN_EMAIL.'\' -t \''.$FORM{user_mail}.'\' -u '.$email_subject.' -xu '.$VARS{userName}.' -xp '.$VARS{userPass}.' -s '.$VARS{smtp_server}.' -m '.$email_message."\""; # TO ACTIVATE IF THE NODES IN CLUSTER FAILS TO COMMUNICATE WITH NET
	print LOG "MESSAGE:$email_message\nCOMMAND:$msg\n";
	chdir $VARS{send_email_dir};
	my $email_system_return = `$msg`;
	unless ($email_system_return =~ /successfully/)    {
		print LOG "send_mail: The message was not sent successfully. system returned: $email_system_return\n";
	}
	
}
#---------------------------------------------
sub update_output_that_run_finished
#---------------------------------------------
{
	my $OutHtmlFile=shift;
	### write to the output that the job has finished and stop reload
	open OUTPUT, "<$OutHtmlFile";
	my @output = <OUTPUT>;
	close OUTPUT;

	open OUTPUT, ">$OutHtmlFile";
	foreach my $line (@output){
		if (($line=~/REFRESH/) or ($line=~/NO-CACHE/))
		{
			next;
		}
		elsif ($line=~/(.*)RUNNING(.*)/)
		{
			print OUTPUT $1."FINISHED".$2;
		}
		else {
			print OUTPUT $line;
		}
	}
	print OUTPUT "<hr>\n";
	print OUTPUT "<h4 class=footer><p align='center'>Questions and comments are welcome! Please <span class=\"admin_link\"><a href=\"mailto:bioSequence\@tauex.tau.ac.il\?subject=Guidance\%20Run\%20Number\%20$VARS{run_number}\">contact us</a></span></p></h4>";
	print OUTPUT "<div id=\"bottom_links\"> <!-- links before the footer -->\n";
	
	print OUTPUT "<span class=\"bottom_link\">\n";
	print OUTPUT "<a href=\"$home_URL\" target=\"_blank\">Home</a>\n";
	print OUTPUT "&nbsp;|&nbsp\n";
	print OUTPUT "<a href=\"$overview_URL\" target=\"_blank\">Overview</a>\n";
	print OUTPUT "&nbsp;|&nbsp;\n";
	print OUTPUT "<a href=\"$gallery_URL\" target=\"_blank\">Gallery</a>\n";
	print OUTPUT "&nbsp;|&nbsp;\n";
	print OUTPUT "<a href=\"http://guidance.tau.ac.il/credits.html\" target=\"_blank\">Credits</a>\n";
	print OUTPUT "</span>\n";
	print OUTPUT "<br />\n";
	print OUTPUT "</div>\n";
	
	print OUTPUT "</body>\n";
	print OUTPUT "</html>\n";
	close OUTPUT;
}
#---------------------------------------------
sub print_message_to_output{
#---------------------------------------------
    my $msg = shift;
    print OUTPUT "\n<ul><li>$msg</li></ul>\n";
}


#---------------------------------------------
sub Bootstrap_Trees
#---------------------------------------------
{
	if ($isServer == 1) {
		print_message_to_output("Constructing bootstrap guide-trees");
	}
	mkdir($VARS{BootStrap_Dir});
	$VARS{Tree_File}="$VARS{dataset}.".$FORM{MSA_Program}.".semphy.outTree";
	$VARS{Semphy_OutFile}="$VARS{dataset}.".$FORM{MSA_Program}.".semphy.out";
	$VARS{Semphy_LogFile}="$VARS{dataset}.".$FORM{MSA_Program}.".semphy.log";
	$VARS{Semphy_StdFile}="$VARS{dataset}.".$FORM{MSA_Program}.".semphy.std";
	my $cmd="";
	my $MSA_depth=MSA_Depth("$VARS{WorkingDir}$VARS{Alignment_File}");
	my $verbose_level=8;
	if (uc $VARS{BBL} eq "YES")
	{
		$semphy_prog=$semphy_prog." -n ";
		$verbose_level=1;
	}
	if (($FORM{Seq_Type} eq "AminoAcids") or ($FORM{Seq_Type} eq "Codons")){
		if ($MSA_depth>150) # use JC for distance estimation 
		{
			$cmd=$semphy_prog." -a 20 --aaJC -H -J -v $verbose_level --BPrepeats=".$FORM{Bootstraps}." -s ".$VARS{WorkingDir}.$VARS{Alignment_File}." -o ".$VARS{BootStrap_Dir}.$VARS{Semphy_OutFile}." -T ".$VARS{BootStrap_Dir}.$VARS{Tree_File}." -l ".$VARS{BootStrap_Dir}.$VARS{Semphy_LogFile}." > ".$VARS{BootStrap_Dir}.$VARS{Semphy_StdFile};#." 2>\&1";
		}
		else # use JTT for distance estimation
		{
			$cmd=$semphy_prog." -a 20 --jtt -H -J -v $verbose_level --BPrepeats=".$FORM{Bootstraps}." -s ".$VARS{WorkingDir}.$VARS{Alignment_File}." -o ".$VARS{BootStrap_Dir}.$VARS{Semphy_OutFile}." -T ".$VARS{BootStrap_Dir}.$VARS{Tree_File}." -l ".$VARS{BootStrap_Dir}.$VARS{Semphy_LogFile}." > ".$VARS{BootStrap_Dir}.$VARS{Semphy_StdFile};#." 2>\&1";
		}
	}
	elsif ($FORM{Seq_Type} eq "Nucleotides") {
		if ($MSA_depth>150) # use JC for distance estimation 
		{
			$cmd=$semphy_prog." -a 4 --nucjc -J -H -v $verbose_level --BPrepeats=".$FORM{Bootstraps}." -s ".$VARS{WorkingDir}.$VARS{Alignment_File}." -o ".$VARS{BootStrap_Dir}.$VARS{Semphy_OutFile}." -T ".$VARS{BootStrap_Dir}.$VARS{Tree_File}." -l ".$VARS{BootStrap_Dir}.$VARS{Semphy_LogFile}." > ".$VARS{BootStrap_Dir}.$VARS{Semphy_StdFile};#." 2>\&1";
		}
		else # use HKY for distance estimation
		{
			$cmd=$semphy_prog." -a 4 --hky -J -H -v $verbose_level --BPrepeats=".$FORM{Bootstraps}." -s ".$VARS{WorkingDir}.$VARS{Alignment_File}." -o ".$VARS{BootStrap_Dir}.$VARS{Semphy_OutFile}." -T ".$VARS{BootStrap_Dir}.$VARS{Tree_File}." -l ".$VARS{BootStrap_Dir}.$VARS{Semphy_LogFile}." > ".$VARS{BootStrap_Dir}.$VARS{Semphy_StdFile};#." 2>\&1";
		}
	}
	print LOG "Bootstrap_Trees: $cmd\n";
#	if (!-e "$VARS{BootStrap_Dir}$VARS{Semphy_LogFile}") # JUST TO SAVE TIME...
#	{
		system ($cmd);
#	}
	if (-z "$VARS{BootStrap_Dir}$VARS{Semphy_LogFile}" or !-e "$VARS{BootStrap_Dir}$VARS{Semphy_LogFile}")
	{
		exit_on_error("sys_error","Bootstrap_Trees: '$VARS{BootStrap_Dir}$VARS{Semphy_LogFile}' is empty/does not exist\n");
	}
}

#---------------------------------------------
sub Align{
#---------------------------------------------
	if ($isServer == 1) {
		print_message_to_output("Generating the base alignment");
	}
	$VARS{Alignment_File}="$VARS{dataset}.".$FORM{MSA_Program}.".aln";           
	$VARS{Core_Alignment_File}="$VARS{dataset}.".$FORM{MSA_Program}.".CORE.aln";
	if ($FORM{MSA_Program} eq "MAFFT")
	# align with mafft
	{
		# create the core alignment
		my $cmd="";
		if ($VARS{align_param}=~/addfragments/)
		{
			# create the core alignment
			my $aln_param_for_core_msa=$VARS{align_param};
			my @tmp=split (/\Q\-\-\E/,$aln_param_for_core_msa);
			my $tmp_size=@tmp;
			if ($tmp_size>=1) # command line type
			{
				for (my $i=0;$i<$tmp_size;$i++)
				{
					if ($tmp[$i]=~/addfragments/){delete $tmp[$i];}
					elsif ($tmp[$i]=~/multipair/){delete $tmp[$i];}
					elsif ($tmp[$i]=~/6merpair/){delete $tmp[$i];}
				}
				$aln_param_for_core_msa=join("\\\-\\\-",@tmp);
				#print "1:$aln_param_for_core_msa\n";<STDIN>;
			}
			else
			{
				#for servers
				my @tmp=split (/\Q--\E/,$aln_param_for_core_msa);
				my $tmp_size=@tmp;
				for (my $i=0;$i<$tmp_size;$i++)
				{
					if ($tmp[$i]=~/addfragments/){delete $tmp[$i];}
					elsif ($tmp[$i]=~/multipair/){delete $tmp[$i];}
					elsif ($tmp[$i]=~/6merpair/){delete $tmp[$i];}
				}
				$aln_param_for_core_msa=join ("\\\-\\\-",@tmp);
				#print "2:$aln_param_for_core_msa\n";<STDIN>;
			}
			if (($FORM{Seq_Type} eq "AminoAcids") or ($FORM{Seq_Type} eq "Codons")){
				$cmd=$mafft_prog." $aln_param_for_core_msa --amino  --quiet --thread ".$VARS{proc_num}." ".$VARS{WorkingDir}.$VARS{codded_seq_fileName}." > ".$VARS{WorkingDir}.$VARS{Core_Alignment_File};#." 2> ".$VARS{WorkingDir}.$VARS{Alignment_File}.".std";
			}
			elsif ($FORM{Seq_Type} eq "Nucleotides") {
				$cmd=$mafft_prog." $aln_param_for_core_msa --nuc --quiet  --thread ".$VARS{proc_num}." ".$VARS{WorkingDir}.$VARS{codded_seq_fileName}." > ".$VARS{WorkingDir}.$VARS{Core_Alignment_File};#." 2> ".$VARS{WorkingDir}.$VARS{Alignment_File}.".std";
			}
			print LOG "Core Align: $cmd\n";
			system ($cmd);
			if ((!-e "$VARS{WorkingDir}$VARS{Core_Alignment_File}") or (-z "$VARS{WorkingDir}$VARS{Core_Alignment_File}"))
			{
				exit_on_error("sys_error","Align: '$VARS{WorkingDir}$VARS{Core_Alignment_File}' is empty/does not exist\n");
			}
			## create the full alignment with fragments
			if (($FORM{Seq_Type} eq "AminoAcids") or ($FORM{Seq_Type} eq "Codons")){
				$cmd=$mafft_prog." $VARS{align_param} --amino  --quiet ".$VARS{WorkingDir}.$VARS{Core_Alignment_File}." > ".$VARS{WorkingDir}.$VARS{Alignment_File};#." 2> ".$VARS{WorkingDir}.$VARS{Alignment_File}.".std";
			}
			elsif ($FORM{Seq_Type} eq "Nucleotides") {
				$cmd=$mafft_prog." $VARS{align_param} --nuc --quiet ".$VARS{WorkingDir}.$VARS{Core_Alignment_File}." > ".$VARS{WorkingDir}.$VARS{Alignment_File};#." 2> ".$VARS{WorkingDir}.$VARS{Alignment_File}.".std";
			}
		}	
		else
		{
			if (($FORM{Align_Order} eq "aligned") and ($VARS{align_param}!~/reorder/)) {$VARS{align_param}.=" --reorder";}
			if (($FORM{Seq_Type} eq "AminoAcids") or ($FORM{Seq_Type} eq "Codons")){
				$cmd=$mafft_prog." $VARS{align_param} --amino  --quiet ".$VARS{WorkingDir}.$VARS{codded_seq_fileName}." > ".$VARS{WorkingDir}.$VARS{Alignment_File};#." 2> ".$VARS{WorkingDir}.$VARS{Alignment_File}.".std";
			}
			elsif ($FORM{Seq_Type} eq "Nucleotides") {
				$cmd=$mafft_prog." $VARS{align_param} --nuc --quiet ".$VARS{WorkingDir}.$VARS{codded_seq_fileName}." > ".$VARS{WorkingDir}.$VARS{Alignment_File};#." 2> ".$VARS{WorkingDir}.$VARS{Alignment_File}.".std";
			}
		}
		print LOG "Align: $cmd\n";
		system ($cmd);
	}
	elsif ($FORM{MSA_Program} eq "MAFFT_LINSI")
	# align with mafft
	{
		my $cmd="";
		if (($FORM{Seq_Type} eq "AminoAcids") or ($FORM{Seq_Type} eq "Codons")){
			$cmd=$mafft_prog." --localpair --maxiterate 1000 --amino --quiet ".$VARS{WorkingDir}.$VARS{codded_seq_fileName}." > ".$VARS{WorkingDir}.$VARS{Alignment_File};#." 2> ".$VARS{WorkingDir}.$VARS{Alignment_File}.".std";
		}
		elsif ($FORM{Seq_Type} eq "Nucleotides") {
			$cmd=$mafft_prog." --localpair --maxiterate 1000 --nuc --quiet ".$VARS{WorkingDir}.$VARS{codded_seq_fileName}." > ".$VARS{WorkingDir}.$VARS{Alignment_File};#." 2> ".$VARS{WorkingDir}.$VARS{Alignment_File}.".std";
		}
		print LOG "Align: $cmd\n";
		system ($cmd);
	}
	elsif ($FORM{MSA_Program} eq "PAGAN")
	{
		my $cmd="";
		my $RoughTree_MAFFT="RoughTree_MAFFT_globalpair.tree";
        # builed estimated tree with mafft
		system ("$mafft_prog --retree 0 --treeout --globalpair --reorder $VARS{WorkingDir}$VARS{codded_seq_fileName} >/dev/null");
		# fix the tree (add semi-colon and remove the ___)
		system ("mv $VARS{WorkingDir}$VARS{codded_seq_fileName}.tree $VARS{WorkingDir}$RoughTree_MAFFT");
		fix_MAFFT_RoughTree("$VARS{WorkingDir}$RoughTree_MAFFT");
		# pagan cmd
		$cmd=$pagan_prog." --seqfile $VARS{WorkingDir}$VARS{codded_seq_fileName} --treefile $VARS{WorkingDir}$RoughTree_MAFFT --outfile $VARS{WorkingDir}$VARS{Alignment_File}";
		print LOG "Align: $cmd\n";
		system ($cmd);
		move("$VARS{WorkingDir}$VARS{Alignment_File}.fas", "$VARS{WorkingDir}$VARS{Alignment_File}");
	}
	elsif ($FORM{MSA_Program} eq "PRANK")
	{
		my $PRANK_VERSION=0;
		# if PRANK find out it version
		my $prank_help=`$prank_prog`;
		if ($prank_help=~/.*prunedata.*/){$PRANK_VERSION="121218";}
		elsif ($prank_help=~/.*showanc.*/){$PRANK_VERSION="120626";}
		# print "PRANK VERSION:>=$PRANK_VERSION\n";
		my $cmd="";
		if (($FORM{Seq_Type} eq "AminoAcids") or  ($FORM{Seq_Type} eq "Nucleotides") or  ($FORM{Seq_Type} eq "Codons")){
			if ($PRANK_VERSION==0)
			{
				$cmd=$prank_prog." $VARS{align_param} -quiet -d=".$VARS{WorkingDir}.$VARS{codded_seq_fileName}." -o=".$VARS{WorkingDir}.$VARS{Alignment_File}." -noxml > ".$VARS{WorkingDir}.$VARS{Alignment_File}.".std";# 2>\&1"; # -noxml is not supported anymore
			}
			else # version >=120626
			{
				$cmd=$prank_prog." $VARS{align_param} -quiet -d=".$VARS{WorkingDir}.$VARS{codded_seq_fileName}." -o=".$VARS{WorkingDir}.$VARS{Alignment_File}." > ".$VARS{WorkingDir}.$VARS{Alignment_File}.".std";# 2>\&1"; # -noxml is not supported anymore
			}
		}
		print LOG "Align: $cmd\n";
		if (!-e "$VARS{WorkingDir}$VARS{Alignment_File}.2.fas") # Just to save time
		{
			system ($cmd);
		}
		if ($PRANK_VERSION eq "121218") # The output file now named best.fas
		{
			system ("cp $VARS{WorkingDir}$VARS{Alignment_File}.best.fas $VARS{WorkingDir}$VARS{Alignment_File}");
		}
		else
		{
			system ("cp $VARS{WorkingDir}$VARS{Alignment_File}.2.fas $VARS{WorkingDir}$VARS{Alignment_File}");
		}

		if ($FORM{Align_Order} eq "as_input")
		{
			MSA_parser::sort_alignment("$VARS{WorkingDir}$VARS{Alignment_File}","fasta");
			$VARS{Alignment_File_NOT_SORTED}=$VARS{Alignment_File};
			$VARS{Alignment_File}=$VARS{Alignment_File}.".Sorted";
		}
	}
	elsif ($FORM{MSA_Program} eq "CLUSTALW")
	{
		my $cmd="";
		if (($FORM{Seq_Type} eq "AminoAcids")  or ($FORM{Seq_Type} eq "Codons")){
			$cmd=$clustalw_prog." -quiet -infile=".$VARS{WorkingDir}.$VARS{codded_seq_fileName}." -outfile=".$VARS{WorkingDir}.$VARS{Alignment_File}." -TYPE=PROTEIN > ".$VARS{WorkingDir}.$VARS{Alignment_File}.".std";# 2>\&1";	
	}
		elsif ($FORM{Seq_Type} eq "Nucleotides") {
			$cmd=$clustalw_prog." -quiet -infile=".$VARS{WorkingDir}.$VARS{codded_seq_fileName}." -outfile=".$VARS{WorkingDir}.$VARS{Alignment_File}." -TYPE=DNA > ".$VARS{WorkingDir}.$VARS{Alignment_File}.".std";# 2>\&1";
		}
		print LOG "Align: $cmd";
		system ($cmd);
		MSA_parser::convert_msa_format("$VARS{WorkingDir}$VARS{Alignment_File}","clustalw","$VARS{WorkingDir}$VARS{Alignment_File}.fs","fasta");
		system ("mv $VARS{WorkingDir}$VARS{Alignment_File} $VARS{WorkingDir}$VARS{Alignment_File}.orig");
		system ("mv $VARS{WorkingDir}$VARS{Alignment_File}.fs $VARS{WorkingDir}$VARS{Alignment_File}");
		
		if ($FORM{Align_Order} eq "as_input")
		{
			print "MSA_parser::sort_alignment($VARS{WorkingDir}$VARS{Alignment_File},fasta);\n";
			my @ans=MSA_parser::sort_alignment("$VARS{WorkingDir}$VARS{Alignment_File}","fasta");
			print join("",@ans);
			$VARS{Alignment_File_NOT_SORTED}=$VARS{Alignment_File};
			$VARS{Alignment_File}=$VARS{Alignment_File}.".Sorted";
		}
	}
	elsif ($FORM{MSA_Program} eq "MUSCLE")
	{
		my $cmd="";
		if (($FORM{Seq_Type} eq "AminoAcids")  or ($FORM{Seq_Type} eq "Codons")){
			$cmd=$muscle_prog." -quiet -in ".$VARS{WorkingDir}.$VARS{codded_seq_fileName}." -out ".$VARS{WorkingDir}.$VARS{Alignment_File}." -seqtype protein > ".$VARS{WorkingDir}.$VARS{Alignment_File}.".std";# 2>\&1";
		}
		elsif ($FORM{Seq_Type} eq "Nucleotides") {
			$cmd=$muscle_prog." -quiet -in ".$VARS{WorkingDir}.$VARS{codded_seq_fileName}." -out ".$VARS{WorkingDir}.$VARS{Alignment_File}." -seqtype dna > ".$VARS{WorkingDir}.$VARS{Alignment_File}.".std";# 2>\&1";
		}
		print LOG "MUSCLE Align: $cmd";
		system ($cmd);
		if ($FORM{Align_Order} eq "as_input")
		{
			print "MSA_parser::sort_alignment($VARS{WorkingDir}$VARS{Alignment_File},fasta);\n";
			my @ans=MSA_parser::sort_alignment("$VARS{WorkingDir}$VARS{Alignment_File}","fasta");
			print join("",@ans);
			$VARS{Alignment_File_NOT_SORTED}=$VARS{Alignment_File};
			$VARS{Alignment_File}=$VARS{Alignment_File}.".Sorted";
		}
	}
	if ((!-e "$VARS{WorkingDir}$VARS{Alignment_File}") or (-z "$VARS{WorkingDir}$VARS{Alignment_File}"))
	{
		exit_on_error("sys_error","Align: '$VARS{WorkingDir}$VARS{Alignment_File}' is empty/does not exist\n");
	}

}

#---------------------------------------------
sub exit_on_error{
#---------------------------------------------
    my $which_error = shift;
    my $error_msg = shift;
    my $error_definition = "<font size=+1 color='red'>ERROR! GUIDANCE session has been terminated:</font><br />\n";
    my $syserror = "<font size=+1 color='red'>A SYSTEM ERROR OCCOURED!</font><br />Plesae try to run GUIDANCE again in a few minutes.<br />We apologize for the inconvenience.<br />\n";
	if ($isServer == 0) {
		$syserror = "Guidance error\n";
		$error_definition = "Guidance error: ";
	}

   	if ($isServer == 1) {
		open OUTPUT, ">>$VARS{WorkingDir}/$VARS{output_page}";
		if ($which_error eq 'user_error'){
			print LOG "\nEXIT on error:\n$error_msg\n";
			print OUTPUT  $error_definition."$error_msg";
	        # print $error_msg to the screen
		}
		elsif($which_error eq 'sys_error'){
			send_administrator_mail_on_error($error_msg);
			print LOG "\n$error_msg\n";
			print OUTPUT $syserror;
		}
		close OUTPUT;
		# finish the output page
		sleep 10;
		open OUTPUT, "$VARS{WorkingDir}/$VARS{output_page}";
		my @output = <OUTPUT>;
		close OUTPUT;
		# remove the refresh commands from the output page
		open OUTPUT, ">$VARS{WorkingDir}/$VARS{output_page}";
		foreach my $line (@output){
			if (($line=~/TTP-EQUIV="REFRESH"/) or ($line=~/CONTENT="NO-CACHE"/))
			{
				next;
			}
			elsif ($line=~/(.*)RUNNING(.*)/)
			{
				print OUTPUT $1."FAILED".$2;
			}
			else {
				print OUTPUT $line;
			}
        }
		print OUTPUT "<hr> <h4 class=footer><p align='center'>\nQuestions and comments are welcome! Please <span class=\"admin_link\"><a href=\"mailto:bioSequence\@tauex.tau.ac.il\?subject=GUIDANCE\%20Run\%20Number\%20$VARS{run_number}\">contact us</a></span></p></h4>\n<div id=\"bottom_links\"> <!-- links before the footer --><span class=\"bottom_link\"> <a href=\"$home_URL\" target=\"_blank\">Home</a> &nbsp;|&nbsp<a href=\"$overview_URL\" target=\"_blank\">Overview</a> &nbsp;|&nbsp;<a href=\"$gallery_URL\" target=\"_blank\">Gallery</a> &nbsp;|&nbsp;<a href=\"/credits.html\" target=\"_blank\">Credits</a> </span> <br /> </div>";
		print OUTPUT "</body>\n";
		print OUTPUT "</html>\n";	
		close OUTPUT;
		if ($FORM{user_email} ne "")
		{
			send_mail_on_error();
		}

		print LOG "\nExit Time: ".(BIOSEQUENCE_FUNCTIONS::printTime)."\n";
		close LOG;
		chmod (0755, $VARS{WorkingDir});
	} else {
		if ($which_error eq 'user_error'){
			print LOG "\nEXIT on error:\n$error_msg\n";
			print "ERROR: $error_msg\n";
	        # print $error_msg to the screen
		}
		elsif($which_error eq 'sys_error'){
			print LOG "\n$error_msg\n";
			print "ERROR: $error_msg\n";
			print $syserror."\n";
		}
	}
	
	if (($FORM{PROGRAM} eq "GUIDANCE") and ($isServer == 1)) # Zip BP dir on server
	{
        # tar and remove the BP dir
		$cmd="cd $VARS{WorkingDir};tar -czf $VARS{Output_Prefix}_BP_Dir.tar.gz ./BP";
        #print "$cmd\n";
		system ($cmd);
		system ("rm -r -f $VARS{BootStrap_Dir}") if (-e "$VARS{Output_Prefix}_BP_Dir.tar.gz");
	}
	exit;
}

#---------------------------------------------
sub send_mail_on_error
#---------------------------------------------
{
        my $email_subject;
        my $HttpPath=$VARS{run_url}.$VARS{output_page};
        $email_subject = "'Your GUIDANCE run for $FORM{usrSeq_File} FAILED'";
        my $email_message = "'Hello,\\n\\nUnfortunately your GUIDANCE run (number ".$VARS{run_number}.") has failed.\\nPlease have a look at ".$HttpPath." for further details\\n\\nSorry for the inconvenience\\nGUIDANCE Team'";

#		my $msg = "ssh bioseq\@lecs \"cd $VARS{send_email_dir};".'perl sendEmail.pl -f \''.GENERAL_CONSTANTS::ADMIN_EMAIL.'\' -t \''.$FORM{user_email}.'\' -u '.$email_subject.' -xu '.$VARS{userName}.' -xp '.$VARS{userPass}.' -s '.$VARS{smtp_server}.' -m '.$email_message."\""; # ACTIVATE IN CASE THE CLUSTER NODE FAILS TO COMMUNICATE WITH THE NET
        my $msg = './sendEmail.pl -f \''.GENERAL_CONSTANTS::ADMIN_EMAIL.'\' -t \''.$FORM{user_email}.'\' -u '.$email_subject.' -xu '.$VARS{userName}.' -xp '.$VARS{userPass}.' -s '.$VARS{smtp_server}.' -m '.$email_message;
		
        #if ($attach ne ''){$msg.=" -a $attach"; print LOG "sending $msg\n";}
       print LOG "MESSAGE:$email_message\nCOMMAND:$msg\n";
       chdir $VARS{send_email_dir};
       my $email_system_return = `$msg`;
       unless ($email_system_return =~ /successfully/)    {
            print LOG "send_mail: The message was not sent successfully. system returned: $email_system_return\n";
            }
}

#---------------------------------------------
sub send_administrator_mail_on_error
#---------------------------------------------
{
        my $message=shift;
        my $email_subject;
        $email_subject = "'SYSTEM ERROR has occurred on GUIDANCE: $VARS{run_url}'";
        my $email_message = "'Hello,\\n\\nUnfortunately a system SYSTEM ERROR has occurred on GUIDANCE: $VARS{run_url}.\\nERROR: $message.'";
        my $Admin_Email=GENERAL_CONSTANTS::ADMIN_EMAIL;
#		my $msg = "ssh bioseq\@lecs \" cd $VARS{send_email_dir};".'perl sendEmail.pl -f \'bioSequence@tauex.tau.ac.il\' -t \''."bioSequence\@tauex.tau.ac.il".'\' -u '.$email_subject.' -xu '.$VARS{userName}.' -xp '.$VARS{userPass}.' -s '.$VARS{smtp_server}.' -m '.$email_message."\""; # ACTIVATE IN CASE THE CLUSTER NODE FAILS TO COMMUNICATE WITH THE NET
		my $msg = $VARS{send_email_dir}.'/sendEmail.pl -f \'bioSequence@tauex.tau.ac.il\' -t \''."bioSequence\@tauex.tau.ac.il".'\' -u '.$email_subject.' -xu '.$VARS{userName}.' -xp '.$VARS{userPass}.' -s '.$VARS{smtp_server}.' -m '.$email_message;
		#  print LOG "MESSAGE:$email_message\nCOMMAND:$msg\n";
        chdir $VARS{send_email_dir};
        my $email_system_return = `$msg`;
	}

sub Convert_to_CSV
{
    my $In=shift;
    my $Out=shift;
    open (IN,$In) or return ("Can't Open In file: '$In' $!");
    open (OUT,">$Out") or return ("Can't Open Output file: '$Out' $!");
    while (my $line=<IN>)
	{
		chomp($line);
		$line=trim($line);
		my @line=split(/\s+/,$line);
		my $new_line="";
		foreach my $element (@line)
		{
			$new_line.=$element.",";
		}
		chop ($new_line);
		$line=~s/\#//;
		if ($line !~/END/)
		{
			print OUT $new_line,"\n";
		}
	}
    close (IN);
    close (OUT);
    return "OK";
}

sub Convert_to_CSV_Col_Scores
  {
    my $In=shift;
    my $Out=shift;
    my $last_pos=0;
    open (IN,$In) or return ("Can't Open In file: '$In' $!");
    open (OUT,">$Out") or return ("Can't Open Output file: '$Out' $!");
    while (my $line=<IN>)
	{
		chomp($line);
		$line=trim($line);
		my @line=split(/\s+/,$line);
		my $new_line="";
		if ($line[0]=~/[0-9]+/)
		{
			my $current_pos=$line[0];
			if ($current_pos-1==$last_pos)
			{
				print OUT "$current_pos,$line[1]";
				$last_pos=$current_pos;
			}
			elsif ($current_pos-1>$last_pos)
			{
				while ($current_pos-1>$last_pos)
				{
					print OUT "$current_pos,1\n";
					$last_pos++;
				}
			}
		}
		$line=~s/\#//;
		if ($line !~/END/)
		{
			print OUT $new_line,"\n";
		}
	}
    close (IN);
    close (OUT);
    return "OK";
}
sub trim($)
{
	my $string = shift;
	$string =~ s/^[\s\t]+//;
	$string =~ s/[\s\t]+$//;
	return $string;
}
sub print_remove_site_selection_mask
{
    my $sp_res_file=shift;
    my $msa_length=shift;
    
    my $return_str="";
    my %Cutoff_Removed_Pos_Hash; # hash that will hold for each cutoff the precentage of positions remained
    open (SP_RES_FILE,$sp_res_file); #|| exit_on_error(...)
	my $c1=0;
    while (my $line=<SP_RES_FILE>)
	{
		$c1++;
		$line=trim($line);
		if ($line=~/^[0-9]+/)
		{
			my @line=split(/\s+/,$line);
			my $score=sprintf("%.3f",$line[2]);
			
			if (defined ($Cutoff_Removed_Pos_Hash{$score})){$Cutoff_Removed_Pos_Hash{$score}=$Cutoff_Removed_Pos_Hash{$score}+1;}
			else {$Cutoff_Removed_Pos_Hash{$score}=1;}
		}
	}
    
    $return_str.="<select name=\"cutoff\" id=\"cutoff\">";
    foreach my $cutoff (sort {$a<=>$b} keys %Cutoff_Removed_Pos_Hash){
		my $res_bellow_Cutoof=0;
		foreach my $score (keys %Cutoff_Removed_Pos_Hash)
		{
			if ($score<$cutoff)
			{
				$res_bellow_Cutoof+=$Cutoff_Removed_Pos_Hash{$score};
			}
		}
		my $res_remain_precent=sprintf("%.3f",(1-($res_bellow_Cutoof/$c1)));
		$res_remain_precent=$res_remain_precent*100;
		$return_str.="<option value=\"$cutoff\">$cutoff ($res_remain_precent\% of residues remain)</option>" if (($res_remain_precent !=100) and ($res_remain_precent !=0));
    }
	
    close (SP_RES_FILE);
    return $return_str;
}

sub print_remove_site_selection_box
{
    my $sp_sites_file=shift;
    my $msa_length=shift;
    
    my $return_str="";
    my %Cutoff_Removed_Pos_Hash; # hash that will hold for each cutoff the precentage of positions remained
    open (SP_COL_FILE,$sp_sites_file); #|| exit_on_error(...)
    while (my $line=<SP_COL_FILE>)
	{
		$line=trim($line);
		if ($line=~/^[0-9]+/)
		{
			my @line=split(/\s+/,$line);
			my $score=sprintf("%.3f",$line[1]);
			
			if (defined ($Cutoff_Removed_Pos_Hash{$score})){$Cutoff_Removed_Pos_Hash{$score}=$Cutoff_Removed_Pos_Hash{$score}+1;}
			else {$Cutoff_Removed_Pos_Hash{$score}=1;}
		}
	}
    $return_str.="<input type=\"hidden\" name=\"run_num\" value=\"$VARS{run_number}\">\n";
    $return_str.="<input type=\"hidden\" name=\"VARS\" value=\"$stored_data_file\"/>\n";
    $return_str.="<select name=\"Col_Cutoff\" id=\"Col_Cutoff\">\n";
    foreach my $cutoff (sort {$a<=>$b} keys %Cutoff_Removed_Pos_Hash){
		my $Cols_bellow_Cutoof=0;
		foreach my $score (keys %Cutoff_Removed_Pos_Hash)
		{
			if ($score<$cutoff)
			{
				$Cols_bellow_Cutoof+=$Cutoff_Removed_Pos_Hash{$score};
			}
		}
		my $Col_remain_precent=sprintf("%.3f",(1-($Cols_bellow_Cutoof/$msa_length)));
		$Col_remain_precent=$Col_remain_precent*100;
		$return_str.="<option value=\"$cutoff\">$cutoff ($Col_remain_precent\% of columns remain)</option>" if (($Col_remain_precent !=100) and ($Col_remain_precent !=0));
    }
    $return_str.="<input type=\"submit\" value=\"Remove Columns\"/>\n";
    $return_str.="</select><font size =\"-1\">(see <a href=\"$overview_URL#Remove_COL\" target=\"_blank\">help</a></font>)\n";
    $return_str.="</form>\n";
    close (SP_COL_FILE);
    return $return_str;
}

sub print_remove_seq_selection_box
{
    my $sp_seq_file=shift;
    my $msa_depth=shift;
    
    my $return_str="";
    my %Cutoff_Removed_Seq_Hash; # hash that will hold for each cutoff the precentage of positions remained
    open (SP_SEQ_FILE,$sp_seq_file); #|| exit_on_error(...)
    while (my $line=<SP_SEQ_FILE>)
	{
		$line=trim($line);
		my @line=split(/\s+/,$line);
		if ($line[1]=~/^[0-9]+/)
		{
			my $score=sprintf("%.3f",$line[1]);
			if (defined ($Cutoff_Removed_Seq_Hash{$score})){$Cutoff_Removed_Seq_Hash{$score}=$Cutoff_Removed_Seq_Hash{$score}+1;}
			else {$Cutoff_Removed_Seq_Hash{$score}=1;}
		}
	}
    $return_str.="<input type=\"hidden\" name=\"run_num\" value=\"$VARS{run_number}\">\n";
    $return_str.="<input type=\"hidden\" name=\"VARS\" value=\"$stored_data_file\"/>\n";
	$return_str.="<input type=\"hidden\" name=\"FORM\" value=\"$stored_form_data\"/>\n";
    $return_str.="<select name=\"Seq_Cutoff\" id=\"Seq_Cutoff\">\n";
    foreach my $cutoff (sort {$a<=>$b} keys %Cutoff_Removed_Seq_Hash){
		my $Seq_bellow_Cutoof=0;
		foreach my $score (keys %Cutoff_Removed_Seq_Hash)
		{
#	  print "CUT:$cutoff\tSCORE:$score\t$Cutoff_Removed_Seq_Hash{$score}\n";
	  if ($score<$cutoff)
	    {
	      $Seq_bellow_Cutoof+=$Cutoff_Removed_Seq_Hash{$score};
	    }
	}
#      print "\n===========\n";
		my $Seq_remain_precent=sprintf("%.3f",(1-($Seq_bellow_Cutoof/$msa_depth)));
		$Seq_remain_precent=$Seq_remain_precent*100;
		$return_str.="<option value=\"$cutoff\">$cutoff ($Seq_remain_precent\% of sequences remain)</option>" if (($Seq_remain_precent !=100) and ($Seq_remain_precent !=0));
    }
    $return_str.="<input type=\"submit\" value=\"Remove Seqs\"/>\n";
	
    $return_str.="</select><font size =\"-1\">(see <a href=\"$overview_URL#Remove_Seq\" target=\"_blank\">help</a></font>)\n";
    $return_str.="</form>\n";
    close (SP_SEQ_FILE);
    return $return_str;
}

sub Prepare_rerun_param
{
	my %MAFFT_max_iterates=(
							0 => '0',
							1 => '1',
							2 => '2',
							5 => '3',
							10 => '4',
							20 => '5',
							50 => '6',
							80 => '7',
							100 => '8',
							1000 => '9',
							);
	
	open (PARAM,">$VARS{WorkingDir}/rerun_param");
	print PARAM	"<?php\n";
	if ($FORM{PROGRAM} eq "GUIDANCE"){
		print PARAM "\$PROGRAM=0;\n";
		if ($FORM{MSA_Program} eq "MAFFT") {print PARAM "\$MSA_Prog=0;\n";}
		elsif ($FORM{MSA_Program} eq "PRANK") {print PARAM "\$MSA_Prog=1;\n";}
		elsif ($FORM{MSA_Program} eq "CLUSTALW") {print PARAM "\$MSA_Prog=2;\n";}
		elsif ($FORM{MSA_Program} eq "MUSCLE") {print PARAM "\$MSA_Prog=3;\n";}
		elsif ($FORM{MSA_Program} eq "PAGAN") {print PARAM "\$MSA_Prog=4;\n";}
	}
	elsif ($FORM{PROGRAM} eq "HoT") {
		print PARAM "\$PROGRAM=1;\n";
		if ($FORM{MSA_Program} eq "MAFFT") {print PARAM "\$MSA_Prog=0;\n";}
		elsif ($FORM{MSA_Program} eq "PRANK") {print PARAM "\$MSA_Prog=1;\n";}
		elsif ($FORM{MSA_Program} eq "CLUSTALW") {print PARAM "\$MSA_Prog=2;\n";}
	}
	if ($FORM{MSA_Program} eq "MAFFT") {
		if ($FORM{MAFFT_maxiterate} eq "") {print PARAM "\$MAFFT_MAX_ITERATES=0;\n";}
		else {print PARAM "\$MAFFT_MAX_ITERATES=$MAFFT_max_iterates{$FORM{MAFFT_maxiterate}};\n";}
		print PARAM "\$MAFFT_REFINE=0;\n" if ($FORM{MAFFT_refinement} eq "");
		print PARAM "\$MAFFT_REFINE=1;\n" if ($FORM{MAFFT_refinement} eq "localpair");
		print PARAM "\$MAFFT_REFINE=2;\n" if ($FORM{MAFFT_refinement} eq "genafpair");
		print PARAM "\$MAFFT_REFINE=3;\n" if ($FORM{MAFFT_refinement} eq "globalpair");
	}
	else
	{
	    print PARAM "\$MAFFT_MAX_ITERATES=0;\n";
	    print PARAM "\$MAFFT_REFINE=0;\n";
	}
	if ($FORM{MSA_Program} eq "PRANK") {
		print PARAM "\$PRANK_INSERTION=0;\n" if ($FORM{PRANK_F} eq "+F");
		print PARAM "\$PRANK_INSERTION=1;\n" if ($FORM{PRANK_F} eq "-F");
	}
	else
	{
		print PARAM "\$PRANK_INSERTION=0;\n";
	}
	print PARAM "\$Bootstraps=$FORM{Bootstraps};\n";
	print PARAM "\$SP_COL_CUTOFF=$VARS{SP_COL_CUTOFF};\n";
	print PARAM "\$SP_SEQ_CUTOFF=$VARS{SP_SEQ_CUTOFF};\n";
	if ($FORM{Align_Order} eq "as_input") {print PARAM "\$Align_Order=0;\n";}
	elsif ($FORM{Align_Order} eq "aligned") {print PARAM "\$Align_Order=1;\n";}
	
	if ($FORM{Seq_Type} eq "AminoAcids") {
		print PARAM "\$Seq_Type=0;\n";
		print PARAM "\$CodonTable=0;\n"; # DEFAULT VALUE
	}
	elsif ($FORM{Seq_Type} eq "Nucleotides") {
		print PARAM "\$CodonTable=0;\n"; # DEFAULT VALUE
		print PARAM "\$Seq_Type=1;\n";
	}
	elsif ($FORM{Seq_Type} eq "Codons") 
	{
		print PARAM "\$Seq_Type=2;\n";
		if ($FORM{CodonTable}==1){print PARAM "\$CodonTable=0;\n";}
		elsif ($FORM{CodonTable}==15){print PARAM "\$CodonTable=1;\n";}
		elsif ($FORM{CodonTable}==6){print PARAM "\$CodonTable=2;\n";}
		elsif ($FORM{CodonTable}==10){print PARAM "\$CodonTable=3;\n";}
		elsif ($FORM{CodonTable}==2){print PARAM "\$CodonTable=4;\n";}
		elsif ($FORM{CodonTable}==5){print PARAM "\$CodonTable=5;\n";}
		elsif ($FORM{CodonTable}==3){print PARAM "\$CodonTable=6;\n";}
		elsif ($FORM{CodonTable}==13){print PARAM "\$CodonTable=7;\n";}
		elsif ($FORM{CodonTable}==9){print PARAM "\$CodonTable=8;\n";}
		elsif ($FORM{CodonTable}==14){print PARAM "\$CodonTable=9;\n";}
		elsif ($FORM{CodonTable}==4){print PARAM "\$CodonTable=10;\n";}
	}
	print PARAM "?>\n";
}
sub count_lines
{
    my $file=shift;
    print "FILE:NAME:$file\n";
    my $count=-1;
    open (IN,$file) || return $count;
    while (my $line=<IN>)
	{
		$count++;
	}
    close (IN);
    return $count;
}
sub store_data{
    print LOG "store_data : storing hashes to files $stored_data_file\n";
    store \%VARS, "$stored_data_file";
    chmod 0600, "$stored_data_file";
}

sub extract_seq_from_MSA{
        my $inMSA=shift;
        my $SeqFile=shift;
        my $line_num=0;

        open (IN,$inMSA) || exit_on_error('sys_error', "extract_seq_from_MSA:Can't open '$inMSA': $!");
        open (OUT,">$SeqFile") || exit_on_error('sys_error', "extract_seq_from_MSA: Can't open seqs File: '$SeqFile': $!");
        my $seqNum=1;
        while (my $line=<IN>)
        {
                $line_num++;
                if ($line!~/>/)
                {
                        chomp ($line);
                        $line=~s/-//g;
                        print OUT $line;
                }
                else
                {
                        if ($seqNum>1){print OUT "\n$line";}
                        else {print OUT "$line";}
                        $seqNum++;
                }
        }
        close (OUT);
        close (IN);
}

sub names_according_CoS{
	my $MSA=shift;
	open (MSA,"$MSA");
	my @MSA=<MSA>;
	close (MSA);
	open (MSA,">$MSA");
	foreach my $line (@MSA)
    {
		if ($line=~/^>([0-9]+)/)
		{
			my $seq_num=$1;
			$seq_num=$seq_num-1;
			if ($seq_num<10)
			{
				print MSA ">seq000$seq_num\n";
			}
			elsif ($seq_num<100)
			{
				print MSA ">seq00$seq_num\n";
			}
			elsif ($seq_num<1000)
			{
				print MSA ">seq0$seq_num\n";
			}
			else
			{
				print MSA ">seq$seq_num\n";
			}
		}
		else
		{
			print MSA $line;
		}
    }
	close (MSA);
}

sub MSA_Depth{
  my $inMSA=shift;
  my $in  = Bio::AlignIO->new( '-format' => 'fasta' , -file => $inMSA) or exit_on_error ("sys_error","MSA_Depth: Can't open MSA: '$inMSA' $!"); # die "Can't open $inMSA: $!";
  my $aln = $in->next_aln;
  my $MSA_Depth=0;
  foreach my $seq ($aln->each_seq)
  {
      $MSA_Depth++;
  }
  return $MSA_Depth;
}

sub Round_Scores_File
{
	my $Score_File=shift;
	open (IN,$Score_File);
	my @In=<IN>;
	close (IN);
	open (OUT,">$Score_File");
	foreach my $line (@In)
	{
		$line=trim($line);
		my @line=split(/\s+/,$line);
		if ($line[1]=~m/[0-9\.]+/)
		{
			print OUT "$line[0]\t",sprintf("%.3f", $line[1]),"\n";
		}
		else
		{
			print OUT "$line[0]\t$line[1]\n";
		}
	}
	close (OUT);
}
sub extract_MEAN_RES_PAIR_SCORE
{
	my $MSA_score_file=shift;
	open (MSA_SCR,$MSA_score_file);
	my $MEAN_RES_PAIR_SCORE=0;
	while (my $line=<MSA_SCR>)
	{
		if ($line=~/^\#MEAN_RES_PAIR_SCORE ([0-9.]+)/)
		{
			$MEAN_RES_PAIR_SCORE=$1;
		}
	}
	close (MSA_SCR);
	return $MEAN_RES_PAIR_SCORE;
}
sub count_files_on_dir
{
	my $dir=shift; 	# Path to directory you want to count files in
	my $file_count=0;
	opendir(DIR, $dir) || die "Can't open DIR:'$dir' $!\n";
	while(my $FILE = readdir(DIR)) 
	{
		if($FILE =~ /^\.\.?/)
		{
			next;
		}
		elsif ((-e "$dir$FILE") and (!-d "$dir$FILE") and (-s "$dir$FILE">0)) {
			$file_count++;
		}
	}
	closedir(DIR);
	return $file_count;
}

sub convert_fs_to_lower_case
{
	my $File=shift;
	open (FILE,$File) || return "convert_fs_to_lower_case: Fail to open '$File' for reading: $!";
	my @file=<FILE>;
	close (FILE);
	open (OUT_FILE,">$File") || return "convert_fs_to_lower_case: Fail to open '$File' for writing: $!";;
	foreach my $line (@file)
	{
		if ($line=~/^>/)
		{
			print OUT_FILE $line;
		}
		else
		{
			print OUT_FILE lc($line);
		}
	}
	close (OUT_FILE);
}
sub convert_fs_to_upper_case
{
	my $File=shift;
	open (FILE,$File) || return "convert_fs_to_upper_case: Fail to open '$File' for reading: $!";
	my @file=<FILE>;
	close (FILE);
	open (OUT_FILE,">$File") || return "convert_fs_to_upper_case: Fail to open '$File' for writing: $!";;
	foreach my $line (@file)
	{
		if ($line=~/^>/)
		{
			print OUT_FILE $line;
		}
		else
		{
			print OUT_FILE uc($line);
		}
	}
	close (OUT_FILE);
}

sub Codon_Aln_to_AA_Aln
{
	my $codon_aln=shift;
	my $AA_aln=shift;
	my $codon_table=shift;
	my $xcodon_MSAfile=shift;
	open (CODON_ALN,$codon_aln) || return ('sys_error', "Codon_Aln_to_AA_Aln:Can't open '$codon_aln': $!");
	open (OUT,">$AA_aln")       || return ('sys_error', "Codon_Aln_to_AA_Aln:Can't open '$AA_aln': $!");
	my $seq="";
	my $seq_name="";
	my $Warnning="";
	my $x_flag;
	my $StopCodon_Warnning="";
	while (my $line=<CODON_ALN>)
	{
		chomp ($line);
		$line=~ s/^\s+|\s+$//g;
		if (($line!~/>/)and ($line ne ""))
		{
			$seq=$seq.$line;
		}
		elsif ($line=~/^>(.*)/)
		{
			# translate prev seq
			if (($seq eq "") and ($seq_name ne ""))
			{
				return ("The sequence named '$seq_name' is missing<br>");
			}
			if (($seq ne "") and ($seq_name ne ""))
			{
				my $AA_Seq="";
				my $StopCodon_in_Seq="";
				($x_flag,$AA_Seq,$StopCodon_in_Seq)=translate_sequence($seq,$seq_name,$codon_table,"no","$VARS{WorkingDir}$xcodon_MSAfile");
				if ($AA_Seq ne "")
				{
					print OUT ">$seq_name\n",$AA_Seq,"\n";
				}
				else
				{
					return ("Seq:$seq_name is empty or without any leagal codons");
				}
				$Warnning="Unknown codons on the alignment file were treated as 'X' see more detailes <A href=\"$xcodon_MSAfile\">here</A>\n" if ($x_flag eq "yes");
				#$StopCodon_Warnning=$StopCodon_Warnning.",$StopCodon_in_Seq" if($StopCodon_in_Seq ne "");
				#$StopCodon_Warnning="Stop codons were removed from all the sequences" if($StopCodon_in_Seq ne "");
			} 
			# Start new seq
			if ($line=~/^>(.*)/)
			{
				$seq_name=$1;
				$seq_name=~ s/^\s+|\s+$//g ;
				$seq="";
			}
		}
	}
	# validate last seq
	if (($seq eq "") and ($seq_name ne ""))
	{
		return ("The sequence named '$seq_name' is missing<br>");
	}
	else
	{
		my $AA_Seq="";
		my $StopCodon_in_Seq="";
		($x_flag,$AA_Seq,$StopCodon_in_Seq)=translate_sequence($seq,$seq_name,$codon_table,"no","$VARS{WorkingDir}$xcodon_MSAfile");
		if ($AA_Seq ne "")
		{
			print OUT ">$seq_name\n",$AA_Seq,"\n";
		}
		else
		{
			return ("Seq:$seq_name is empty or without any leagal codons");
		}
		$Warnning="Unknown codons on the alignment file were treated as 'X' see more detailes <A href=\"$xcodon_MSAfile\">here</A>\n" if ($x_flag eq "yes");
		#	$StopCodon_Warnning=$StopCodon_Warnning.",$StopCodon_in_Seq" if($StopCodon_in_Seq ne "");
		#   $StopCodon_Warnning="Stop codons were removed from all the sequences" if($StopCodon_in_Seq ne "");
	}
	close (OUT);
	close (CODON_ALN);
	$Warnning=$Warnning.$StopCodon_Warnning;
	return ("OK",$Warnning);
}

sub translate_sequence{

    my $DNASequence = shift;
    my $DNASequenceName = shift;
    my $codonTableIndex = shift;
    my $xFlag = shift;
    my $xCodonfile = shift;
    my ($codon, $AA,$StopCodon_in_Seq);

    # The actual translation
    my $codonTable_obj  = Bio::Tools::CodonTable -> new ( -id => $codonTableIndex );
    my $seq_length = length($DNASequence);
    my $i =0;
    my $AASeq = "";


    while ($i<$seq_length-2){
        $codon = substr($DNASequence, $i, 3);
        if ($codon eq '---'){
            $AA = '-';
        }
        else{
            $AA = $codonTable_obj->translate($codon);
            # if the AA is X we print the codon to a file and later inform the user
            if ($AA eq "X"){
                $xFlag = "yes";
                open X_CODON, ">>$xCodonfile";
                print X_CODON "<tr><td>$DNASequenceName<\/td><td>".($i+1)."<\/td><td>$codon<\/td></tr>\n";
                close X_CODON;
            }
        }
        $AASeq.= $AA;
        $i+=3;
    }
#	# if the last charachter in the sequence is '*', we cut it from the sequence. we put a flag, so we can inform the user if needed. # UPDATE: STOP CODONS ARE NOT ALLOWED ANY MORE IN CODON ALIGNMENTS (We give the user to choose if he wants to put a gap or remove it from ALN).
#	if ($AASeq =~ m/\*/){
#		$AASeq=~s/\*//;
#		$StopCodon_in_Seq=$DNASequenceName;
#	}
    return ($xFlag, $AASeq,$StopCodon_in_Seq); 
}

sub fix_MAFFT_RoughTree
{
	my $TreeFile=shift;
	my $TreeOrig=$TreeFile.".orig";
	copy($TreeFile,$TreeOrig);
	open (TREE,$TreeFile) || die "Can't open MAFFT RoughTree:'$TreeFile' for reading: $!";
	my $Tree=<TREE>;
	close (TREE);
	chomp $Tree;
	$Tree=$Tree.";" if ($Tree!~/;$/);
	$Tree=~s/([0-9]+)_([0-9]+)/$1/g;
	$Tree=~s/[_]+//g;
	open (TREE,">$TreeFile") || die "Can't open MAFFT RoughTree:'$TreeFile' for writing: $!";
	print TREE "$Tree\n";
	close (TREE);
}
sub prune_tree
{
	my $InTree=shift;
	my $TaxToRemove=shift;
	my $OutTree=shift;
	my $prune_tree_prog=shift;
	system ("$prune_tree_prog $InTree $TaxToRemove $OutTree");
}

# Returns the filename without directory
sub getFilename{
        my $fullFile = pop @_;
        if ($fullFile =~ m/.*[\\\/](.*)$/) {
                return $1;
        } else {return $fullFile}

}
# Returns the basename of the file (no directory)
sub getBasename{
        my $fullFile = pop @_;
        my $filename = getFilename($fullFile);
        if ($filename =~ m/(.*)\..*$/) {
                return $1;
        } else {return $filename}

}
