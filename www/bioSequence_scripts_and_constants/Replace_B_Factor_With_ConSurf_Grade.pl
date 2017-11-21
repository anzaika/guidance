use lib "/bioseq/pupkoSVN/trunk/www/bioSequence_scripts_and_constants/";
use lib "/bioseq/bioSequence_scripts_and_constants/";

use GENERAL_CONSTANTS;
use cp_rasmol_gradesPE_and_pipe;
my $PDB_ID=shift;
my $chain=shift;
my $zipped_PDB_File=shift;
my $Grades_File=shift;

my $Output_Dir=GENERAL_CONSTANTS::BIOSEQ_TEMP;
my $out1=$Output_Dir.$PDB_ID."_".$chain."_With_Conservation.pdb";

my $Unzipped_PDB=$Output_Dir.$PDB_ID."_".$chain;
if (-e $zipped_PDB_File)
{
	system ("zcat $zipped_PDB_File > $Unzipped_PDB");
	my %Rate4Site_Grades=();
	print "Calling:cp_rasmol_gradesPE_and_pipe::read_Rate4Site_gradesPE($Grades_File,\%Rate4Site_Grades)\n";
	cp_rasmol_gradesPE_and_pipe::read_Rate4Site_gradesPE($Grades_File,\%Rate4Site_Grades);
	print "Calling:cp_rasmol_gradesPE_and_pipe::replace_tempFactor(\$Unzipped_PDB,$chain,\%Rate4Site_Grades,$out1)\n";
	my $ans=cp_rasmol_gradesPE_and_pipe::replace_tempFactor($Unzipped_PDB,$chain,\%Rate4Site_Grades,$out1);
	if ($ans eq "nothing")
	{
		print "$out1";
	}
	else
	{
		print "FAILED";#$ans;
	}
	
	#print "cp_rasmol_gradesPE_and_pipe::ReplaceTempFactConSurf($chain,$Unzipped_PDB,$Grades_File,$out1,$out2);";
	#my $ans=cp_rasmol_gradesPE_and_pipe::ReplaceTempFactConSurf($chain,$Unzipped_PDB,$Grades_File,$out1,$out2);
	#print "$ans";
}
