use lib "/bioseq/pupkoSVN/trunk/www/bioSequence_scripts_and_constants/";
use lib "/bioseq/bioSequence_scripts_and_constants/";

use GENERAL_CONSTANTS;
use cp_rasmol_gradesPE_and_pipe;
my $PDB_ID=shift;
my $chain=shift;
my $zipped_PDB_File=shift;
my $Grades_File=shift;

my $Output_Dir=GENERAL_CONSTANTS::BIOSEQ_TEMP;
my $out=$Output_Dir.$PDB_ID."_".$chain."_With_ConSurf.pdb";
my $out_isd=$Output_Dir.$PDB_ID."_".$chain."_With_ConSurf_isd.pdb";

my $Unzipped_PDB=$Output_Dir.$PDB_ID."_".$chain;
if (-e $zipped_PDB_File)
{
	system ("zcat $zipped_PDB_File > $Unzipped_PDB");
#	my $ans=cp_rasmol_gradesPE_and_pipe::replace_tempFactor($Unzipped_PDB,$chain,\%Rate4Site_Grades,$out1);
	print  "Calling: cp_rasmol_gradesPE_and_pipe::ReplaceTempFactConSurfScore($chain,$Unzipped_PDB,$Grades_File,$out,$out_isd);\n";
	my @ans=cp_rasmol_gradesPE_and_pipe::ReplaceTempFactConSurfScore($chain,$Unzipped_PDB,$Grades_File,$out,$out_isd);
	if ($ans[0] eq "OK") 
	{
		print "$out"
	}
	else
	{
		print "FAILED";
	}
	
	#print "cp_rasmol_gradesPE_and_pipe::ReplaceTempFactConSurf($chain,$Unzipped_PDB,$Grades_File,$out1,$out2);";
	#my $ans=cp_rasmol_gradesPE_and_pipe::ReplaceTempFactConSurf($chain,$Unzipped_PDB,$Grades_File,$out1,$out2);
	#print "$ans";
}

