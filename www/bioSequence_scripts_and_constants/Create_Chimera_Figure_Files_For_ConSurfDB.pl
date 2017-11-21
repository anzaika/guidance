use lib "/bioseq/pupkoSVN/trunk/www/bioSequence_scripts_and_constants/";
use lib "/bioseq/bioSequence_scripts_and_constants/";

use GENERAL_CONSTANTS;
use cp_rasmol_gradesPE_and_pipe;
my $PDB_ID=shift;
my $chain=shift;
my $zipped_PDB_File=shift;
my $Grades_File=shift;
my $OUT_DIR_URL=shift;

my $Output_Dir=GENERAL_CONSTANTS::BIOSEQ_TEMP;
my $ATOMS_with_ConSurf_Scores=$Output_Dir.$PDB_ID."_".$chain."_With_ConSurf.pdb";
my $ATOMS_with_ConSurf_Scores_isd=$Output_Dir.$PDB_ID."_".$chain."_With_ConSurf_isd.pdb";
my $ATOMS_with_ConSurf_Scores_no_Path=$PDB_ID."_".$chain."_With_ConSurf.pdb";
my $ATOMS_with_ConSurf_Scores_isd_no_Path=$PDB_ID."_".$chain."_With_ConSurf_isd.pdb";
my $chimerax_script_for_figure=$Output_Dir.$PDB_ID."_".$chain."_Figure.chimerax";
my $chimerax_script_for_figure_isd=$Output_Dir.$PDB_ID."_".$chain."_Figure_isd.chimerax";


# Check For insufficient_data
my $insufficient_data="no"; # Are There insuficient data?
if (is_ther_insufficient_data_ConSurf_gradesPE($Grades_File) eq "yes")
{
	$insufficient_data="yes";
}

my $Unzipped_PDB=$Output_Dir.$PDB_ID."_".$chain;
if (-e $zipped_PDB_File)
{
	system ("zcat $zipped_PDB_File > $Unzipped_PDB");
#	my $ans=cp_rasmol_gradesPE_and_pipe::replace_tempFactor($Unzipped_PDB,$chain,\%Rate4Site_Grades,$out1);
	my @ans;
	$ans[0]="OK";
	if ((!-e $ATOMS_with_ConSurf_Scores) and !-e ($ATOMS_with_ConSurf_Scores_isd))
	{
		print  "Calling: cp_rasmol_gradesPE_and_pipe::ReplaceTempFactConSurfScore($chain,$Unzipped_PDB,$Grades_File,$ATOMS_with_ConSurf_Scores,$ATOMS_with_ConSurf_Scores_isd);\n";
		@ans=cp_rasmol_gradesPE_and_pipe::ReplaceTempFactConSurfScore($chain,$Unzipped_PDB,$Grades_File,$ATOMS_with_ConSurf_Scores,$ATOMS_with_ConSurf_Scores_isd);
	}
	if ($ans[0] eq "OK") 
	{
		print "Calling: cp_rasmol_gradesPE_and_pipe::create_chimera_image_script($chimerax_script_for_figure,$ATOMS_with_ConSurf_Scores_no_Path,$OUT_DIR_URL)\n";
        	@ans=cp_rasmol_gradesPE_and_pipe::create_chimera_image_script($chimerax_script_for_figure,$ATOMS_with_ConSurf_Scores_no_Path,$OUT_DIR_URL);
        	if ($insufficient_data eq "yes")
        	{
	                print "Calling: cp_rasmol_gradesPE_and_pipe::create_chimera_image_script($chimerax_script_for_figure_isd,$ATOMS_with_ConSurf_Scores_isd_no_Path,$OUT_DIR_URL)\n";
                	@ans=cp_rasmol_gradesPE_and_pipe::create_chimera_image_script($chimerax_script_for_figure_isd,$ATOMS_with_ConSurf_Scores_isd_no_Path,$OUT_DIR_URL);
        	}
		if ($ans[0] ne "OK") 
		{
			print "FAILED";
		}
	}
	else
	{
		print "FAILED";
	}
	
	#print "cp_rasmol_gradesPE_and_pipe::ReplaceTempFactConSurf($chain,$Unzipped_PDB,$Grades_File,$out1,$out2);";
	#my $ans=cp_rasmol_gradesPE_and_pipe::ReplaceTempFactConSurf($chain,$Unzipped_PDB,$Grades_File,$out1,$out2);
	#print "$ans";
}


sub is_ther_insufficient_data_ConSurf_gradesPE{
# the routine matches each position in the gradesPE file its ConSurf grade. In case there was a grade mark with *, we put it in a seperate hash with the grade 0.
# the routine returns "yes" if a * was found and "no" otherwise
    my $gradesPE_file = shift;
    my $insufficient = "no";
    
    open GRADES, $gradesPE_file;
    while (<GRADES>){
        if (/^\s*\d+\s+\w/ ){
		    my @grades=split;            
		    $grades[2] =~ s/[a-z\:]//gi;
            if ($grades[4] =~/\d\*?/){
		# if it is insufficient color
		if ($grades[4] =~/(\d)\*/){
			$insufficient = "yes";
			last;
		}
            }
        }
    }
    close GRADES; 
    return $insufficient;
}
