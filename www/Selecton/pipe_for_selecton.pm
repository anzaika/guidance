#!/usr/local/bin/perl -w

package pipe_for_selecton;

use lib "/bioseq/bioSequence_scripts_and_constants";
use design_pipe;
use strict;

sub create_pipe{
    my $run_number = shift;
    my $log_file = shift;
    my $WorkingDir = shift;
    my $rasmol_file = shift;
    my $pipe_file = shift;
    my $pdb_input_file = shift;
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    $year += 1900;
    $mon += 1;   
    my $run_date = $year . '-' . $mon . '-' . $mday;   
    my $run_begin_time;  # will be read from the LOG file, after the words: 'Begin time: '
    my $run_end_time = $hour . ':' . $min . ':' . $sec;
    my $pdb_id = shift; #  #FORM  should be id code or "UPLOADED";
    my $chain = shift; #  #FORM  blank for name "none"
    my $identical_chains = ""; #blank for name "none"
    
    my $DNA_Filename = shift;  #FORM : since this should be the user's input, I think I should save the user's file name as recived in the cgi. it is saved in the var $upload_unaligned_file_dna or $upload_MSA_file_dna
    
    my $query_seq = shift;  #FORM : if not given in the form, the program chooses it: than write the sequence on which the anlysis was made
    
    my $Model = shift; #FORM  Evolutionary Model : pointer to 1 cell array which holds the string that describes the model
    my $distribute_categories = shift; #FORM  a number between 4 to 14
    my $optimizeBL = shift;     #FORM  change to TRUE if "y", FALSE if "n"
    my $genetic_code = shift;  #FORM   : pointer to 1 cell array which holds the string that describes the genetic code
    my $pipe_error  = shift;
    my $pipe_empirical = shift;
    my $pipe_precision = shift;

        
    my $header_found = "no";
    my $header_line;
    my $title_found = "no";
    my $compnd_found = "no";
    my $title_lines = "";
    my $compnd_lines = "";
    
    
    open LOG, ">>".$log_file;
    
    unless (open ERROR, ">".$WorkingDir.$pipe_error){
	print LOG "pipe_for_selecton.pl : could not open $WorkingDir"."$pipe_error file for writing $!\n";
	close LOG;
	exit;
    }
    close LOG;
    
    unless (open PDB, $pdb_input_file){
	&exitOnError("pipe_for_selecton::can't open file $pdb_input_file for reading.\n", $pipe_error);
    }
    
    while (<PDB>)
    {
       if(/^HEADER/)   
       {
	    $_ =~ s/\s+$//;
	    $header_found = "yes";
	    $header_line = $_."\n";
       }
       elsif (/^TITLE\s+/)
	{
	    $_ =~ s/\s+$//;
	    $_ =~ /^TITLE\s+\d*\s(.*)\s*/;
	    $title_lines.=$1;
	    $title_found = "yes";
	}
	elsif (/^COMPND\s+\d*\s(.*)/)
	{
	    $compnd_lines.=$1;
	    $compnd_found = "yes";
	}
	elsif (/^SOURCE/ || /^KEYWDS/ || /^AUTHOR/ || /^SEQRES/ || /^ATOM/) # no nead to go over all the pdb
	{
	  last;
	}
    }
    close PDB;
    
    unless (open LOG, "<$log_file")
    {
	&exitOnError("pipe_for_selecton:: Can't open $log_file for reading.\n", $pipe_error);
    }
    
    while(<LOG>){
	if(/Begin time: (\d+:\d+:\d+),\s.*/)   {
	    $run_begin_time = $1;
	    last;
       }
    }
    close LOG;
    
    unless (open PIPE, ">".$WorkingDir.$pipe_file){
       &exitOnError("pipe_for_selecton:: can't open file $WorkingDir"."$pipe_file for writing\n", $pipe_error);
    }
    
    if ($header_found eq "yes"){
       print PIPE $header_line;
    }
    else{
       print PIPE "HEADER                                 [THIS LINE ADDED FOR JMOL COMPATIBILITY]\n";
    }
    
    print PIPE "!! ====== IDENTIFICATION SECTION ======\n";
    #This section defines javascript variables that describe the Selecton job.
    print PIPE "!js.init\n";
    print PIPE "! selecton_version = \"2.2\";\n";
    print PIPE "! selecton_run_number = \"$run_number\";\n";
    print PIPE "! selecton_run_date = \"$run_date\";\n"; 
    print PIPE "! selecton_run_submission_time = \"$run_begin_time\";\n";
    print PIPE "! selecton_run_completion_time = \"$run_end_time\";\n";
    print PIPE "!\n";
    if ($DNA_Filename =~ m/^SELECTON_UN(.+)/){
		print PIPE "! selecton_dna_filename = \"$1\";\n";
		print PIPE "! selecton_codon_aligned_filename = \"\";\n";
    }
    elsif ($DNA_Filename =~ m/^SELECTON_MSA(.+)/){
		print PIPE "! selecton_dna_filename = \"\";\n";
		print PIPE "! selecton_codon_aligned_filename = \"$1\";\n";
    }
    else{
		print PIPE "! selecton_dna_filename = \"\";\n";
		print PIPE "! selecton_codon_aligned_filename = \"\";\n";
    }
    print PIPE "! selecton_query_sequence_name = \"$query_seq\";\n";
    print PIPE "!\n";
    print PIPE "! selecton_pdb_id = \"$pdb_id\";\n";
    ($chain =~ /none/i) ? print PIPE "! selecton_chain = \"\";\n" : print PIPE "! selecton_chain = \"$chain\";\n";
    print PIPE "! selecton_identical_chains = \"$identical_chains\";\n";
    print PIPE "!\n";
    print PIPE "! selecton_model = \"".$Model->[0]."\";\n";
    print PIPE "! selecton_aa_matrix = \"$pipe_empirical\";\n";
    print PIPE "! selecton_distribution_categories = \"$distribute_categories\";\n";
    print PIPE "! selecton_tree_filename = \"\";\n";
    print PIPE "! selecton_optimize_branch_lengths = \"$optimizeBL\";\n";
    print PIPE "! selecton_genetic_code = \"".$genetic_code->[0]."\";\n";
    print PIPE "! selecton_precision_level = \"$pipe_precision\";\n";

    print PIPE "!\n";
    print PIPE "!! ====== CONTROL PANEL OPTIONS SECTION ======\n";
    print PIPE "!js.init\n";
    print PIPE "! pipe_title = \"<i>Selecton View:</i> $pdb_id chain $chain.\"\n";
    print PIPE "!! pipe_subtitle is from TITLE else COMPND\n";
    print PIPE "!!\n";
    print PIPE "! pipe_subtitle = \n";
    if ($title_found eq "yes"){
	print PIPE (design_pipe::print_string_max80($title_lines, "!"))."\";\n";
    }
    elsif ($compnd_found eq "yes"){
        print PIPE (design_pipe::print_string_max80($compnd_lines, "!"))."\";\n";
    }
    print PIPE "! pipe_title_enlarged = false;\n";
    print PIPE "! pipe_background_color = \"white\";\n";
    print PIPE "!\n";
    print PIPE "!! Specify the ConSurf control panel (also used for Selecton)\n";
    print PIPE "!!\n";
    print PIPE "! pipe_cp1 = \"consurf\/consurf.htm\";\n";
    print PIPE "!\n";
    print PIPE "!! If you want the frontispiece to be reset every time you enter this\n";
    print PIPE "!! page, use false. If this is a one-page presentation (no contents)\n";
    print PIPE "!! and you want to be able to return from QuickViews without resetting\n";
    print PIPE "!! the view, use true.\n";
    print PIPE "!!\n";
    print PIPE "! frontispiece_conditional_on_return = true;\n";
    print PIPE "!\n";
    print PIPE "!! Open the command input slot\/message box to 30% of window height.\n";
    print PIPE "!!\n";
    print PIPE "! pipe_show_commands = true;\n";
    print PIPE "! pipe_show_commands_pct = 30;\n";
    print PIPE "!\n";
    print PIPE "!! Do not show the PiPE presentation controls in the lower left frame.\n";
    print PIPE "!!\n";
    print PIPE "! pipe_hide_controls = true;\n";
    print PIPE "!\n";
    print PIPE "!! Hide development viewing mode links at the bottom of the control panel.\n";
    print PIPE "!!\n";
    print PIPE "! pipe_tech_info = false;\n";
    print PIPE "!\n";
    print PIPE "!! pipe_start_spinning = true; \/\/ default is PE\'s Preference setting.\n";
    print PIPE "!! top.nonStopSpin = true; \/\/ default: spinning stops after 3 min.\n";
    print PIPE "!!\n";
    print PIPE "!! ====== COLORS SECTION ======\n";
    print PIPE "!!\n";
    print PIPE "!color color_carbon C8C8C8\n";
    print PIPE "!color color_sulfur FFC832\n";
    print PIPE "!\n";
    print PIPE "!! Seven Selecton color grades follow:\n";
    print PIPE "!!\n";
    print PIPE "!color color_grade0 FFFF96 insufficient data\n";
    print PIPE "!color color_grade1 FFBD00\n";
    print PIPE "!color color_grade2 FFFF78\n";
    print PIPE "!color color_grade3 FFFFFF\n";
    print PIPE "!color color_grade4 FCEDF4\n";
    print PIPE "!color color_grade5 FAC9DE\n";
    print PIPE "!color color_grade6 F07DAB\n";
    print PIPE "!color color_grade7 A02560\n";
    print PIPE "!\n";
    print PIPE "!\n";
    my $chain_to_send;
    if ($chain =~ /none/i){  
        $chain_to_send = "protein\n";}
    else{
        $chain_to_send = ":$chain\n";}
    
    &read_rasmol($WorkingDir.$rasmol_file, $chain_to_send);
        
    print PIPE "!! ====== END OF selecton PiPE BLOCK ======\n";    
    
    
    unless (open PDB, $pdb_input_file){
	&exitOnError("pipe_for_selecton::can't open file $pdb_input_file for reading.\n", $pipe_error);
    }
    while (<PDB>){   # We copy all lines from original PDB, apart from the header line and lines that don't unclude info
	if(/^HEADER/) {}
        elsif(/^\s+$/) {}
	else{
	    print PIPE $_;
	}
    }
    close PDB;    
    close PIPE;
    close ERROR;
}

#--------------------------------------------
# the routine which extract the information from rasmol script and prints it to the pipe in the desired format
# the routine prints two scripts: one for isd - InSufficient Data - and the other for classic coloring.
# REMARK: when the InSufficient Data part will be implemented - the color: grade0 should be added. currently the colors are marked grade1-grade7
sub read_rasmol{
    my $rasmol_file = shift;
    my $chain = shift;
    my $last_char;
    my @residues;   # array to hold splitted lines
    my $all_residues = "";  # line that will hold all the residues for each color
    my $no_of_color_res = 0;
    my @AoA = ();	
        $AoA[0] = [0, ""];  # Array of 8 cells which hold : AoA[i][0]:number of residues for color i, AoA[i][1]:names of residues for color i
	my $i;
	for ($i=0; $i<7; $i++)
		{$AoA[$i] = ["",""];}	
        
    open RASMOL, $rasmol_file or die "cannot open file $rasmol_file $!";
    while (<RASMOL>){
        chomp;
        if (/select\s+[A-Za-z]{3}\d+/ or /select selected or/){
            @residues = split(' ', $_);
            foreach (@residues){
                if (/([A-Za-z]{3}\d+):?\w?/){
                    $no_of_color_res++;
                    $all_residues.="$1, ";
                }
            }
        }    
        elsif (/color \[160,37,96\]/){
            $AoA[7] = [$no_of_color_res, $all_residues];
        }
        elsif (/color \[240,125,171\]/){
            $AoA[6] = [$no_of_color_res, $all_residues];
        }
        elsif (/color \[250,201,222\]/){
            $AoA[5] = [$no_of_color_res, $all_residues];
        }
        elsif (/color \[252,237,244\]/){
            $AoA[4] = [$no_of_color_res, $all_residues];
        }
        elsif (/color \[255,255,255\]/){
            $AoA[3] = [$no_of_color_res, $all_residues];
        }
        elsif (/color \[255,255,120\]/){
            $AoA[2] = [$no_of_color_res, $all_residues];
        }
        elsif (/color \[255,190,0\]/){
            $AoA[1] = [$no_of_color_res, $all_residues];
        }    
        elsif (/spacefill/){
            $all_residues = "";
            $no_of_color_res = 0;
        }
        #print $all_residues."\n";
    }
    close RASMOL;
    $i = 0;
    my $color_string_to_print = "";
    for ($i=0; $i<7; $i++){
        if ($AoA[$i][0] eq "") {
            $color_string_to_print.= "0,";}
        else{
            $color_string_to_print.= $AoA[$i][0].",";}
    }
    if ($AoA[7][0] eq ""){
        $color_string_to_print.= "0\);\n";}
    else{
        $color_string_to_print.= $AoA[7][0]."\);\n";} #last number will not come with ","
    
    print PIPE "!! ====== RESULTS SECTION ======\n";
    print PIPE "!!----------------------------------------\n";
    print PIPE "!!\n";
    print PIPE "!js.init\n";
    print PIPE "! selecton_grade_freqs_isd = Array\(".$color_string_to_print;
    print PIPE "! selecton_grade_freqs = Array\(".$color_string_to_print;
    print PIPE "!\n";
    print PIPE "!! ====== SCRIPTS SECTION ======\n";
    print PIPE "!!----------------------------------------\n";
    print PIPE "!!\n";
    print PIPE "!spt \#name=select_and_chain\n";
    print PIPE "! select selected and $chain";
    print PIPE "!\n";
    print PIPE "!!----------------------------------------\n";
    print PIPE "!!\n";
    print PIPE "!spt \#name=view01\n";
    print PIPE "! \@spt selecton_view_isd\n";
    print PIPE "!\n";
    print PIPE "!!----------------------------------------\n";
    print PIPE "!!\n";
    print PIPE "!spt \#name=hide_all\n";
    print PIPE "! restrict none\n";
    print PIPE "! ssbonds off\n";
    print PIPE "! hbonds off\n";
    print PIPE "! dots off\n";
    print PIPE "! list * delete\n";
    print PIPE "!\n";
    print PIPE "!!----------------------------------------\n";
    print PIPE "!! common_spt uses CPK carbon gray (or phosphorus yellow) for backbones.\n";
    print PIPE "!!\n";
    print PIPE "!spt \#name=common_spt\n";
    print PIPE "! \@spt hide_all\n";
    print PIPE "! select all\n";
    print PIPE "! color [xC8C8C8] \# rasmol/chime carbon gray\n";
    print PIPE "! select nucleic\n";
    print PIPE "! color [xFFA500] \# phosphorus orange\n";
    print PIPE "! select hetero\n";
    print PIPE "! color cpk\n";
    print PIPE "! select not hetero\n";
    print PIPE "! backbone 0.4\n";
    print PIPE "! javascript top.water=0\n";
    print PIPE "! \n";
    print PIPE "! ssbonds 0.3\n";
    print PIPE "! set ssbonds backbone\n";
    print PIPE "! color ssbonds \@color_sulfur\n";
    print PIPE "! \n";
    print PIPE "! select hetero and not water\n";
    print PIPE "! spacefill 0.45\n";
    print PIPE "! wireframe 0.15\n";
    print PIPE "! dots 50\n";
    print PIPE "! \n";
    print PIPE "! select protein\n";
    print PIPE "! center selected\n";
    print PIPE "! \n";
    print PIPE "!!----------------------------------------\n";
    print PIPE "!!\n";
    print PIPE "!spt \#name=selecton_view_isd\n";
    print PIPE "! \@spt common_spt\n";
    print PIPE "! \@for \$=0, 7\n";
    print PIPE "! \@spt select_isd_grade\$\n";
    print PIPE "! \@spt select_and_chain\n";
    print PIPE "! color \@color_grade\$\n";
    print PIPE "! spacefill\n";
    print PIPE "! \@endfor\n";
    print PIPE "! zoom 115\n";
    print PIPE "!\n";
    for ($i=0; $i<8; $i++){
        print PIPE "!!----------------------------------------\n";
        print PIPE "!!\n";
        print PIPE "!spt \#name=select_isd_grade".$i."\n";
        print PIPE "!\n";
        if ($AoA[$i][1] eq ""){
            print PIPE "! select none\n";
        }
        else{
            print PIPE &print_string_max80($AoA[$i][1]);
            print PIPE "\n";
        }
        print PIPE "!\n!\n";
    }
    for ($i=0; $i<8; $i++){
        print PIPE "!!----------------------------------------\n";
        print PIPE "!!\n";
        print PIPE "!spt \#name=select_grade".$i."\n";
        print PIPE "!\n";
        if ($AoA[$i][1] eq ""){
            print PIPE "! select none\n";
        }
        else{
            print PIPE &print_string_max80($AoA[$i][1]);
            print PIPE "\n";
        }
        print PIPE "!\n!\n";
    }
}
#--------------------------------------------
# edit the lines so it will not exceed 80 chars
sub print_string_max80($)
{
    my $input_string = shift; #"Ala24:H, Ser25:H, Tyr27:H, Thr30:H, Tyr32:H, Thr74:H, Ser75:H, Lys76:H, Ser77:H";
    
    my $final_string;
    my @string;
    my $line_length=0;
    my $last_char;
    

    @string = split(/ /,$input_string);
    $final_string = "! select ";
    $line_length = length $final_string;
    foreach (@string)
    {
        $line_length += ((length $_)+1); # add 1, since it will be the length of the added " "
        if ($line_length > 80)
        {
            chop($final_string);
            chop($final_string);
            $final_string.= "\n! select selected or $_ ";
            $line_length = 22 + length $_;
        }
        else
        {
            $final_string.= $_." ";
        }
        
    }
    chop($final_string);
    $last_char = chop($final_string);    
    unless($last_char eq ","){
        $final_string.=$last_char;
    }
    
    return "$final_string";    
}

#------------------------------------------
sub exitOnError()
{
   my $errorMessage = shift;
   my $pipe_error = shift;
   print ERROR $errorMessage;
   close ERROR;
   chmod 0755, $pipe_error;
   exit;
}

1;