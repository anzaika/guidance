#!/usr/bin/perl -w

package new_rate4site_routines;

use lib "/bioseq/bioSequence_scripts_and_constants";
use GENERAL_CONSTANTS;
use lib "/bioseq/ConSurf";
use CONSURF_CONSTANTS;

my %ColorScale = (0 => 9,
                1 => 8,
                2 => 7,
                3 => 6,
                4 => 5,
                5 => 4,
                6 => 3,
                7 => 2,
                8 => 1);

####################################################################################
# There are some tests to see if rate4site failed.
# Since we can't trust only one of them, we do all of them. If one of them is tested to be true - than a flag will get TRUE value
# 1. the .res file might be empty.
# 2. if the run failed, it might be written to the log file of r4s.
# 3. in a normal the r4s.log file there will lines that describe the grades. if it fail - we won't see them
# In one of these cases we try to run the slower version of rate4site.
# We output this as a message to the user.
sub check_if_rate4site_failed
{
    my $res_flag = shift; #path of the file, including working dir
    my $r4s_log = shift;    
    my $return_ans = "no";
    my $print_to_html = "";
    my $print_to_log = "";
    my $r4s_process_id = "";
    my $error_found = "no";
    
    if(!-e $res_flag){
        $print_to_log = "new_rate4site_routines::check_if_rate4site_failed : the file $res_flag does not exsits. \n";
        $error_found = "yes";
    }
    elsif (-e $res_flag && -z $res_flag) #1
    {
        $print_to_log = "new_rate4site_routines::check_if_rate4site_failed : the file $res_flag was found to be of size 0. \n";        
        $error_found = "yes";
    }
    if(-e $res_flag){
        unless (open R4S_RES, $res_flag){
            $print_to_log = "new_rate4site_routines::check_if_rate4site_failed : can not open file: $res_flag. aborting.\n";
            $error_found = "yes";
        }
        while(<R4S_RES>){
            if(/In the tree file there is the name: (.+) that is not found in the sequence file/){
                $print_to_html.= "The sequence name $1 was found in the tree file, but was not found in your MSA.<br>\nPlease correct your tree file, so it will include the same names as they appear in the MSA and re-run your query.<br>\n";
                $error_found = "yes";
                $print_to_log = "new_rate4site_routines::check_if_rate4site_failed : sequence name $1 was found in the tree file, was not found in the MSA\n";
                last;
            }
        }
        close R4S_RES;
    }
    if (-e $r4s_log && !(-z $r4s_log)) #2,3
    {
        unless (open R4SLOG, $r4s_log) {
            $print_to_log = "new_rate4site_routines::check_if_rate4site_failed : can not open file: $r4s_log. aborting.\n";
            $error_found = "yes";
        }
        while (<R4SLOG>)
        {
            if (/^.Process_id= (\d+)/){
                $r4s_process_id = $1;
            }
            if ($_ =~ m/likelihood of pos was zero/){
                $print_to_log = "new_rate4site_routines::check_if_rate4site_failed : the line: \"likelihood of pos was zero\" was found in $r4s_log.\n";
                $error_found = "yes";
                last;
            }
            if ($_ =~ m/rate of pos\:\s\d\s=/){ #if we see this line, we change the flag{                
                $return_ans = "no";
                last;
            }
            if($_ =~ m/The amino-acid sequences contained the character: (.+)/){
                $print_to_html .= "The illegal character $1 was found in your MSA. Please make sure your MSA contains only the following characters:<br />\nA, B, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V, X, Z, -";
                $print_to_log = "new_rate4site_routines::check_if_rate4site_failed : illegal character $1 was found in the MSA\n";
                $error_found = "yes";
                last;
            }
            if($_ =~ m/Could not find a sequence that matches the sequence name/){
                my $seq_name = <R4SLOG>;
                $print_to_log = "new_rate4site_routines::check_if_rate4site_failed : the submitted query sequence name $seq_name was not found in MSA";
                $print_to_html .= "The query sequence name $seq_name you have submitted was not found in the uploaded MSA. Please note that the sequence name should be exactly as it is in the MSA file<br>";
                $error_found = "yes";
                last;
            }
            if($_ =~ m/The sequence name: (.+)was found in the tree file but not found in the sequence file/){
                my $seq_name = $1;
                $print_to_log = "new_rate4site_routines::check_if_rate4site_failed : the sequence name $1 was found in the tree file, but not in the MSA";
                $print_to_html .= " The tree file is inconsistant with the uploaded MSA. The sequence: \'$1\' was found in the tree file, but was not found in the MSA.<br>";
                $error_found = "yes";
                last;
            }
            if ($_ =~ m/Bad format in tree file/){
                $print_to_log = "new_rate4site_routines::check_if_rate4site_failed : ";
                $print_to_html .= " There is an error in the tree file format. Please check that your tree is in the <a href = \"".GENERAL_CONSTANTS::CONSURF_TREE_FAQ."\">requested format</a> and reupload it to the server.<br>";
                $error_found = "yes";
                last;
            }
            if ($_ =~ m/not all sequences are of the same lengths/){
                $print_to_log = "new_rate4site_routines::check_if_rate4site_failed : problem with the MSA : not all sequences are of the same lengths";
                $error_found = "yes";
                last;
            }
        }
        close R4SLOG;
    }
    $return_ans = "yes" if ($error_found eq "yes");
    return ($return_ans, $print_to_log, $r4s_process_id, $print_to_html);
}
##------------------------------
##create @output
#sub read_R4S_res{
##------------------------------
#    my ($r4s_res_file, $ref_to_output_array) = @_;
#    my $ret = "OK";
#    # open file
#    unless (open RATE4SITE, $r4s_res_file){
#        $ret = "new_rate4site_routines::read_R4S_res : can't open '$r4s_res_file' for reading $!";
#    }
#    else{
#        my $line;
#        my $i = 0;
#        while (<RATE4SITE>) {
#            $line = $_;
#            chomp $line;
#            # Baysean
#            if ($line =~ /^\s+(\d+)\s+(\w)\s+(\S+)\s+\[\s*(\S+),\s*(\S+)\]\s+\S+\s+(\d+)\/(\d+)/){
#                    $ref_to_output_array->[$i]{POS} = $1;
#                    $ref_to_output_array->[$i]{SEQ} = $2;
#                    $ref_to_output_array->[$i]{GRADE} = $3;
#                    $ref_to_output_array->[$i]{INTERVALLOW} = $4;
#                    $ref_to_output_array->[$i]{INTERVALHIGH} = $5;
#                    $ref_to_output_array->[$i]{MSA_NUM}=$6;
#                    $ref_to_output_array->[$i]{MSA_DENUM}=$7;
#                    $i++;
#            }
#            # Maximum likelihood
#            elsif($line=~m/(\d+)\s+(\w)\s+(\S+)\s+(\d+)\/(\d+)/){
#                $ref_to_output_array->[$i]{POS} = $1;
#                $ref_to_output_array->[$i]{SEQ} = $2;
#                $ref_to_output_array->[$i]{GRADE} = $3;
#                $ref_to_output_array->[$i]{INTERVALLOW} = $3;
#                $ref_to_output_array->[$i]{INTERVALHIGH} = $3;
#                $ref_to_output_array->[$i]{MSA_NUM}=$4;
#                $ref_to_output_array->[$i]{MSA_DENUM}=$5;
#                $i++;
#            }
#        }
#        close RATE4SITE;
#    }    
#    return $ret;
#}
##------------------------------
##calculates layer unit
#sub calc_layer_unit{
##------------------------------
#    my $ref_to_output = shift;
#    my @Output = @$ref_to_output;
#    my $element;
#    my $max_cons = $Output[0]{GRADE};
#    my $ConsColorUnity; #unity of conservation to be colored
#    foreach $element (@Output){
#        if ($$element{GRADE} < $max_cons) {$max_cons = $$element{GRADE};}
#    }
#    $ConsColorUnity = $max_cons / 4.5 * -1; 
#    if ($max_cons !~ /^\-/){$ConsColorUnity = $max_cons;}
#    return ($max_cons, $ConsColorUnity);
#}
##------------------------------
##calculates the grades for each color
#sub calc_layers{
#    my ($MaxCons,$ConsColorUnity, $ref_to_ColorLayers) = @_;
#    my $i;
#    my $NoLayers = 9;
#    for ($i = 0; $i <= $NoLayers; $i++) {
#        $ref_to_ColorLayers->[$i] = $MaxCons + ($i * $ConsColorUnity);
#    }
#}
##------------------------------
##gives the color to the interval
#sub assign_colors_interval{
#    my ($ref_to_output, $ref_to_layers) = @_;
#    my @Output = @$ref_to_output;
#    my @ColorLayers = $ref_to_layers;
#    my $i;
#    my $element;
#    my $Count = 0;
#
#    foreach $element (@Output) {
#        
#        for ($i = 0; $i <= $#ColorLayers; $i++){
#            if( ($i==$#ColorLayers) and !exists $$element{INTERVALLOWCOLOR}){
#                $$element{INTERVALLOWCOLOR} = 8;
#            }
#            elsif ($$element{INTERVALLOW} >= $ColorLayers[$i] and $$element{INTERVALLOW} < $ColorLayers[$i + 1]) {
#                $$element{INTERVALLOWCOLOR} =$i;
#            } 
#            elsif ( ($$element{INTERVALLOW} < $ColorLayers[$i]) and !exists $$element{INTERVALLOWCOLOR}){
#                $$element{INTERVALLOWCOLOR} = 0;
#            } 
#            if (($i == $#ColorLayers)  and !exists $$element{INTERVALHIGHCOLOR}){
#                $$element{INTERVALHIGHCOLOR} = 8;
#            }
#            elsif ($$element{INTERVALHIGH} >= $ColorLayers[$i] and $$element{INTERVALHIGH} < $ColorLayers[$i + 1]) {
#                $$element{INTERVALHIGHCOLOR} =$i;
#            } 
#            elsif ( ($$element{INTERVALHIGH} < $ColorLayers[$i]) and !exists $$element{INTERVALHIGHCOLOR}){
#                $$element{INTERVALHIGHCOLOR} = 0;
#            } 
#        } # END FOR
#    } # END FOREACH
#} 
##------------------------------
##give the color for each position based on the grades
#sub assign_colors{    
#    my ($ref_to_output, $ref_to_layers) = @_;
#    my @Output = @$ref_to_output;
#    my @ColorLayers = $ref_to_layers;
#    my $bayesInterval = CONSURF_CONSTANTS::BAYES_INTERVAL;
#    
#    # match the colors to the grades
#    foreach my $element (@Output){ 
#        for (my $i = 0; $i <= $#ColorLayers; $i++) {
#            if ($i == $#ColorLayers) {
#                $$element{COLOR} = $ColorScale{$i-1};
#            }
#            elsif ($$element{GRADE} >= $ColorLayers[$i] && $$element{GRADE} < $ColorLayers[$i + 1]) {
#                $$element{COLOR} = $ColorScale{$i};         
#                last;
#            }            
#        }
#        if ((($$element{INTERVALHIGHCOLOR}-$$element{INTERVALLOWCOLOR})>$bayesInterval ) or ($$element{MSA_NUM} <= 5)){
#            $$element{ISD}=1;
#        }
#        else{$$element{ISD}=0;}
#    }  
#}
##------------------------------
#sub read_residue_variety
#{
#	my ($msa_filename, $ref_to_residue_variety_hash) = @_;
#	my $line;
#	my $position = 1;
#	my $data_length;
#	my @elements_to_remove = ();
#
#	# open msa file
#	unless (open MSA, $msa_filename){
#        return ("new_rate4site_routines::read_residue_variety : can't open '$msa_filename' for reading $!");
#
#	# read header lines
#	$line = <MSA>;
#	$line = <MSA>;
#	while ($line = <MSA>)
#	{
#	   if ($line =~ /^(.+) +([KRHDEYWSTFLIMCNQAVPG-]+)$/)
#	    {
#	    	my $seq_name = $1;
#	        my $seq_data = $2;
#	        $seq_name =~ s/\s+$//; # remove spaces
#	        $data_length = length($seq_data);
#	        for (my $position_offset = 0;$position_offset < length($seq_data);$position_offset++)
#	        {
#	            # get residue at this position
#	            my $current_aa = substr($seq_data,$position_offset,1);
#	            
#	            # if this is the query seq and no aa is found we need to save
#	            # this position in order to erase it in the end
#	            if ($current_aa =~ /^-$/ and $seq_name eq $pdb_id.$chain_id)
#	            {
#	            	push @elements_to_remove,$position+$position_offset;
#	            }
#	            # skip if no residue found
#	            if ($current_aa =~ /^-$/)
#	            {
#	                next;
#	            }
#
#	            # save residue value    
#	            if (!exists $ref_to_residue_variety_hash->{$position+$position_offset})
#	            {
#	                $ref_to_residue_variety_hash->{$position+$position_offset} = $current_aa;
#	            }
#	            elsif (index($ref_to_residue_variety_hash->{$position+$position_offset},$current_aa) == -1)
#	            {
#	                $ref_to_residue_variety_hash->{$position+$position_offset} .= ",$current_aa";
#	            }
#	        }
#	    }
#	    else
#	    {
#	    	#inc position and read next (blank) line
#	       	$position += $data_length;
#	       	$line = <MSA>;
#	    }
#	}
#	close MSA;
#	
#	# remove all positions where the query seq was skipped
#	my $index = pop(@elements_to_remove);
#	while (defined($index))
#	{
#		residues_shift_left($index, $ref_to_residue_variety_hash);
#		$index = pop(@elements_to_remove);
#	}
#    return "OK";
#}
##------------------------------
#sub residues_shift_left
#{
#    my ($start_position,$ref_to_residue_variety_hash) = @_;
#    my $index = $start_position;
#    
#    while (exists($ref_to_residue_variety_hash->{$index+1}))
#    {
#        $ref_to_residue_variety_hash->{$index} = $ref_to_residue_variety_hash->{$index+1};
#        $index++;
#    }
#    delete($ref_to_residue_variety_hash->{$index});
#    
#}
##------------------------------
##printing the reults into a file
#sub print_gradesPE_output{
#    my ($output_file, $ref_to_output_array) = @_;
#    my@Output = @$ref_to_output_array;
#    
#    # open file
#    unlness (open PE, ">$output_file"){
#        return ("new_rate4site_routines::print_gradesPE_output : can't open '$output_file' for writing $!");}
#    print PE "\t Amino Acid Conservation Scores\n";
#    print PE "\t===============================\n\n";
#    print PE "- POS: The position of the AA in the SEQRES derived sequence.\n";
#    print PE "- SEQ: The SEQRES derived sequence in one letter code.\n";
#    print PE "- 3LATOM: The ATOM derived sequence in three letter code, including the AA's positions as they appear in the PDB file and the chain identifier.\n";
#    print PE "- SCORE: The normalized conservation scores.\n";
#    print PE "- COLOR: The color scale representing the conservation scores (9 - conserved, 1 - variable).\n";
#    print PE "- CONFIDENCE INTERVAL: When using the bayesian method for calculating rates, a confidence interval is assigned to each of the inferred evolutionary conservation scores.\n"; 
#    print PE "- CONFIDENCE INTERVAL COLORS: When using the bayesian method for calculating rates. The color scale representing the lower and upper bounds of the confidence interval.\n"; 
#    print PE "- MSA DATA: The number of aligned sequences having an amino acid (non-gapped) from the overall number of sequences at each position.\n";
#    print PE "- RESIDUE VARIETY: The residues variety at each position of the multiple sequence alignment.\n\n";
#    print PE " POS\t SEQ\t    3LATOM\tSCORE\t\tCOLOR\tCONFIDENCE INTERVAL\tCONFIDENCE INTERVAL COLORS\tMSA DATA\tRESIDUE VARIETY\n";
#    print PE "    \t    \t        \t(normalized)\t        \t               \n";
#    foreach my $elem (@Output){
#        printf (PE "%4d", "$$elem{POS}");
#            printf (PE "\t%4s", "$$elem{SEQ}");
#        printf (PE "\t%10s", "atom_res");
#        printf (PE "\t%6.3f", "$$elem{GRADE}");
#        if($$elem{ISD}==1){
#            printf (PE "\t\t%3d", "$$elem{COLOR}");
#            printf (PE "%1s", "*");
#        }
#        else{printf (PE "\t\t%3d", "$$elem{COLOR}");}
#        printf (PE "\t%6.3f", "$$elem{INTERVALLOW}");
#        printf (PE "%1s", ",");
#        printf (PE "%6.3f", "$$elem{INTERVALHIGH}");
#        printf (PE "\t\t\t%5d", "$ColorScale{$$elem{INTERVALLOWCOLOR}}");
#        printf (PE "%1s", ",");
#        printf (PE "%1d\t\t", "$ColorScale{$$elem{INTERVALHIGHCOLOR}}");
#        printf (PE "\t%8s", "$$elem{MSA_NUM}\/$$elem{MSA_DENUM}");
#        printf (PE "\t%-18s\n", "X");
#    }
#    print PE "\n\n*Below the confidence cut-off - The calculations for this site were performed on less than 6 non-gaped homologue sequences,\nor the confidence interval for the estimated score is equal to- or larger than- 4 color grades.\n";
#    close(PE);
#}

#------------------------------
#------------------------------
1;