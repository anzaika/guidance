#!/usr/bin/perl -w

package design_pipe;

use strict;


# we should create string that will be of format:
# (1, 1, "Tyr32:H,Tyr27:H,,Ser25:H,Ala24:H,Ser77:H", "1 Tyr-Tyr32:H,2 Tyr-Tyr27:H,3 Arg,4 Asn-Ser25:H,5 Ala-Ala24:H,6 Ser-Ser77:H")
# item1: from 1 to no. of peptides found in the file peptides.txt
# item2: from 1 to no. of alignments found in the file sig_path_align that match that peptide (from item1)
# item3, item4: designed here, according to text inside sig_path_align that match that peptide
# 
# the file pipe_txt contains the info to produce all the relevant strings to the pipe regarding clusters and their peptides.
# read this file and extract the relevant info from it using the routine create_path_alignment_string that parses the peptides alignments files

#output: a list of 2 items. ("OK", "") : if no error was found
#  ("error", <description>) : if an error was found
sub create_cluster_lines_for_pipe
{
    my $pipe_txt = shift; # The text file which helps numbering of clusters
    my $SigPathAln_dir = shift;
    my $pipe_pdb = shift; # This is the pipe pdb file, which is the output file for this routine
    
    my $no_of_clusters;
    my $no_of_clusters_line;
    my $cluster_id;
    my $no_of_peptides;
    my $no_of_peptides_line;
    my $peptide_counter=0;  # Assuming the number to be reported as peptide is its serial number AND NOT its ID.
            # if we wish to change it to be its ID, Than we can get this info from the third elsif
    my @pepsurf_peptides;  # To hold the info of format : ! pepsurf_peptide[1] = "YYRNAS";
    my @alignments; #To hold the strings returned from create_path_alignment_string
    my @record_alignments; # To hold the info of format: ! pepsurf_record_alignment(1, 1, "Tyr32:H,Tyr27:H,,Ser25:H,Ala24:H,Ser77:H", "1 Tyr-Tyr32:H,2 Tyr-Tyr27:H,3 Arg,4 Asn-Ser25:H,5 Ala-Ala24:H,6 Ser-Ser77:H")
    
    
    unless (open PIPE_TXT, "<$pipe_txt"){
        return ("error", "design_pipe::create_cluster_lines_for_pipe  : can't open $pipe_txt for reading.\n");
    }
    while (<PIPE_TXT>)
    {
        chomp;
        if(/^! pepsurf_cluster_count = (\d+)/)
        {
            $no_of_clusters = $1;
            $no_of_clusters_line = $_;
        }
        elsif(/^CLUSTER (\d+)/)
        {
            $cluster_id = $1;
        }
        elsif(/^(\d+)_(\w+)/) # $1: peptide_id  $2: peptide sequence        
        {
            $peptide_counter++;
            $pepsurf_peptides[$peptide_counter-1] = "! pepsurf_peptide[".$1."] = \"$2\";\n";
            # open alignment file according to peptide id
            @alignments = &create_path_alignment_string($SigPathAln_dir.($1-1)."_significantPaths.txt", $2);
            if ($alignments[0] eq "Error")
            {
                return ("error", $alignments[1]);
            }
            else
            {
                $record_alignments[$peptide_counter-1] = "! pepsurf_record_alignment($cluster_id, $1, \n".&print_string_max80($alignments[0], "!,")."\",\n".&print_string_max80($alignments[1], "!,")."\");\n";
            }
        }
        elsif(/^! pepsurf_peptide_count = (\d+)/)
        {
            $no_of_peptides = $1;
            $no_of_peptides_line = $_;
        }
    }#finish while, closing file:
    close PIPE_TXT;
    
    unless (open PIPE_PDB, ">>$pipe_pdb")
    {
        return ("error", "design_pipe::create_cluster_lines_for_pipe  : can't open $pipe_pdb for reading.\n");
    }
    print PIPE_PDB $no_of_clusters_line."\n";
    foreach (@pepsurf_peptides)
    {
        print PIPE_PDB $_;
    }
    print PIPE_PDB "!\n";
    print PIPE_PDB "!! PEPTIDE ALIGNMENTS FOR CLUSTER 1:\n";
    print PIPE_PDB "!! FIRST TWO PARAMETERS: CLUSTER [1] PEPTIDE [1].\n";
    print PIPE_PDB "!! THIRD PARAMETER: ALIGNED RESIDUES \"PATH\" (POSSIBLY SOME NOT IN THE CLUSTER)\n";
    print PIPE_PDB "!! FOURTH PARAMETER: ALIGNMENT \n!!\n";
    foreach (@record_alignments)
    {
        print PIPE_PDB $_;
    }
    close PIPE_PDB;
    return ("OK", "");
}

#1.  sub routine:
#    convert amino acids codes: 1 letter to 3 letters in format Capital-lower-lower
#2.  convert aligmnment to 2 lines:
#    a. peptide (3 letters next to numbers)
#    b. alignment (3 letters against 3 letters)


# open sigPathsAln directory.
# open first file (begins with 0
# while !EOF
# read line begins in "path:"
# convert it to 3 letters   --> $converted_aligned_path
# read after Alignemnt:  --> $alignment_final_string
#/data/bental-data/PepSurf/pepSurf-Results/1160243370/sigPathsAln



# input: 1. sigPathsAln file 2. the peptide which was read from the "cluster" file
# When the routine reads a peptide, it first checks that the peptide sequence matched the sequence which
# was read in the cluster file. If there is no match - aborting.
#  *** The sub routine assumes the file is not empty!! ***
# output:list of 2 strings
# 1. the conversion of a path string :    Y32H Y27H S25H A24H S77H
#   to array with 3 letter (xxx) instead of 1, followed by residue no. and chain. format: XXXNO:Chain
# 2. the conversion of the alignment of the peptide with the cluster in 3 letters, according to PiPe requested format
#

sub create_path_alignment_string
{
    my $sigPathsAln = shift;
    my $peptide_sequence_from_cluster = shift;

    my %amino_converter = ('a' => 'Ala', 'A' => 'Ala', 'r' => 'Arg', 'R' => 'Arg', 'n' => 'Asn', 'N' => 'Asn', 'd' => 'Asp', 'D' => 'Asp',
                             'b' => 'Asx', 'B' => 'Asx', 'c' => 'Cys', 'C' => 'Cys', 'e' => 'Glu', 'E' => 'Glu', 'q' => 'Gln', 'Q' => 'Gln',
                             'z' => 'Glx', 'Z' => 'Glx', 'g' => 'Gly', 'G' => 'Gly', 'h' => 'His', 'H' => 'His', 'i' => 'Ile', 'I' => 'Ile',
                             'l' => 'Leu', 'L' => 'Leu', 'k' => 'Lys', 'K' => 'Lys', 'm' => 'Met', 'M' => 'Met', 'f' => 'Phe', 'F' => 'Phe',
                             'p' => 'Pro', 'P' => 'Pro', 's' => 'Ser', 'S' => 'Ser', 't' => 'Thr', 'T' => 'Thr', 'w' => 'Trp', 'W' => 'Trp',
                             'y' => 'Tyr', 'Y' => 'Tyr', 'v' => 'Val', 'V' => 'Val');    
    # vars which read from the input file:
    my $peptide;  # to hold the peptide sequence, 1 letter code
    my $aligned_path; # to hold the path, 1 letter code
    my $aligned_peptide; # to hold the alignment of the path to the peptide (including "-" signs)
    #my $alignments_counter = 0; # to hold the number of aligned_peptide for that peptide
    
    # vars which cerated by the routine
    my $line;
    my @aligned_path_array;
    my $converted_aligned_path= ""; # item3 in Eric's request
    my @peptide_array;
    my $alignment_final_string = ""; #item4 in Eric's request
    my $aligned_path_counter = 0;
    
    unless (open SigPathAln, "<$sigPathsAln"){
        return ("Error", "design_pipe::create_path_alignment_string  : can't open file $sigPathsAln for reading.\n");
    }
    
    # According to spesification we only read the first alignment from the file. Therefore the 
    while (<SigPathAln>)
    {  
        chomp;
        if (/^peptide:\s+(\w+)/)
        {
            $peptide = $1;
            if($peptide_sequence_from_cluster ne $1){
                return ("Error", "design_pipe::create_path_alignment_string  : The sequence $peptide_sequence_from_cluster which was sent from the cluster, does not match the sequence $1 which is found in the file $sigPathsAln. Aborting\n.");
            }
        }
        elsif (/^path:\s+(.+)/)
        {
            $aligned_path = $1;
        }
        elsif (/^Alignment:/)
        {
            $peptide = <SigPathAln>;  #after we varified that the peptide is equal to the one found in the cluster, we change it to the piece of peptide found in the alignment. (since the aligned one might be shorter that the original peptide)
            chomp($peptide);
            $line = <SigPathAln>;
            chomp($line);
            $aligned_peptide = $line;
        }
        elsif(/^_+/)
        {
            last;
        }
    }
    close SigPathAln;
    #my $aligned_path = "Y32H Y27H S25H A24H S77H"; #this info is read from the file
    # put the $path in an array, seperate by space delimiter
    
    @aligned_path_array = split(/\s+/, $aligned_path);
    foreach (@aligned_path_array)
    {
    #    if ($_ =~ m/(\w)(\d+)(\w)/)
    #    {
    #        $_ = $amino_converter{$1}.$2.":".$3;        
    #    }
    #}
    
        if ($_ =~ m/(\w)(\d+)(.?)/){
            my $maybe_chain = $3;
            $_ = $amino_converter{$1}.$2;
            if($maybe_chain =~ /\w/) {# in case there is no chain, we don't write anything after the ':'
                $_.=":".$maybe_chain;
            }
            #else{
            #    $_.=$maybe_chain.":";
            #}
        }
    }    
    
     
    #my $peptide = "YYRNAS";  #this info is read from the file
    # put the $peptide in an array, translated
    
    for (my $key=0; $key<length($peptide); $key++)
    {
        $peptide_array[$key] = $amino_converter{substr($peptide, $key, 1)};
        #print $peptide_array[$key]."\n";
    }
    
    #read the aligned line parralel to a counter.
    
    #my $aligned_peptide = "YY-SAS"; #this info is read from the file

    for (my $peptide_counter=0; $peptide_counter<length($peptide); $peptide_counter++)
    {
        $alignment_final_string.= ($peptide_counter+1)." ".$peptide_array[$peptide_counter];        
        if (substr($aligned_peptide, $peptide_counter, 1) eq "-")
        {
            if($peptide_counter==0) # if the alignment begins with "-" we add an extra "," to the beginning
            {
                $converted_aligned_path.=",|,"; 
            }
            elsif($peptide_counter+1 == length($peptide)) # if the alignment ends with "-" we add an extra "," to the end
            {
                $converted_aligned_path.="|,|,";
            }
            else
            {                
                $converted_aligned_path.="|,"; #this is an identifier, for the next string th know it had a gap.
            }
            $alignment_final_string.=","; #this is an identifier, for the next string th know it had a gap.
        }
        else
        {
            $alignment_final_string.="-".$aligned_path_array[$aligned_path_counter].",";
            $converted_aligned_path.=$aligned_path_array[$aligned_path_counter].",";
            $aligned_path_counter++;
        }
    }
    $aligned_path_counter = 0;
    chop($alignment_final_string); #removing the last ","
    chop($converted_aligned_path);
    
    return ($converted_aligned_path, $alignment_final_string);
}

#---------------------------------------------------------
#
    # breaking string "part 1 part 2 part  3" in one of 2 ways, depending on the $begin_sign:
    # ! "part 1" +
    # ! "part 2" +
    # ! "part 3"
    # or
    # ! select part1
    # ! select selected or part 2
    # ! select selected or part 3
    
    # in case ! the delimiter is ,. therfore when ,, are in the string, they are converted (in advance) to ,| 

sub print_string_max80($$)
{
    my $input_string = shift; #"Ala24:H, Ser25:H, Tyr27:H, Thr30:H, Tyr32:H, Thr74:H, Ser75:H, Lys76:H, Ser77:H";
    my $begin_sign = shift; #"select";
    
    my $final_string;
    my @string;
    my $line_length=0; 
    
    if ($begin_sign eq "!,")
    {
        @string = split(/,/ , $input_string);
        $final_string = "! \"";
        $line_length = length $final_string;
        foreach (@string)
        {
            if ($_ =~ m/\|/)
            {
                $_=""
            } 
            $line_length += ((length $_)+1); # add 1, since it will be the length of the added ","
            if ($line_length > 73)
            {
                $final_string.= "\" +\n! \"$_,";
                $line_length = (length $_);
            }
            else
            {
                $final_string.= $_.",";
                #if($final_string =~ m/,,,$/){
                #    chop($final_string);
                #}
            }        
        }
        chop($final_string);
    }
    elsif ($begin_sign eq "!")
    {
        @string = split(/ /, $input_string);
        $final_string = "! \"";
        $line_length = length $final_string;
        foreach (@string)
        {
            unless ((length $_) ==0)
            {
                $line_length += ((length $_)+1); # add 1, since it will be the length of the added ","
                if ($line_length > 73)
                {
                    $final_string.= "\" +\n! \"$_ ";
                    $line_length = (length $_);
                }
                else
                {
                    $final_string.= $_." ";
                }
            }
        }
        chop($final_string);
    }    
    elsif($begin_sign eq "select")
    {
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
    }
    
    return "$final_string";    
}
1;