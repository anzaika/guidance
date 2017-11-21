#!/usr/bin/perl -w

# 1. choose homolougs from blast - according to some parameters
# in case the blast input is psi-blast with few rounds: runs only on the last round hits
# 2. send these homolougs to cd-hit
# 3. filter cd-hit output and choose final homolougs to send to Muscle
# 4. runs muscle on chosen homolougs, create MSA in clustAlW format


use lib "/db1/System/var/www/html/bio/test/external_scripts";
use parseFiles;
use strict;

my $input_details = $ARGV[0];

#-----------------------
# reading the input file
#-----------------------

my ($query, $blast_output, $working_dir, $redundancyRate, $frag_overlap, $min_length_percent, $MaxNumHomol);
my $Max_num_given = 0;

if ($input_details eq "") {die "The script should be given the in_file.txt path as an argument!\n";}
open IN_FILE, $input_details or die "cannot open file $input_details for reading $!";
while(<IN_FILE>){
    chomp;
    if (/^QUERY_FASTA =\s?(.+)/) {$query = $1;}
    elsif (/^BLAST_IN =\s?(.+)/) {$blast_output = $1;}
    elsif (/^OUTPUT_DIR =\s?(.+)/) {$working_dir = $1;}
    elsif (/^REDUNDANCY_RATE =\s?(.+)/) {$redundancyRate = $1;}
    elsif (/^FRAGMENT_OVERLAP =\s?(.+)/) {$frag_overlap = $1;}
    elsif (/^MINIMUM_LENGTH =\s?(.+)/) {$min_length_percent = $1;}
    elsif (/^MAX_NUM_HOMOLOUGS =\s?(.+)/) {$MaxNumHomol = $1;}
}
close IN_FILE;

#--------------------------------------
# checking that the input file is legal
#--------------------------------------
# check that all the input was given
if (!($query)){
    die "The full path to the query file was not found in the input file in_file.txt.\nPlease write it in the line that begins with: QUERY_FASTA and re-run the script\n";}
elsif (!($blast_output)){
    die "The full path to the blast output file was not found in the input file in_file.txt.\nPlease write it in the line that begins with: BLAST_IN and re-run the script\n";}
elsif(!($working_dir)){
    die "The full path to the output-files directory was not found in the input file in_file.txt.\nPlease write it in the line that begins with: OUTPUT_DIR and re-run the script\n";}
elsif(!($redundancyRate)){
    die "The redundancy rate value was not found in the input file in_file.txt.\nPlease write it in the line that begins with: REDUNDANCY_RATE and re-run the script\n";}
elsif(!($frag_overlap)){
    die "The maximum value for overlapping was not found in the input file in_file.txt.\nPlease write it in the line that begins with: FRAGMENT_OVERLAP and re-run the script\n";}
elsif(!($min_length_percent)){
    die "The minimum length for a homoloug sequence was not found in the input file in_file.txt.\nPlease write it in the line that begins with: MINIMUM_LENGTH and re-run the script\n";}
elsif(!($MaxNumHomol)){
    print "The maximum number of homolougs to collect was not found in the input file in_file.txt.\nPlease note that the constructed MSA may be very large\n";}

#check that all the input is correct

$redundancyRate =~ s/\s//g;
$frag_overlap =~ s/\s//g;
$min_length_percent =~ s/\s//g;
if ($MaxNumHomol){$MaxNumHomol=~ s/\s//g;}
if(!($working_dir=~/\/$/)){$working_dir.="/";}

if($redundancyRate<0 || $redundancyRate>100 || (!($redundancyRate=~m/^\d+$/))){
    die "The redundancy rate value given is: \'$redundancyRate\'\nPlease correct this value in the input file in_file.txt, as it should be an integer between 0-100\n";}
elsif($frag_overlap<0 || $frag_overlap>=1){
    die "The maximum value for overlapping given is: \'$frag_overlap\'\nPlease correct this value in the input file in_file.txt, as it should be a decimal fraction between 0.00-0.99\n";}
elsif($min_length_percent<0 || $min_length_percent>1){
    die "The minimum length for a homoloug sequence given is: \'$min_length_percent\'\nPlease correct this value in the input file in_file.txt, as it should be a decimal fraction between 0.00-1.00\n";}

elsif(!($MaxNumHomol)){
    print "The maximum number of homolougs to collect was not found in the input file in_file.txt.\nPlease note that the constructed MSA may be very large\n";}
elsif((!($MaxNumHomol =~ m/^\d+$/)) || $MaxNumHomol<6){
    die "The maximum number of homolougs given is: \'$MaxNumHomol\'\nPlease correct this value in the input file in_file.txt, as it should be an integer greater than 6\n";}
elsif ($MaxNumHomol ne ""){
    $Max_num_given = 1;}

#-----------------------------------
# defining variables for the program
#-----------------------------------

# OUTPUT FILES
my $blast_last_round = "last_round.blast";  # a copy of the last round of the given blast input
my $fasta_output = "legal_sequences.fasta";  # the final homolouges file, to write to
my $rejected_seqs = "rejected_sequences.txt"; # file that will have the rejected sequences and the reason for the rejection
my $log = "run.log"; # log file for the run of this script
my $cd_hit_output = "cd_hit_chosen_sequences.fasta"; # the chosen sequences according to cd-hit
my $final_homolougs_file = "final_sequences.fasta"; # the final sequences to send to the MSA
my $muscle_output = "MSA.aln";
# EXTERNAL PROGRAMS PATHS
my $cd_hit_dir = "/db1/Local/src/cd-hit_redundency/";         # cd-hit program
my $muscle =  "/usr/local/bin/muscle";

my $min_num_of_homolougs = 5; # the minimum number of homolougs we must collect
my %blast_hash;
my %cd_hit_hash;
my $blast_rounds = "none";
my $round_flag = 0;
my $error;

open LOG, ">".$working_dir.$log or die "cannot open the LOG file \'".$working_dir.$log."\' for writing $!";
print LOG "*** input parameters read from in_file:\n";
print LOG "query : \'$query\'\nblast_output : \'$blast_output\'\nworking_dir : \'$working_dir\'\n";
print LOG "redundancyRate : \'$redundancyRate\'\nfrag_overlap : \'$frag_overlap\'\n";
print LOG "min_length_percent : \'$min_length_percent\'\n";
if ($Max_num_given==1) {print LOG "MaxNumHomol : \'$MaxNumHomol\'\n";}
else {print LOG "MaxNumHomol was not given\n\n";}

#------------------------------------------------------------
# check if 1. hits were found in the blast
# 2. the blast file is composed of few rounds, if so -
# copies only the last round to a new file
#------------------------------------------------------------
unless (open BLAST, $blast_output){
    $error = "cannot open the file \'$blast_output\' for reading $!";
    &print_to_log_and_exit($error);
}
while(<BLAST>){
    if(/No hits found/){
        &print_to_log_and_exit("No hits found is the Blast output file. Cannot continue calculation");
    }
    elsif(/^Searching/){
        unless(/done/){
            &print_to_log_and_exit("No hits found is the Blast output file. Cannot continue calculation");
        }
    }
    elsif(/Results from round (\d+)/){
        $blast_rounds = $1;
    }
}
close BLAST;


# in case there are few rounds: create a new file only with the last round
if ($blast_rounds ne "none"){
    print LOG "going to read only blast results from round $blast_rounds\n";
    unless(open BLAST, $blast_output){ 
        $error = "cannot open the file \'$blast_output\' for reading $!";
        &print_to_log_and_exit($error);
    }
    unless (open BLAST_OUT, ">".$working_dir.$blast_last_round){
        $error = "cannot open the file \'".$working_dir.$blast_last_round."\' for writing $!";
        &print_to_log_and_exit($error);
    }
    while(<BLAST>){
        if(/Results from round (\d+)/){
            if($blast_rounds == $1){
                $round_flag=0;
            }
            else
                {$round_flag=1;}            
        }
        if ($round_flag==0) {print BLAST_OUT $_;}
    }
    close BLAST_OUT;
    close BLAST;
    if ((!(-e $working_dir.$blast_last_round)) || (-z $working_dir.$blast_last_round)){
        $error = "Could not create a new blast file that will include only the last round hits. abort. $!";
        &print_to_log_and_exit($error);    
    }
    else {
        $blast_output = $working_dir.$blast_last_round;}
}


#--------------------------------------------------------------
# choosing homolougs, create fasta file for all legal homolougs
#--------------------------------------------------------------

print LOG "running parseFiles::choose_homologoues_from_blast\n";
my @ans = parseFiles::choose_homologoues_from_blast($working_dir, $query, $redundancyRate, $frag_overlap, $min_length_percent, $min_num_of_homolougs, $blast_output, $fasta_output, $rejected_seqs, \%blast_hash);
if ($ans[0] eq "sys" || $ans[0] eq "user"){
    $error = $ans[1];
    &print_to_log_and_exit($error);
}

#-----------------------------------------------------------
# screen homolougs by redundancy rate, according to clusters
#-----------------------------------------------------------

print LOG "running parseFiles::create_cd_hit_output\n";
$redundancyRate = ($redundancyRate/100);
@ans = parseFiles::create_cd_hit_output($working_dir, $fasta_output, $cd_hit_output, $redundancyRate, $cd_hit_dir, \%cd_hit_hash);
if ($ans[0] eq "sys" || $ans[0] eq "user"){
    $error = $ans[1];
    &print_to_log_and_exit($error);
}

#--------------------------------------------------
# Check that the number of homolougs found is legal
#--------------------------------------------------
my $final_num_homolougs = (keys %cd_hit_hash);
my $final_num_blast_homolougs = (keys %blast_hash);
my $message;
if ($final_num_homolougs == 1){
    $error = "There is only 1 unique hit. The minimal number of sequences required for building a MSA is 5.\n";
    &print_to_log_and_exit($error);
}
elsif ($final_num_homolougs <= $min_num_of_homolougs){
    $error = "There are only $final_num_homolougs unique hits. The minimal number of sequences required for the calculation in order to build a MSA is $min_num_of_homolougs.\n";
    &print_to_log_and_exit($error);
}
if ($Max_num_given == 1){
    if (($final_num_blast_homolougs > $final_num_homolougs) && ($final_num_homolougs> $MaxNumHomol)){
        $message = "There are $final_num_blast_homolougs PSI-BLAST unique sequences, of which $final_num_homolougs are considered \"clean\".\nThe calculation is performed on the $MaxNumHomol clean sequences with the lowest E-value.\n";
        print LOG $message;
        print $message;
    }
    elsif ($final_num_homolougs> $MaxNumHomol){
        $message = "There are $final_num_homolougs unique PSI-BLAST hits.\nThe calculation is performed on the $MaxNumHomol clean sequences with the lowest E-value.\n";
        print LOG $message;
        print $message;
    }
}

#-------------------------------------------------
# extracting query info and choose final homolougs
#-------------------------------------------------

my ($s_name, $s_aa_sq, $query_AAseq, $query_name, $counter);
unless (open QUERY, $query){
    $error = "could not open the file \'$query\' for reading $!";
    &print_to_log_and_exit($error);
}
while (<QUERY>){
    chomp;
    if (/^>(.+)/){
        $query_name = $1;
    }
    else{
        $query_AAseq.= $_;
    }
}
close QUERY;

print LOG "going to choose the best $MaxNumHomol hits for the MSA\n";
$counter = 1;
open FINAL, ">".$working_dir.$final_homolougs_file;
# print query details first
print FINAL ">$query_name\n"."$query_AAseq\n";
foreach $s_name (sort { $blast_hash{$a} <=> $blast_hash{$b} } keys %blast_hash ){
    if ($Max_num_given==1 && $counter>$MaxNumHomol) {last;}    
    if (defined($cd_hit_hash{$s_name})){
        $s_aa_sq = $cd_hit_hash{$s_name};
        $s_name =~ m/(\S+)\|.+(_\d+_\d+)/;
        print FINAL ">$1"."$2\n".$s_aa_sq."\n";
        $counter++;
    }    
}
close FINAL;

#-----------
# run muscle
#-----------

print LOG "going to run MUSCLE\n";
my $cmd = "$muscle -in $working_dir"."$final_homolougs_file -out $working_dir".$muscle_output."  -clwstrict  -quiet";
`$cmd`;

if((!(-e $working_dir.$muscle_output)) || (-z $working_dir.$muscle_output)){
    $error = "For some reason MUSCLE did not construct the MSA for you.\nYou can try to run it yourself. this is the command line:\n".$cmd;
    &print_to_log_and_exit($error);
}
close LOG;


#----------------------------------
# sub routine print_to_log_and_exit
#----------------------------------
sub print_to_log_and_exit(){
    my $error = shift;
    
    print LOG $error;
    close LOG;
    die($error);
}