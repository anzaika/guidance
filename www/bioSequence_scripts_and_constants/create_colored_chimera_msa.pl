#! /usr/bin/perl -w

use lib "/bioseq/bioSequence_scripts_and_constants";

use rasmol_gradesPE_and_pipe;
use strict;


my $dir = "/bioseq/data/results/ConSurf/1217408435/";
my $query_name = "1php_A";
my $msa_file = $dir."pdb1php.aln";
my $rate4site_filename = $dir."bayesR4s.res";
my $chimera_script_name = $dir."pdb1php_aln.scf";
my $headers = $dir. "pdb1php_aln.hdr";
my $chimera_script_ins = $dir."pdb1php_aln_isd.scf";
my $headers_ins = $dir. "pdb1php_aln_isd.hdr";

my @r4s_position_grade = ();
my %aln_position_grade = ();
my %aln_position_grade_isd = ();

my %consurf_rasmol_colors=("0"=>" 255 255 150","1"=>" 16 200 209","2"=>" 140 255 255","3"=>" 215 255 255","4"=>" 234 255 255","5"=>" 255 255 255","6"=>" 252 237 244","7"=>" 250 201 222","8"=>" 240 125 171","9"=>" 160 37 96");


my @ans = rasmol_gradesPE_and_pipe::assign_colors_according_to_r4s_layers($rate4site_filename, \@r4s_position_grade);
if ($ans[0] ne "OK") {die $ans[0];}

my $ref_r4s_position_grade = \@r4s_position_grade;


#open MSA file, read only the query sequence lines. assign each position with increasing number. if it is not a '-' : read the grade from the array @r4s_position_grade and put it in the hash %aln_position_grade.
#if there is insufficient data: save it to a saparate hash

my $position_in_msa = 0;
my $position_in_r4s = 0;
open MSA, $msa_file or die "cannot open $msa_file $!";
while(<MSA>){
    chomp;
    if (/^(.+)\s+([A-Za-z-]+)$/){
        my $seq_name = $1;
        my $seq = $2;
        $seq_name =~ s/\s+$//;
        next if ($seq_name ne $query_name);        
        # we are in the query name
        $seq =~ s/\s+$//; # remove spaces
        my $seq_length = length($seq);
        my @sequence_block = split //, $seq;
        foreach(@sequence_block){
            unless(/-/){
                if ($ref_r4s_position_grade->[$position_in_r4s]{ISD} == 1){
                    $aln_position_grade_isd{$position_in_msa} =0;
                }
                else{
                    $aln_position_grade_isd{$position_in_msa} = $ref_r4s_position_grade->[$position_in_r4s]{COLOR};
                }
                $aln_position_grade{$position_in_msa} = $ref_r4s_position_grade->[$position_in_r4s]{COLOR};
                $position_in_r4s++;
            }
            $position_in_msa++;
        }
    }
}
close MSA;

print_chimera_scripts($chimera_script_name, $headers, \%aln_position_grade);
print_chimera_scripts($chimera_script_ins,$headers_ins, \%aln_position_grade_isd);


sub print_chimera_scripts{
    my $chimera_script_name = shift;
    my $chimera_headers = shift;
    my $ref_aln_position_grade = shift;
    
    open SCRIPT, ">".$chimera_script_name;
    open HEADER, ">".$chimera_headers;
    print HEADER "name: Conservation scores\nstyle: character\n";
    my $j=0;
    for(my $i=0; $i<=$position_in_msa; $i++){
        $j=$i+1;
        if (exists $aln_position_grade{$i}){        
            print SCRIPT "$i $i 0 0 $consurf_rasmol_colors{$ref_aln_position_grade->{$i}}\n";
            print HEADER "\t$j\t$ref_aln_position_grade->{$i}\tblack\n";
        }
    }
    close SCRIPT;
    
    print HEADER "#  ConSurf data shown as histogram\nname: ConSurf histogram\nstyle: numeric\n";
    $j=0;
    my $histo;
    for(my $i=0; $i<=$position_in_msa; $i++){
        $j=$i+1;
        if (exists $aln_position_grade{$i}){
            if ($ref_aln_position_grade->{$i} == 9){$histo = "1.0";}
            else{$histo = ".".($ref_aln_position_grade->{$i}+1);}
            
            print HEADER "\t$j\t$histo\tblack\n";
        }
    }
    close HEADER;
}
