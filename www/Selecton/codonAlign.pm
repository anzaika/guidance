#!/usr/local/bin/perl -w

package codonAlign;

# checking if sequence input is correct for aligned/ un-aligned DNA sequences input
# translating DNA input according to specific codon-table
# aligning the DNA sequences
# web sites with relevant info for bioPerl modoules in use:
# http://aspn.activestate.com/ASPN/CodeDoc/bioperl/Bio/SeqIO.html
# http://doc.bioperl.org/releases/bioperl-1.4/Bio/Tools/CodonTable.html
# http://www.bioperl.org/Core/Latest/bptutorial.html#iii_3_manipulating_sequences
#   search: Translating a nucleotide sequence from start to end

use Bio::SeqIO;
use Bio::Perl;
use Bio::Tools::CodonTable;
use Math::BigInt;
use strict;

# input: $DNAinput : DNA sequences in FASTA format. We assume that all are in Fasta format.
# output: $DNA_aligned_file : a file where all the input DNA sequences are aligned.


# all routines return list. first argument: "sys"/"user"/"ok". sys - a system error, should not be displayed in html, only in log file.
#    user - user error, should be displayed in html and in log file.  ok - no problems.


sub DNA_align{
    #### INPUT PATHs
    my $DNAinput = shift;	    # the DNA unaligned file
    my $muscle_output = shift;	    # a file that will hold the MSA created by Muscle
    my $workingDir = shift;
    my $DNA_aligned_file = shift;   # The final output, DNA aligned file
    my $outputHtml = shift;         # the html file, used to report of errors to the user
    my $muscle = shift;             # the path to muscle on the running machine
    my $WWWdir = shift;
    
    $muscle_output = $workingDir.$muscle_output;
    $DNA_aligned_file = $workingDir.$DNA_aligned_file;
    
    #### OTHER INPUTS:
    my $codonTableIndex = shift;    # Number that denotes codon table    
    my $ref_seq_name = shift;       # an pointer to the sequence name array, to store DNA sequences names
    
    my %DNA_AA_seq = ();
    my $translatedAA_file = $workingDir."trans_DNA.txt";
    my $muscle_log = $workingDir."muscle.log";    # a file that will hold the peptides sequences
    my $xCodonfile = $workingDir."xCodons.html";    # an html file that will hold all the codons that were not translated
    
    my @ans = &translate_DNA_to_AA($DNAinput, $translatedAA_file, $codonTableIndex, $xCodonfile, $outputHtml, \%DNA_AA_seq, $WWWdir, $ref_seq_name);
    if ($ans[0] eq "ok" ){
        @ans = "";
        # after the AA file is ready, submit it to muscle
        system 'echo "('.$muscle.' -in '.$translatedAA_file.' -out '.$muscle_output.' -log '.$muscle_log.' -quiet)" | /bin/tcsh';
        if (-z $muscle_output){
            @ans = ("sys", "\ncodonAlign::DNA_align : muscle output file $muscle_output is empty. look in $muscle_log for more information\n");
        }
        @ans = &AA_to_DNA_aligned($muscle_output, $DNA_aligned_file, \%DNA_AA_seq);
    }
    else{
        return @ans;
    }
    return @ans;
}
#-------------------------
# gets the DNA aligned file, check if it is legal, removes '*' from the end of sequences.
# translate each legal sequence and write it to AA file
# make sure that the input file have write permission!
sub DNA_checkLegal_and_crate_AAFile{
    
    my $DNAinput = shift;
    my $DNAFinal = shift;
    my $workingDir = shift;
    my $codonTableIndex = shift;
    my $ref_seq_name = shift;
    my $AA_file = shift;
    my $outputHtml = shift;
    my $WWWdir = shift;
    my $print_last_stop_codon = "no";
    
    $DNAFinal = $workingDir.$DNAFinal;
    $AA_file = $workingDir.$AA_file;
    
    my $counter = 0;
    my ($DNASequence, $ter_mark_found, $DNASequenceName, @legalDNA, $SeqLength, $AASeq);
    my $inFile  = Bio::SeqIO->new('-file' => "$DNAinput" , '-format' => 'Fasta');
    my $firstSeq = 0;
    my $xFlag = "no";
    
    my $xCodonfile = $workingDir."xCodons.html";    # an html file that will hold all the codons that were not 
    
    open X_CODON, ">$xCodonfile";
    print X_CODON "<html><table width=50%>\n<tr><td>Sequence Name<\/td><td>Codon Position<\/td><td>Codon<\/td><\/tr>\n";
    close X_CODON;
    
    while ( my $seqObj = $inFile->next_seq() ) {    
        $counter++;
        $DNASequence = 	$seqObj->seq();
        # marking the first sequence, in roder to compare its length to all other sequences.
        if ($firstSeq == 0) {
            $firstSeq = 1;
            $SeqLength = length($DNASequence);            
        }
        # if the last charachter in the sequence is '*', we cut it from the sequence. we put a flag, so we can inform the user if needed.
        if ($DNASequence =~ m/\*$/){
            chop($DNASequence);
            $ter_mark_found = "yes";
        }
        $DNASequenceName = $seqObj->display_id();
        if ($DNASequenceName eq ""){
            $DNASequenceName = "Seq_".$counter;	
        }
        # checking from the second sequence onward that the length matches the first sequence length
        if ($firstSeq==1){
            if ($SeqLength != length($DNASequence)){
                unlink $xCodonfile;
                return ("user", "When submitting <b>codon-aligned</b> file, All your sequences must be of the same length. The sequence $DNASequenceName length is: ".length($DNASequence).". Previous sequence was found to have length of ".$SeqLength.". Please correct your input file ans resubmit your query.\n");
            }
        }
        #inserting the DNASequenceName to the sequence name array.
        #first we check that this seq. name doesn't already exists.
        for (my $i=1; $i<=$#$ref_seq_name; $i++){
            if ($DNASequenceName eq $ref_seq_name->[$i]){
                unlink $xCodonfile;                
                return ("user", "The Sequence name $DNASequenceName appears in your DNA input file more than once. Please make sure that each sequence name in your input files is unique.\n");
            }
        }
        $ref_seq_name->[$counter] = $DNASequenceName;
        #check that the DNA seq is legal. if not: output message to user
        # if a codon stop was found -  We output a message to the user.
        @legalDNA = &check_DNA_seq($DNASequence, $DNASequenceName, $ter_mark_found, $codonTableIndex, "yes");
        unless ($legalDNA[0] eq "yes"){
            if ($legalDNA[0] eq "no"){
                unlink $xCodonfile;
                return ("user", $legalDNA[1]);
            }
            elsif ($legalDNA[0] eq "fix"){
                if ($print_last_stop_codon eq "no"){
                    open HTML_OUT, ">>$outputHtml";
                    print HTML_OUT "\n<p><ul><li><font color='red'><b>Please note:</b></font> $legalDNA[1]. <br>\nThe calculation continues nevertheless.</li></ul></p>\n";
                    close HTML_OUT;
                    $print_last_stop_codon = "yes";
                }
                $DNASequence = substr $DNASequence, 0, ((length $DNASequence) -3);
            }
        }
        # dna Seq is legal - we print it to the final DNA aligned file
        # the reason we open and close this file each iteration, is that other wise some nasty remarks are sent to the file
        unless (open OUT, ">>$DNAFinal"){
            unlink $xCodonfile;            
            return ("sys", "\ncodonAlign::DNA_checkLegal_and_crate_AAFile : can't open file $DNAFinal for writing\n");
        }
        print OUT '>'.($counter)."\n$DNASequence\n";
        close OUT;
        
        # creating file with all translated dna sequences
        ($xFlag, $AASeq) = &translate_sequence($DNASequence, $DNASequenceName, $codonTableIndex, $xFlag, $xCodonfile);
        unless(open AA_FILE, ">>".$AA_file){
            unlink $xCodonfile; 
            unlink $DNAFinal;             
            return ("sys", "\ncodonAlign::DNA_checkLegal_and_crate_AAFile : can't open file $AA_file for writing\n");
        }
        print AA_FILE ">$counter\n$AASeq\n";
        close AA_FILE;
        
    }
    $inFile->close();
    open X_CODON, ">>$xCodonfile";
    print X_CODON "<\/table><\/html>\n";
    close X_CODON;
    if ($xFlag eq "yes"){
        chmod 0644, $xCodonfile;
        open HTML_OUT, ">>$outputHtml";
        print HTML_OUT "\n<p><ul><li><font color='red'><b>Warning:</b></font> Unknown codons were found in your query file and were translated to 'X'. Please look <a href=\"$WWWdir"."xCodons.html\" target=\"xcodon\">here</a> for details. Please note that many 'X' signs in your translated DNA sequence have an impact on the alignment's quality. <br>\nThe calculation continues nevertheless.</li></ul></p>\n";;
        close HTML_OUT;
    }
    # if there were no X codons, we delete the xCodon file since it is useless
    else{
        unlink $xCodonfile;
    }    
    #system ("rm -f $DNAinput; mv $outFile $DNAinput");
    return ("ok");
}
#-------------------------

# building a hash that will hold for each sequence name, a pair of DNA-AA seq put in an array
# if a seqName was not given, give a tentative name.
# create an AA file that will be submitted to MUSCLE
sub translate_DNA_to_AA{
    
    my $input_file = shift;
    my $output_file = shift;
    my $codonTableIndex = shift;    # Number that denotes codon table
    my $xCodonfile = shift;
    my $outputHtml = shift;
    my $ref_DNA_AA_seq = shift;     # a reference to hash table to hold DNA and AA seqs
    my $WWWdir = shift;
    my $ref_seq_name = shift;
	my $OutNameFormat=shift; # {num|seqNum} - The format of name in the coded file (num=only the seq number, seqNum=seq%04u [HoT]) default=num;
	if (!defined $OutNameFormat){$OutNameFormat="NUM";}
	else {$OutNameFormat=uc($OutNameFormat);}
	if ($ref_seq_name eq ""){my @ref_seq_names=();$ref_seq_name=\@ref_seq_names;}

    my $counter = 0;
    my ($DNASequence, $DNASequenceName, $AASeq, @legalDNA);
    my $xFlag = "no";
    my $ter_mark_found = "";
    my $print_last_stop_codon = "no";
        
    my $codonTable_obj  = Bio::Tools::CodonTable -> new ( -id => $codonTableIndex );
    unless (open OUT_AA, ">$output_file"){
        return ("sys", "\ncoldonAlign::translate_DNA_to_AA : can't open file $output_file for writing\n");
    }
    open X_CODON, ">$xCodonfile";
    print X_CODON "<html><table width=50%>\n<tr><td>Sequence Name<\/td><td>Codon Position<\/td><td>Codon<\/td><\/tr>\n";
    close X_CODON;
    my $inFile  = Bio::SeqIO->new('-file' => "$input_file" , '-format' => 'Fasta');
    while ( my $seqObj = $inFile->next_seq() ) {    
        $counter++;
        $DNASequence = 	$seqObj->seq();
        # if the last charachter in the sequence is '*', we cut it from the sequence. we put a flag, so we can inform the user if needed.
        if ($DNASequence =~ m/\*$/){
            chop($DNASequence);
            $ter_mark_found = "yes";
        }
        $DNASequenceName = $seqObj->display_id();
        if ($DNASequenceName eq ""){
            $DNASequenceName = "Seq_".$counter;	
        }
        # if a codon stop was found - we remove it, as it interrupts Selecton calculation. 
        @legalDNA = &check_DNA_seq($DNASequence, $DNASequenceName, $ter_mark_found, $codonTableIndex, "no");
        unless ($legalDNA[0] eq "yes"){
            if ($legalDNA[0] eq "no"){
                close OUT_AA;
                unlink $xCodonfile;
                return ("user", $legalDNA[1]);
            }
            elsif($legalDNA[0] eq "fix"){
                if ($print_last_stop_codon eq "no"){
                    open HTML_OUT, ">>$outputHtml";
                    print HTML_OUT "\n<p><ul><li><font color='red'><b>Please note:</b></font> $legalDNA[1]. <br>\nThe calculation continues nevertheless.</li></ul></p>\n";
                    close HTML_OUT;
                    $print_last_stop_codon = "yes";
                }
                $DNASequence = substr $DNASequence, 0, ((length $DNASequence) -3);
                
            }
        }
        # In case the user submitted a DNA sequence with gaps, we exit and output appropriate message
        if($DNASequence=~ /-/){
            return ("user","You have chosen to upload a DNA file which is not codon-aligned. Despite that, the sign \'-\' (which stands for a gap) was found in your input file in sequence \"$DNASequenceName\".<br>\nIf your file is codon-aligned, please use the second box for DNA codon-aligned sequences.<br>\nOtherwise, please remove \'-\' signs from your input file and resubmit your query.");
        }
        # The actual translation
        ($xFlag, $AASeq) = &translate_sequence($DNASequence, $DNASequenceName, $codonTableIndex, $xFlag, $xCodonfile);
        if ($AASeq =~ m/\*$/){
            chop($AASeq);
        }
		if ($OutNameFormat eq "SEQNUM")
		{
			my $seqNum=sprintf('seq%04u',$counter-1);
			$ref_DNA_AA_seq->{($seqNum)} = [$DNASequence, $AASeq];    #DNA_AA_seq{$DNASequenceName}[0] = $DNASequence
                                                                      #$DNA_AA_seq{$DNASequenceName}[1] = $AASeq
		}
		else
		{
			$ref_DNA_AA_seq->{($counter)} = [$DNASequence, $AASeq];    #DNA_AA_seq{$DNASequenceName}[0] = $DNASequence
			                                                           #$DNA_AA_seq{$DNASequenceName}[1] = $AASeq
		}
        for (my $i=1; $i<=($#$ref_seq_name+1); $i++){
            if ((defined $ref_seq_name->[$i]) and ($DNASequenceName eq $ref_seq_name->[$i])){
                return ("user", "The Sequence name \"$DNASequenceName\" appears in your DNA input file more than once. Please make sure that each sequence name in your input files is unique and re-submit your query.\n");
            }
        }                                                                        
        $ref_seq_name->[$counter] = $DNASequenceName;
		if ($OutNameFormat eq "SEQNUM")
		{
			my $seqNum=sprintf('seq%04u',$counter-1);
			print OUT_AA ">".($seqNum)."\n$AASeq\n";     #DNA_AA_seq{$DNASequenceName}[0] = $DNASequence
                                                          #$DNA_AA_seq{$DNASequenceName}[1] = $AASeq
		}
		else
		{
			print OUT_AA ">".($counter)."\n$AASeq\n";        
		}
    }
    $inFile->close();
    close OUT_AA;
    
    open X_CODON, ">>$xCodonfile";
    print X_CODON "<\/table><\/html>\n";    
    close X_CODON;
    if ($xFlag eq "yes"){
        chmod 0644, $xCodonfile;
        open HTML_OUT, ">>$outputHtml";
		if (!defined $WWWdir){$WWWdir="";}
        print HTML_OUT "\n<p><ul><li><font color='red'><b>Warning:</b></font> Unknown codons were found in your query file and were translated to 'X'. Please look <a href=\"$WWWdir"."xCodons.html\" target=\"xcodon\">here</a> for details. Please note that many 'X' signs in your translated DNA sequence have an impact on the alignment's quality. <br>\nThe calculation continues nevertheless.</li></ul></p>\n";;
        close HTML_OUT;
    }
    # if there were no X codons, we delete the xCodon file since it is useless
    else{
        unlink $xCodonfile;
    }
    return ("ok");
}
#-------------------------

# reading Muscle output and returning to the original DNA sequences, print them aligned to the final output file
sub AA_to_DNA_aligned{
    
    my $input_AA_file = shift;
    my $output_DNA_file = shift;
    my $ref_DNA_AA_seq = shift;
    
    my $DNASequenceName;
    my $AA_seq_pointer;
    my $line;
    
    unless (open OUT_ALIGNED, ">$output_DNA_file"){
        return ("sys", "\ncodonAlign::AA_to_DNA_aligned : can't open file $output_DNA_file for writing\n");
    }
    unless (open MUSCLE, $input_AA_file){
        return ("sys", "\ncodonAlign::AA_to_DNA_aligned : can't open file $input_AA_file for reading\n");
    }
    while (<MUSCLE>){
        chomp;
        if (/>(.+)/){
            $DNASequenceName = $1;
            print OUT_ALIGNED $_."\n";
            $AA_seq_pointer = 0;
            $line = 0;
        }
        else{
            $line++;
            my @AA_Muscle = split(//, $_);
            foreach (@AA_Muscle){
                if(/-/){
                    print OUT_ALIGNED "---";
                }
                elsif(/\w/){
                    if($_ eq (substr($ref_DNA_AA_seq->{$DNASequenceName}[1], $AA_seq_pointer, 1)) ){
                        print OUT_ALIGNED substr($ref_DNA_AA_seq->{$DNASequenceName}[0], 0, 3);
                        $ref_DNA_AA_seq->{$DNASequenceName}[0] = substr($ref_DNA_AA_seq->{$DNASequenceName}[0], 3);
                        $AA_seq_pointer++;
                    }
                    else{
                        close MUSCLE;
                        close OUT_ALIGNED;
                        return ("sys", "\ncodonAlign::AA_to_DNA_aligned : In seq name: $DNASequenceName Read from muscle file char: $_ at index: $AA_seq_pointer line: $line. In hash found ".substr($ref_DNA_AA_seq->{$DNASequenceName}[1], $AA_seq_pointer, 1),"\n");                    
                    }
                }
            }
            print OUT_ALIGNED "\n"; 
        }               
    }
    close MUSCLE;
    close OUT_ALIGNED;
    return ("ok");
}
#-------------------------

# Check if DNA sequence is legal:
# 1. it is divisible by 3
# 2. it has no stop codon or * sign in its middle   STOP CODONS are detected according to the chosen codon convert table
# input: DNA sequence
# output: "yes" if all tests are OK, otherwise - a string that describes the input problem

sub check_DNA_seq{
    my $inputDNA = shift;
    my $DNAinputName = shift;
    my $ter_mark = shift;
    my $tableCodonIndex = shift;
    my $is_dna_aligned = shift;
        
    my @ans = ("yes", "yes");
    my $codon;
    
    my $codonTable_obj  = Bio::Tools::CodonTable -> new ( -id => $tableCodonIndex );
    
    my $seq_length = length($inputDNA);
    my $val = Math::BigInt->new($seq_length);    
    if ($val->bmod(3) != 0){
        if ($ter_mark eq "yes"){ # in case an earlier * sign was cut from the sequence, we inform the user
            $ter_mark = '(without the last * sign)';
        }
        @ans = ("no", "The sequence $DNAinputName $ter_mark is of length $seq_length, which is not divisible by 3.");
	return @ans;
    }    
    my $i =0;
    while ($i<$seq_length-2){
	$codon = substr($inputDNA, $i, 3);
	if (((!($codonTable_obj->is_unknown_codon($codon)) && $codonTable_obj->is_ter_codon($codon)) || $codon =~ m/\*/) && $i <= $seq_length-6) {
	    @ans = ("no", "A Stop codon ,\"$codon\", was found in sequence $DNAinputName in position ".($i+1).". Please verify that there are no internal stop-codons in your sequences.");
	    return @ans;
	}
        elsif($i==$seq_length-3 && $codon =~ m/-/){
        }
        # in case the DNA input file was aligned, we ask the user to remove stop codons from the end of the sequences, as some of his seuqneces might have stop codons and some are not - and we don't want to delete it for him (to decide for him whether to remove, or to put a gap etc.)
        elsif($is_dna_aligned eq "yes" && $i==$seq_length-3 && $codon =~ m/[ATCG]{3}/i && $codonTable_obj->is_ter_codon($codon)){
            @ans = ("no", "Please remove the Stop Codon \"$codon\" from you sequence $DNAinputName.");
            return @ans;            
        }
	$i+=3;
    }    
    return @ans;
}

sub translate_sequence{
    
    my $DNASequence = shift;
    my $DNASequenceName = shift;
    my $codonTableIndex = shift;
    my $xFlag = shift;
    my $xCodonfile = shift;
    my ($codon, $AA);
    
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
    return ($xFlag, $AASeq);
}

1;
