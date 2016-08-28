#!/usr/bin/perl

use strict;
use lib "/bioseq/bioSequence_scripts_and_constants";
use GENERAL_CONSTANTS;

#*******************************************************************************
# The script parse the output of PSI-BLAST, eliminate the redundant sequences 
# and write the homologues in a file ordered by their Evalue. 
#*******************************************************************************

my $MaxNumHomol = shift; #max number of homologues to extract
my $Cutoff = shift; # the Evalue cutoff
my $server = shift; # ConSurf or ConSeq
my $iterations = shift; # number of PSI-BLAST iterations
my $workingDir = shift;
my $consurf_unique_seqs_file = $workingDir.shift;
my $hom_seq_file = $workingDir.shift;
my $alout = shift;
my $alout_html = $workingDir.shift;
my $target_file = $workingDir.shift;
my $OutputHtml = shift;
my $out_error_file = $workingDir.shift;
my $SeqName = shift; # when the server is ConSeq

#### variables

my @database = ();
my @sort_database = ();
my $BlastHomol = 0; # the number of all the hits + the original seq. 
my $IdenticalSeq = 0;
my $PairLines = 0;

################################################################################################

### extract the sequences from the BLAST output
&Extract_BlastP_Seq;

### add the query sequence to the array of homologues 
&Add_Target_Seq;

### sort the homologues by their Evalue
&Sort_Evalue;

### eliminate the redundant sequences
&Double_Out;

### print the homologues in fasta format
&Print_Seq;

### Number Blast Hits
my $numbered_blast=number_Blast_Seq($hom_seq_file,$iterations);
rename ($hom_seq_file,"$hom_seq_file.WO_Num");
rename ($numbered_blast,$hom_seq_file);
#######################################################
# extract the sequences from the BLAST output
sub Extract_BlastP_Seq {
    
    unless (open HOM_SEQ_FILE, $hom_seq_file){
        &error_and_exit("sys", "blastpgp_to_fasta.pl : cannot open file $hom_seq_file for reading $!");
    }     
    my $SeqNo = -1;
    my $SeqHomol = 0;
    my $SeqName = "";
    my $Evalue;
    my $Annotation_Line="";
    my @file = <HOM_SEQ_FILE>;
    my $index = 0;

    ### deal with PSI-BLAST iterations
    if ($iterations > 1){	
        # find the last iteration that have results
        for (my $i=0 ; $i <= $#file ; $i++){    
            if ($file[$i] =~ /Results\s+from\s+round/){
                $index = $i;
            }
        }	
    }

    for ($index ; $index <= $#file ; $index++){
               
	if ($file[$index] =~ /^>(?:SPROT:)?(?:\w\w\|)?(\w+)/){
	        $SeqName = $1;
		$Annotation_Line=$file[$index];
		chomp ($Annotation_Line);
		$index++;
		while ($file[$index] !~ /Length =/)
		{
			$Annotation_Line=$Annotation_Line.$file[$index];
			chomp $Annotation_Line;
			$Annotation_Line=~ s/\s+/ /g; 
			$index++;
		}
		
	#if ($file[$index] =~ /^>(\S+)\|(\S+)\w+/){
        #$SeqName = $2;
	    
	    $database[++$SeqNo]{SeqName} = $SeqName;
	    $SeqHomol = 0; #for the same sequence
        }
	elsif ($file[$index] =~ /Expect\s+=\s+(\S+)/){

	    $Evalue = $1;
	    $SeqHomol++;

	    if ($SeqHomol > 1) {

		$database[++$SeqNo]{SeqName} = $SeqName . "_" . $SeqHomol;
	    }
	    $database[$SeqNo]{Evalue} = $Evalue;
	    if ($Annotation_Line =~ /\|([^\|]+)$/)
	    {
		$database[$SeqNo]{SeqDesc}=$1;
	    }
	}
	elsif ($file[$index] =~ /^Sbjct:\s+\d+\s+([\w\-]+)/){
	            
            $database[$SeqNo]{Sequence} .= $1;
	    # delete the gaps
	    $database[$SeqNo]{Sequence} =~ s/-//g;
        }	
    } 
    close (HOM_SEQ_FILE);
}

#########################################################################
# Add the target sequence to the first place in @database
sub Add_Target_Seq {

    my %hash = ();
    unless (open SEQ, $target_file){
        &error_and_exit("sys", "blastpgp_to_fasta.pl : Can\'t open the file $target_file $!");
    } 

    while (<SEQ>){
	# extract the name of the query sequence    
        if ($_ =~ /^>(\S+)\s+/){              
            $hash{SeqName} = $1;
        }    
        # extract the query sequence
        elsif ($_ =~ /^([\w\-]+)\s*$/) {
            $hash{Sequence} .= $1;
        }
    }    my @file = <HOM_SEQ_FILE>;
    my $index = 0;
    close SEQ;

    # for ConSeq server
    if ($server eq "ConSeq"){
        $hash{SeqName} = $SeqName;
    }

    unshift @database, \%hash;
    $BlastHomol = $#database + 1; # The number of homologues found by BLAST
}

#########################################################################
# Eliminate the redundant sequences
sub Double_Out {
    for (my $i = 0; $i< $BlastHomol; $i++) {
        for (my $j = $i+1; $j< $BlastHomol; $j++) {            
            ### to delete identical sequences
            if ($sort_database[$i]{Sequence} eq $sort_database[$j]{Sequence} or $sort_database[$i]{SeqName} eq $sort_database[$j]{SeqName}){
                if ($sort_database[$j]{delete} != 1){
                    $sort_database[$j]{delete} = 1;
                    $IdenticalSeq++;
                }
            } 
        }
    }     
}

##########################################################
# sub to print the names and sequences in fasta format
sub Print_Seq {	

    my $CleanHomol = $BlastHomol - $IdenticalSeq;
    unless (open (UNIQUESEQ,">$consurf_unique_seqs_file")){
        &error_and_exit("sys", "blastpgp_to_fasta.pl : Can\'t open the file $consurf_unique_seqs_file for writing $!");
    } 
    print UNIQUESEQ $CleanHomol;
    close UNIQUESEQ;
    #####################
    ### to set all homol
    if ($MaxNumHomol eq "all" or $MaxNumHomol > $CleanHomol) {	
        $MaxNumHomol = $CleanHomol;
    }
    #PRINTS SEQUENCES
    unless (open HOM_SEQ_FILE_OUT, ">$workingDir".$alout){
        &error_and_exit("sys", "blastpgp_to_fasta.pl : Can\'t open the file ".$workingDir."$alout for writing $!");
    }
    
    # print the target-seq
    print HOM_SEQ_FILE_OUT ">$sort_database[0]{SeqName} \n$sort_database[0]{Sequence} \n";

    # print only the unique sequences 
    my $UniqueHomol = 1;   
    for (my $i = 1; $i<= $BlastHomol; $i++) {
        last if ($UniqueHomol >= $MaxNumHomol);        
        if (($sort_database[$i]{delete} != 1) and ($sort_database[$i]{Evalue} <= $Cutoff)) {
            print HOM_SEQ_FILE_OUT ">$sort_database[$i]{SeqName} | Blast_Hit_$i | $sort_database[$i]{SeqDesc}\n$sort_database[$i]{Sequence} \n";
    	    $UniqueHomol++;
        }	
    }
    close (HOM_SEQ_FILE_OUT);

#Buid the _homol.html file

   #PRINTS SEQUENCES
    if (!open HOM_SEQ_FILE_OUT_HTML, ">$alout_html"){
        &error_and_exit("sys", "blastpgp_to_fasta.pl : Can\'t open the file $alout_html for writing $!");
    }  
       
    print HOM_SEQ_FILE_OUT_HTML "<HTML> \n  <HEAD> \n ";
    print HOM_SEQ_FILE_OUT_HTML  " </HEAD> \n <body> \n <TABLE> \n";

    # print the target-seq
    
    print HOM_SEQ_FILE_OUT_HTML " <TR> <TD><FONT SIZE=2>>$sort_database[0]{SeqName} </FONT> </TD> </TR> <TR> <TD><FONT SIZE=2>$sort_database[0]{Sequence}</FONT> </TD> </TR> \n";

    # print only the unique sequences 
    my $UniqueHomol = 1;   
    for (my $i = 1; $i<= $BlastHomol; $i++) {
        last if ($UniqueHomol >= $MaxNumHomol);
        if (($sort_database[$i]{delete} != 1) and ($sort_database[$i]{Evalue} <= $Cutoff)) {
            print HOM_SEQ_FILE_OUT_HTML " <TR> <TD> ><A HREF=\"http://www.expasy.org/cgi-bin/niceprot.pl?$sort_database[$i]{SeqName} \" ><FONT SIZE=2>$sort_database[$i]{SeqName}</A> | Blast_Hit_$i | $sort_database[$i]{SeqDesc}<BR></FONT> </TD> </TR> <TR> <TD><FONT SIZE=2>$sort_database[$i]{Sequence}</FONT></TD> </TR> \n ";    
            $UniqueHomol++;
        }	
    }
    print HOM_SEQ_FILE_OUT_HTML "</TABLE> \n  </BODY> \n  </HTML>";
    close (HOM_SEQ_FILE_OUT_HTML);

    #*****************************
    ##### PRINT ON output.html NO OF HOMOLOGUES

    # If there are not enough homologues, write error message and exit    
	if ($CleanHomol < 5) {
        open OUTPUT, ">>$OutputHtml";    
        if ($CleanHomol == 1){
            &error_and_exit("<b>There is only 1 <A HREF=$alout TARGET=blast>unique hit</A>. The minimal number of sequences required for the calculation is 5.<br>You can try to run $server with a multiple sequence alignment file of your own.</b></p>\n", "blastpgp_to_fasta: Print_Seq : found only 1 homoloug\n");
        }
        else {    
            &error_and_exit("<b>There are only $CleanHomol <A HREF=$alout TARGET=blast>unique PSI-BLAST hits</A>. The minimal number of sequences required for the calculation is 5.<br>You can try to run $server with a multiple sequence alignment file of your own.</b></p>\n", "blastpgp_to_fasta: Print_Seq : found less than 5 homolougs\n");
        }    
    }

    open OUTPUT, ">>" .$OutputHtml;    
    # There are more homologues than the cutoff
    if ($MaxNumHomol < $CleanHomol) {		
        # There are identical sequences
        if ($IdenticalSeq > 0){
            print OUTPUT "\n<p><ul><li>There are $BlastHomol PSI-BLAST hits, $CleanHomol of them are unique sequences.<br>The calculation is performed on the $MaxNumHomol sequences with the lowest E-value.</li></ul></p>\n";
			print "\n<p><ul><li>There are $BlastHomol PSI-BLAST hits, $CleanHomol of them are unique sequences.<br>The calculation is performed on the $MaxNumHomol sequences with the lowest E-value.</li></ul></p>\n";
        }        
        # All the sequences are unique
        else {
            print OUTPUT "\n<p><ul><li>There are $BlastHomol unique PSI-BLAST hits.<br>The calculation is performed on the $MaxNumHomol sequences with the lowest E-value.</li></ul></p>\n";
        }
    }

    # There are less homologues than the cutoff.
    else {
        # There are identical sequences
        if ($IdenticalSeq > 0){
            # There are less than 10 unique sequences
            if ($CleanHomol < 10){
                print OUTPUT "\n<p><ul><li><font color='red'><b>Warning:</font></b> There are $BlastHomol PSI-BLAST hits, only $CleanHomol of them are unique sequences. The calculation is performed on the $CleanHomol unique sequences, but it is recommended to run $server with a multiple sequence alignment file containing at least 10 homologues.</li></ul></p>\n";
            }
    
            # There are at least 10 unique sequences
            else {
                print OUTPUT "\n<p><ul><li>There are $BlastHomol PSI-BLAST hits, $CleanHomol of them are unique sequences.<br>The calculation is performed on the $CleanHomol unique sequences.</li></ul></p>\n";
            }
        }
    
        # All the sequences are unique
        else {
            # There are less than 10 unique sequences
            if ($CleanHomol < 10){
                print OUTPUT "\n<p><ul><li><font color='red'><b>Warning:</font></b> There are only $CleanHomol PSI-BLAST hits. The calculation continues nevertheless, but it is recommended to run $server with a multiple sequence alignment file containing at least 10 homologues.</li></ul></p>\n";
            }
    
            # There are at least 10 unique sequences
            else {
                print OUTPUT "\n<p><ul><li>There are $BlastHomol unique PSI-BLAST hits.</li></ul></p>\n";
            }
        }
    }
    close OUTPUT;    
}

###############################################################
# The function sorts the homologues array according to the Evalue,
# and deletes the ones with Evalue grater than the cutoff. 
###############################################################
sub Sort_Evalue {

    # move the target-seq (the first one) to the first place in the sorted array.
    push @sort_database, shift(@database);

    # enter the values, sorted by Evalue, to the sorted array
    push @sort_database, (sort { $$a{Evalue} <=> $$b{Evalue} } @database);    

    # delete the sequences with Evalue grater than the cutoff
    # (searching the array from the end to the beginning, except from the first element)
    my $popped = 0;
    my $i = @sort_database - 1;
    for ($i ; $i > 0 ; $i--){
        if ($sort_database[$i]{Evalue} > $Cutoff){	    
            pop @sort_database;
            $popped++;
        }
    }

    # update the total number of sequences
    $BlastHomol -= $popped;
    
}

sub error_and_exit{
    my $h_err = shift;
    my $l_err = shift;
    open ERROR_OUT, ">".$out_error_file;
    print ERROR_OUT "HTML: $h_err\n";
    print ERROR_OUT "LOG: $l_err";
    close ERROR_OUT;
    chmod 0755, $out_error_file;
    exit;
}

sub number_Blast_Seq{
    my $Blast_File = shift;
    my $Iterations = shift;
    my $Numbered_Blast = $Blast_File.".Numbered";
    
   
    unless (open (BLAST,"$Blast_File")){
        &error_and_exit("sys", "blastpgp_to_fasta.pl : Can\'t open the file $Blast_File for reading $!");}
    unless (open (BLAST_NUMBERED,">$Numbered_Blast")){
    	&error_and_exit("sys", "blastpgp_to_fasta.pl : Can\'t open the file $Numbered_Blast for writing $!");}
    
    my $Counter_E_Value_Section=1;
    my $Counter_Align_Section=1;	   
    
    my @file = <BLAST>;
    my $index = 0;     
    ### deal with PSI-BLAST iterations
    if ($iterations > 1){	
        # find the last iteration that have results
        for (my $i=0 ; $i <= $#file ; $i++){
            if ($file[$i] =~ /Results\s+from\s+round/){
                $index = $i;
            }
        }	
    }
    for (my $i=0; $i<=$index;$i++)
    {
    	print BLAST_NUMBERED $file[$i];
    }
    for ($index ; $index <= $#file ; $index++){      
    	my $SeqID=""; 
        if ($file[$index]=~/^([A-Za-z0-9]+)\|(.*)\|/){
		$SeqID=$2;
	    	if ($Counter_E_Value_Section<10) {print BLAST_NUMBERED "$Counter_E_Value_Section   ";}
		elsif ($Counter_E_Value_Section>9 and $Counter_E_Value_Section<100) {print  BLAST_NUMBERED "$Counter_E_Value_Section  ";}
		elsif ($Counter_E_Value_Section>=100) {print BLAST_NUMBERED "$Counter_E_Value_Section ";}
		print BLAST_NUMBERED "$file[$index]";
		$Counter_E_Value_Section++;
            }
        elsif ($file[$index]=~/^>([A-Za-z0-9]+)\|([A-Za-z0-9]+)/){
	     $SeqID=$2;
	     $file[$index]=~s/>//;
	     print BLAST_NUMBERED ">$Counter_Align_Section|".$file[$index]; 
	     $Counter_Align_Section++;
        }
        else
        {
   	  	print BLAST_NUMBERED $file[$index];
	} 
   }
	close (BLAST);
	close (BLAST_NUMBERED);	
	return $Numbered_Blast;
}










