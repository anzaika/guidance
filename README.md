                GUIDANCE:  GUIDe tree based AligNment ConfidencE

GUIDANCE is a software package for aligning biological sequences (DNA or
amino acids) using either MAFFT, PRANK, or CLUSTALW, and calculating
confidence scores for each column, sequence and residue in the alignment.

URL: http://guidance.tau.ac.il/

Authors: Osnat Penn, Eyal Privman, Haim Ashkenazy, Itamar Sela, Giddy Landan, Dan Graur, and Tal Pupko.

When using the GUIDANCE2 algorithm please cite:
-----------------------------------------------------------
Sela, I., Ashkenazy, H., Katoh, K. and Pupko, T. (2015)
GUIDANCE2: accurate detection of unreliable alignment regions accounting for the uncertainty of multiple parameters.
Nucleic Acids Research, 2015 Jul 1; 43: W7-W14.; doi: 10.1093/nar/gkq443

Landan, G., and D. Graur. (2008).
Local reliability measures from sets of co-optimal multiple sequence alignments.
Pac Symp Biocomput 13:15-24

When using the GUIDANCE algorithm please cite:
-----------------------------------------------------------
Penn, O., Privman, E., Landan, G., Graur, D. and Pupko, T. (2010).
An alignment confidence score capturing robustness to guide-tree uncertainty. 
Molecular Biology and Evolution, 2010 Aug;27(8):1759-67; doi:10.1093/molbev/msq066

When using the HoT algorithm please cite:
-----------------------------------------------------------
Landan, G., and D. Graur. (2008).
Local reliability measures from sets of co-optimal multiple sequence alignments.
Pac Symp Biocomput 13:15-24

Installation
============

  1. Unpack the archive by typing:
        % tar -xzf guidance.v2.01.tgz

  2. Compile the package by typing:
        % cd guidance.v2.01
        % make
     (Running `make' takes a while)

  3. Check if you have the desired alignment program installed:

     MAFFT: Type "mafft" and check that you have version 6.712 or newer.
      * Else download and install MAFFT from:  http://mafft.cbrc.jp/alignment/software/

     PRANK: Type "prank" and check that you have it insalled
      * Else download and install PRANK from:  http://www.ebi.ac.uk/goldman-srv/prank/prank/

     CLUSTALW: Type "clustalw" and check that you have it insalled
      * Else download and install CLUSTALW from:  http://www.ebi.ac.uk/Tools/clustalw2/index.html

     MUSCLE: Type "muscle" and check that you have it insalled
      * Else download and install MUSCLE from:  http://www.drive5.com/muscle/index.htm

	 PAGAN: Type "pagan" and check that you have it installed
      * Else download and install PAGAN from:  http://code.google.com/p/pagan-msa/
      * In case PAGAN is used not for user provided alignment (using --msaFile), mafft is also required to be installed.


  4. GUIDANCE also uses Perl, BioPerl and Ruby:
                * Type "perl -v" and check that you Perl installed.
                        Else download and install it from: http://www.perl.org/
                * Type "perl -e 'use Bio::SeqIO'" to check that you have BioPerl.
                        Else download and install it from: http://www.bioperl.org/
                * Type "ruby -version" to check that you have ruby.
                        Else download and install it from: http://www.ruby-lang.org/en/

Usage
=====

Run the Perl script: guidance.v2.01/www/Guidance/guidance.pl
(Note that you cannot move this script out of its directory, because it uses relative paths to other files in other directories. Sorry)
GUIDANCE uses flags in the command line arguments: (for help, type: "perl guidance")

USAGE:perl guidance.pl --seqFile SEQFILE --msaProgram [MAFFT|PRANK|CLUSTALW|MUSCLE] --seqType [aa|nuc|codon] --outDir FULL_PATH_OUTDIR

Required parameters:
  --seqFile          Input sequence file in FASTA format
  --seqType      Sequence type may be either of: nuc (nucleotides), aa (amino acids),
                 or codon (nucleotides that will be treated as whole codons)
  --msaProgram   The alignment program - may be either MAFFT, PRANK, CLUSTALW or  MUSCLE
  --outDir       The output directory were all output files will be created [please provide full and not relative path]
                 (will be created automatically)

Optional parameters:
  --program      The confidence measure may be GUIDANCE2, GUIDANCE or HoT. default=GUIDANCE2
  --bootstraps   Number of bootstrap iterations. default=100
  --genCode      Genetic code for use in codon sequence. default=1
                 1>  Nuclear Standard
                 15> Nuclear Blepharisma
                 6>  Nuclear Ciliate
                 10> Nuclear Euplotid
                 2>  Mitochondria Vertebrate
                 5>  Mitochondria Invertebrate
                 3>  Mitochondria Yeast
                 13> Mitochondria Ascidian
                 9>  Mitochondria Echinoderm
                 14> Mitochondria Flatworm
                 4>  Mitochondria Protozoan
  --outOrder     May be either aligned or as_input. default=aligned
  --msaFile      Input alignment file - not recommended, see documentation online at: guidance.tau.ac.il
  --seqCutoff    Sequence confidence cutoff between 0 to 1. default=0.6
  --colCutoff    Columnd confidence cutoff between 0 to 1. default=0.93
  --mafft        Path to mafft executable. default=mafft
  --prank        Path to prank executable. default=prank
  --clustalw     path to clustalw executable. default=clustalw
  --muscle       path to muscle executable. default=muscle
  --pagan        path to pagan executable. default=pagan
  --ruby         path to ruby executable. default=ruby
  --dataset      Unique name for the Dataset - will be used as prefix to outputs (default=MSA)
  --MSA_Param:   Passing parameters for the alignment program e.g -F to prank. To pass parameter containning '-' in it, add \ before each '-' e.g. \-F for PRANK
  --proc_num:    Number of processors to use (default=1)

EXAMPLES:

>perl guidance.pl --seqFile protein.fas --msaProgram MAFFT --seqType aa --outDir /somedir/protein.guidance
Will align the amino acid sequences in the fasta file "protein.fas" using MAFFT and output all results to the diretory "/somedir/protein.guidance"

>perl guidance.pl --seqFile codingSeq.fas --msaProgram PRANK --seqType codon --outDir /somedir/codingSeq.guidance --genCode 2 --bootstraps 30
Will align the codon sequences in the fasta file "codingSeq.fas" using PRANK after translation using the vertebrate mitochondrial genetic code and output all results to the diretory "/somedir/codingSeq.guidance". Only 30 bootstrap iterations will be done instead of the default 100 (cut run-time by a factor of 3)

Copyrights
==========

    * To modify the code, or use parts of it for other purposes, permission should be requested. Please contact Tal Pupko: talp@post.tau.ac.il
    * Please note that the use of the GUIDANCE program is for academic use only
