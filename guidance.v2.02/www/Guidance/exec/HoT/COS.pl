#!/usr/bin/perl

#
########
use strict;
use warnings;
use File::Spec::Functions qw(rel2abs);
use File::Basename;

########
# context:
my $my_name=basename($0);
my $basedir=rel2abs('./');
my $debug=defined($ENV{DEBUG_LEVEL}) ? $ENV{DEBUG_LEVEL} : 0;
my $t0=time;
$SIG{TERM} = \&my_sigtrap;
my $usage="\n\nCOS v2.05\n\nUsage:\n".
					"\n  $my_name case_id msa_method seq_type input_fasta_file output_dir [MSA_Status File] [tgz] [msa_program_path] --- [Additional MSA program parameters]\n\n".
  				"\tmsa_method= \t'CLW' : ClustalW2\n".
  				"\t            \t'MFT' : mafft\n".
  				"\t            \t'PRK' : prank\n".
  				"\t            \t             Append 'h' to do just HoT, e.g.: 'CLWh'\n\n".
  				"\tseq_type=   \t'aa' | 'nt'\n\n".
	  			"\tSuccessful output at: (output_dir)/(case_id)_cos_(msa_method)/ or (output_dir)/(case_id)_cos_(msa_method).tgz\n".
		  		"\tError output at: (output_dir)_err/(case_id)_cos_(msa_method)_err/ or (output_dir)_err/(case_id)_cos_(msa_method)_err.tgz\n";


#  				"                'MFT':mafft maximal | 'MAF':mafft partial | 'MA0' mafft minimal\n\n    seq_type= 'aa' | 'nt'\n\n".

#####################################################
my %metcmds;
$metcmds{CLW}='CW2';
$metcmds{CW2}='CW2';
$metcmds{CW3}='CW2';
$metcmds{MFT}='MAF';
$metcmds{MFM}='MAF';
$metcmds{PRK}='PRK';
#####################################################
#####################################################
if ($#ARGV < 4) {
	print $usage;
	exit;
}
#####################################################
## inputs constants etc.
#####################################################
my ($idstr,$met,$seqtype_str,$input_file,$output_dir,$status_file,$tgz,$msa_program_path,@param)=@ARGV;
$met=~/^(.{3})(.?)/;
$met=$1;
my $hot=$2."C";
if (!defined $metcmds{$met}) {
	print "\nERROR: Unknown msa method $met\n$usage";
	exit;
}
my $met_init=\&{'met_init_'.$metcmds{$met}};
my $make_guide_tree=\&{'make_guide_tree_'.$metcmds{$met}};
my $align_seq=\&{'align_seq_'.$metcmds{$met}};
my $align_prof=\&{'align_prof_'.$metcmds{$met}};

$input_file=rel2abs($input_file);
if (! -e $input_file) {
	print "\nERROR: File not found: $input_file\n$usage";
	exit;
}

$output_dir=rel2abs($output_dir);
my $seqtype=($seqtype_str=~/^[nN]/ ? 1 : 0);

#FFU: codon model : seqtype 2
#my $seqtype=0;
#$seqtype=1 if ($seqtype_str=~/^[nN]/);
#$seqtype=2 if ($seqtype_str=~/^[cC]/);

$tgz= defined($tgz) ? ($tgz=~/^[tT]/ ? 1 : 0) : 0;
my $param=join(' ',(@param,' '));
$param=~s/---//;
my ($metver,$prkver);
$met_init->();
my ($prefix)= $input_file=~m%([^/.]*)\.[^/.]*$%;
my @extstr=('.fasta','.atsaf');
my @ht=('h','t');
my @hot=('heads','tails');

my $output_prefix=$idstr.'_cos_'.$met;
my $done_tar_file=$output_dir.'/'.$output_prefix.'.tgz';
my $err_tar_file=$output_dir.'_err/'.$output_prefix."_err.tgz";

my ($cmdstr,$rc);

#####################################################
#-----------------------------------------
# create directoried and files
#=========================================
my $msgstr="--- init $my_name pid=$$ ".localtime()."\nmet=$met\ninfile=$input_file\nseqtype=$seqtype\noutdir=$output_dir\ntgz=$tgz\ndebug=$debug\n".
					 "basedir=$basedir\noutput_id=$output_prefix\n---\n";
`mkdir -p $output_prefix 2>&1`;
chdir($output_prefix);
open(LOGFILE,">>$output_prefix.log");
my $old_fh = select(LOGFILE);$| = 1;select($old_fh); #autoflush
log_print(0,1,$msgstr);
open(TIMEFILE,">>times.txt");
$old_fh = select(TIMEFILE);$| = 1;select($old_fh); #autoflush
print TIMEFILE $msgstr; 
my_cmd("find . -size 0 -delete");# for resuming aborted runs with null files
#####################################################
#
#
#-----------------------------------------
# load input
#=========================================
my %hot_seqs;
fasta2hotseqs($input_file);
if ($hot_seqs{notu}<2) {
 	print "\nERROR: $input_file has less than 2 sequences....\n";
	exit;
}
#=========================================
#-----------------------------------------
# make Guide Tree
#-----------------------------------------
#=========================================
my $infile='input.fasta';
my $treefile='guide_tree.nwk';
str_print($infile,@{$hot_seqs{fasta}[0]});
if (-e $treefile) {
  log_print(1,1,"-$treefile exists");
} 
else {
  log_print(0,1,"-Making guide tree ...\n");
  if ($status_file ne "")
  {
	  open (STATUS,">$status_file.0");
	  print STATUS "<ul><li>Making guide tree</li></ul>\n";
	  close (STATUS);
  }
  $make_guide_tree->($infile,$treefile);
} #if
my %subtrees;
tree2split($treefile);
#=========================================
#-----------------------------------------
# do HoT
#-----------------------------------------
#=========================================
$treefile="in.dnd";
str_print($treefile,$subtrees{tree});
my ($outfile,$outfiler,$msar);
for my $i (0..1) {
  $outfile="hot_".uc($ht[$i]);
  if ($msar=MSA_check2tails($outfile.$extstr[0],1)) {
    log_print(1,1,"-$outfile$extstr[0] exists");
    next;
  }
  if ($i==1 && ($msar=MSA_check2tails($outfile.$extstr[1],1))) {
	  log_print(1,1,"-$outfile$extstr[1] exists");
    str_print($outfile.$extstr[0],$msar);
		next;
	}
  $infile='in'.$extstr[$i];
  str_print($infile,@{$hot_seqs{fasta}[$i]});
  log_print(0,1,"-Making $outfile ...\n");
  if ($status_file ne "")
  {
      open (STATUS,">>$status_file.0");
      print STATUS "<ul><li>Making $outfile</li></ul>\n";
      close (STATUS);
  }
  $msar=$align_seq->($infile,$treefile,$outfile.$extstr[$i]);
  str_print($outfile.$extstr[0],$msar) if ($i==1);
}#for i=0:1
if ($hot=~/^h/i) {
  cleanup(0);
  exit(0);
} 
#=========================================
#-----------------------------------------
# do CoS
#-----------------------------------------
#=========================================
my (@pfiles,@tfiles,$fasta);
for my $i (0..$subtrees{nbr}-1) { # splits
	my $isprof= $i<$subtrees{notu} ? 0 : 1 ;

	#===== check resumed runs
  my $skip=1;
	for my $j (0..$isprof) { #left prof
		for my $j1 (0..1) { #right prof
			for my $k (0..1) { #hot
				$outfile=$subtrees{br}[$i][0]{name}.'_'.$ht[$j].$ht[$j1].uc($ht[$k]);
			  if ($msar=MSA_check2tails($outfile.$extstr[0],1)) {
			    log_print(1,1,"-$outfile$extstr[0] exists");
			    next;
			  }
			  if ($k==1 && ($msar=MSA_check2tails($outfile.$extstr[1],1))) {
			    log_print(1,1,"-$outfile$extstr[1] exists");
          str_print($outfile.$extstr[0],$msar);
			    next;
			  }
			  $skip=0;
			} #for k=0..1 hot
		} #for j1=0..1
	} #for j=0..0/1
	#====================================

	if ($skip) {
  	log_print(0,1,"-Skiping branch ".($i+1)." of $subtrees{nbr} ...\n");
    next;
  }
	print "-Making branch ".($i+1)." of $subtrees{nbr} ...\n";
	if ($status_file ne "")
	{
		open (MSA_STATUS,">$status_file");
		print MSA_STATUS "<ul><li>Making branch ".($i+1)." of $subtrees{nbr}</li></ul>\n";
		close (MSA_STATUS);
	}
	#-----------------------------------------
	# do 2 profiles x hot =4
	for my $j (0..1) {  #left-right
		my %prof=%{$subtrees{br}[$i][$j]};
		if ($j==0 && !$isprof) { # single sequence
			$tfiles[$j]='';
			for my $k (0..1) { # HoT single sequence dummy
				$pfiles[$j][$k]="prof$i\_$j";
        str_print($pfiles[$j][$k].$extstr[$k] , (@{$hot_seqs{fasta}[$k]})[@{$prof{otu}}]);
			} #for k=0..1 HoT
		} # if single seq
		else { # notu>1 , profile
			$tfiles[$j]="tree_$i\_$j.dnd";
      str_print($tfiles[$j],$prof{tree});
			for my $k (0..1) { # HoT
			  $pfiles[$j][$k]="prof$i\_$j$k";
			  $outfile=$pfiles[$j][$k].$extstr[$k];
			  $outfiler=$pfiles[$j][$k].$extstr[1-$k];
			  if ($msar=MSA_check2tails($outfile,1)) {
			    log_print(1,1,"-$outfile exists");
          str_print($outfiler,$msar) if (! -e $outfiler);
          next;
			  }
  		  $infile='in'.$extstr[$k];
        str_print($infile , (@{$hot_seqs{fasta}[$k]})[@{$prof{otu}}]);
			  $msar=$align_seq->($infile,$tfiles[$j],$outfile);
        str_print($outfiler,$msar);
			} #for k=0..1 HoT
		} #else prof
	} #for j=0..1 right-left
	#-----------------------------------------
	# do 8 (or 4 for term) hot profile alignment
	$tfiles[2]="tree_$i.dnd";
  str_print($tfiles[2],$subtrees{br}[$i][2]{tree});
	for my $j (0..$isprof) { #left prof
		for my $j1 (0..1) { #right prof
			for my $k (0..1) { #hot
				$outfile=$subtrees{br}[$i][0]{name}.'_'.$ht[$j].$ht[$j1].uc($ht[$k]);
			  if ($msar=MSA_check2tails($outfile.$extstr[0],1)) {
			    log_print(1,1,"-$outfile$extstr[0] exists");
			    next;
			  }
			  if ($k==1 && ($msar=MSA_check2tails($outfile.$extstr[1],1))) {
			    log_print(1,1,"-$outfile$extstr[1] exists");
          str_print($outfile.$extstr[0],$msar);
			    next;
			  }
				$msar=$align_prof->($pfiles[0][$j].$extstr[$k],$pfiles[1][$j1].$extstr[$k],$tfiles[0],$tfiles[1],$tfiles[2],$outfile.$extstr[$k]);
        str_print($outfile.$extstr[0],$msar) if ($k==1);
			} #for k=0..1 hot
		} #for j1=0..1
	} #for j=0..0/1
	my_cmd("rm -f prof* tr* in* *$extstr[1]") if ($debug<3) ;
	#-----------------------------------------
} #for i=0..nbr-1 splits
#-----------------------------------------
# done CoS
#=========================================
cleanup(0);
exit(0);

################################################
# main end ------------------------------------#
################################################


#################################################
##---  mafft functions -------------------------#
#################################################


sub met_init_MAF {
	if ($seqtype==2) {
  		print "ERROR: Codon models not supported by MAFFT, use NT instead.\n";
		exit;
	}
	$metver=$msa_program_path.'-profile';
	my $rc=my_cmd("which $metver 2>&1");
	if ($rc=~/which: no/) {
  		print "ERROR: Could not find $metver, please make sure that $metver is on your path and try again.\n";
		exit;
	}
	$metver=$msa_program_path;
	$rc=my_cmd("which $metver 2>&1");
	if ($rc=~/which: no/) {
  		print "ERROR: Could not find $metver, please make sure that $metver is on your path and try again.\n";
		exit;
	}

	#	$metver=$metver.' --localpair --maxiterate 1000' if (substr($met,2,1) ne '0');
	#	$metver=$metver.' --localpair --maxiterate 1000' if ($met!~/[01]$/);

	$metver=$metver.' '.$param if (substr($met,2,1) eq 'T');
	$metver=$metver.' --localpair --maxiterate 1000' if (substr($met,2,1) eq 'M');

	$metver=$metver.($seqtype==0 ? " --amino" : " --nuc");
	$metver=$metver." --quiet" if ($debug<2) ;
	log_print(1,0,"metvar: ".$metver);
	return;
}#sub met_init_MAF


################################################
#----------------------------------------------#
################################################


sub make_guide_tree_MAF {
	my ($infile,$treefile)=@_;
	my $cmdstr="($metver --treeout $infile  | tr  '[:lower:]'  '[:upper:]' | sed 's/>SEQ/>seq/' >hot_H$extstr[0] 2>&1 ) 2>&1 ";
	my $rc=my_cmd($cmdstr);# produces input.fasta.tree
	if ($rc=~/err/i) {
    log_print(0,2,"ERROR: $my_name $$ : mafft error:\n$cmdstr\n---\n$rc\n---\n");
	  cleanup(1);
	}
	my @rtime= $rc=~m/\nreal ([0-9\.]+).*\nuser ([0-9\.]+).*\nsys ([0-9\.]+)/;
	print TIMEFILE "$treefile ".join(',',@rtime)."\n";
	MSA_check2tails("hot_H$extstr[0]",0);
	#$cmdstr="sed 's/_*//g;s/[0-9]*s/s/g;s/\$/;/' $infile.tree >$treefile;rm $infile.tree"; v2.05: mafft changed treefile format
	$cmdstr="sed 's/_*//g;s/[0-9]*s/s/g' $infile.tree >$treefile;rm $infile.tree";

	$rc=my_cmd($cmdstr);# makes guide_tree.nwk
	if ($rc=~/err/i) {
    log_print(0,2,"ERROR: $my_name $$ : mafft guide tree sed error:\n$cmdstr\n---\n$rc\n---\n");
	  cleanup(1);
	}
	return;
}#sub make_guide_tree_MAF


################################################
#----------------------------------------------#
################################################


sub align_seq_MAF {
	my ($infile,$treefile,$outfile)=@_;
	my $cmdstr="($metver --treein $treefile $infile  | tr  '[:lower:]'  '[:upper:]' | sed 's/>SEQ/>seq/' >$outfile 2>&1 ) 2>&1";
	my $rc=my_cmd($cmdstr);
	if ($rc=~/err/i) {
    log_print(0,2,"ERROR: $my_name $$ : mafft error:\n$cmdstr\n---\n$rc\n---\n");
	  cleanup(1);
	}
	my @rtime= $rc=~m/\nreal ([0-9\.]+).*\nuser ([0-9\.]+).*\nsys ([0-9\.]+)/;
	print TIMEFILE "$outfile ".join(',',@rtime)."\n";
	return MSA_check2tails($outfile,0);
}#sub align_seq_MAF


################################################
#----------------------------------------------#
################################################


sub align_prof_MAF {
	my ($pfile1,$pfile2,$tfile1,$tfile2,$tfile3,$ofile)=@_;
#	my $cmdstr = ($tfile1 eq '') ? $pfile1 : "--seed $pfile1 /dev/null"  ;# terminal/internal branch
#	$cmdstr="(time -p $metver --treein $tfile3 --seed $pfile2 $cmdstr | sed 's/>_seed_/>/'  | tr  [:lower:]  [:upper:] | sed 's/>SEQ/>seq/'  >$ofile 2>&1 )2>&1";
	# switched order on purpose! - to accomodate seed vs. seed and seed vs. seq and mafft trees
#	$cmdstr=~s/--localpair// if (substr($met,2,1) eq 'F'); #if both localpair and treein, mafft USED TO says: "Both structure and user tree have been given. Not yet supported!"
#	$cmdstr="(time -p mafft-profile $pfile2 $pfile1 2>&1 >$ofile  )2>&1" if (substr($met,2,1) eq '1'); 
#	print "$met \n$cmdstr\n";
#	exit;
	$cmdstr="($msa_program_path-profile $pfile2 $pfile1  | tr  '[:lower:]'  '[:upper:]' | sed 's/>SEQ/>seq/' 2>&1 >$ofile  )2>&1"; 
	my $rc=my_cmd($cmdstr);
	if ($rc=~/err/i) {
		log_print(0,2,"ERROR: $my_name $$ : mafft error:\n$cmdstr\n---\n$rc\n---\n");
		cleanup(1);
	}
	my @rtime= $rc=~m/\nreal ([0-9\.]+).*\nuser ([0-9\.]+).*\nsys ([0-9\.]+)/;
	print TIMEFILE "$ofile ".join(',',@rtime)."\n";
	return MSA_check2tails($ofile,0);
}#sub align_prof_MAF

################################################
# end mafft functions  ------------------------#
################################################

#################################################
##---  prank functions -------------------------#
#################################################


sub met_init_PRK {
	#$metver='prank';
	$metver=$msa_program_path;
	my $rc=my_cmd("which $metver 2>&1");
	if ($rc=~/which: no/) {
  	print "ERROR: Could not find $metver, please make sure that $metver is on your path and try again.\n";
	  exit;
	}
	$rc=my_cmd("$metver 2>&1");
	# added in ver 2.04: prank interface change
	$prkver= $rc=~/showtree/ ? 1 : 0 ; 
	$metver=$metver.' '.$param; 
	$metver=$metver." -noxml -nopost" if ($prkver==0);
	$metver=$metver." -quiet" if ($debug<2) ;
	$metver=$metver." -codon" if ($seqtype==2) ;
	log_print(1,0,"metver: ".$metver."\nprank version is ".$prkver);
	print("metver: ".$metver."\nprank version is ".$prkver);
	return;
}#sub met_init_PRK


################################################
#----------------------------------------------#
################################################


sub make_guide_tree_PRK {
	my ($infile,$treefile)=@_;
	my $cmdstr="($metver -d=$infile -o=prank_gdt 2>&1) 2>&1;mv prank_gdt.2.fas hot_H$extstr[0];mv prank_gdt.2.dnd $treefile";
	if ($prkver==1) {
		$cmdstr="($metver -d=$infile -o=prank_gdt -showtree 2>&1) 2>&1;mv prank_gdt.2.fas hot_H$extstr[0];mv prank_gdt.2.dnd $treefile";
	} 
	my $rc=my_cmd($cmdstr);# produces input.fasta.tree
	if ($rc=~/err/i) {
    log_print(0,2,"ERROR: $my_name $$ : prank error:\n$cmdstr\n---\n$rc\n---\n");
	  cleanup(1);
	}
	my @rtime= $rc=~m/\nreal ([0-9\.]+).*\nuser ([0-9\.]+).*\nsys ([0-9\.]+)/;
	print TIMEFILE "$treefile ".join(',',@rtime)."\n";
	MSA_check2tails("hot_H$extstr[0]",0);
	my_cmd("rm -f prank_*") if ($debug<3) ;
	return;
}#sub make_guide_tree_PRK


################################################
#----------------------------------------------#
################################################


sub align_seq_PRK {
	my ($infile,$treefile,$outfile)=@_;
	my $cmdstr="($metver -d=$infile -t=$treefile -o=prank_$outfile -notree 2>&1) 2>&1;mv prank_$outfile.1.fas $outfile"; # HERE was 2.fas
	if ($prkver==1) {
		$cmdstr="($metver -d=$infile -t=$treefile -o=prank_$outfile 2>&1) 2>&1;mv prank_$outfile.2.fas $outfile";
	}
	my $rc=my_cmd($cmdstr);
	if ($rc=~/err/i) {
    log_print(0,2,"ERROR: $my_name $$ : prank error:\n$cmdstr\n---\n$rc\n---\n");
	  cleanup(1);
	}
	my @rtime= $rc=~m/\nreal ([0-9\.]+).*\nuser ([0-9\.]+).*\nsys ([0-9\.]+)/;
	print TIMEFILE "$outfile ".join(',',@rtime)."\n";
	my $rmsa=MSA_check2tails($outfile,0);
	my_cmd("rm -f prank_*") if ($debug<3) ;
	return $rmsa;
}#sub align_seq_PRK


################################################
#----------------------------------------------#
################################################


sub align_prof_PRK {
	my ($pfile1,$pfile2,$tfile1,$tfile2,$tfile3,$ofile)=@_;
	$cmdstr="sed 's/^>.*\$/& group_a/' $pfile1 >prank_$ofile\_inp;sed 's/^>.*\$/& group_b/' $pfile2 >>prank_$ofile\_inp";
	$rc=my_cmd($cmdstr);# makes prank profile input
	#exit;
	if ($rc=~/err/i) {
    log_print(0,2,"ERROR: $my_name $$ : prank input sed error:\n$cmdstr\n---\n$rc\n---\n");
	  cleanup(1);
	}
	my $cmdstr="($metver  -partaligned -d=prank_$ofile\_inp -t=$tfile3 -o=prank_$ofile -notree 2>&1) 2>&1;mv prank_$ofile.0.fas $ofile";
	if ($prkver==1) {
		$cmdstr="($metver  -partaligned -d=prank_$ofile\_inp -t=$tfile3 -o=prank_$ofile 2>&1) 2>&1;mv prank_$ofile.0.fas $ofile";
	}
	my $rc=my_cmd($cmdstr);
	if ($rc=~/err/i) {
    log_print(0,2,"ERROR: $my_name $$ : prank error:\n$cmdstr\n---\n$rc\n---\n");
	  cleanup(1);
	}
	my @rtime= $rc=~m/\nreal ([0-9\.]+).*\nuser ([0-9\.]+).*\nsys ([0-9\.]+)/;
	print TIMEFILE "$ofile ".join(',',@rtime)."\n";
	my $rmsa=MSA_check2tails($ofile,0);
	my_cmd("rm -f prank_*") if ($debug<3) ;
	return $rmsa;
}#sub align_prof_PRK

################################################
# end prank functions  ------------------------#
################################################

################################################
#  ClustalW2 functions ------------------------#
################################################

sub met_init_CW2 {
	if ($seqtype==2) {
  		print "ERROR: Codon models not supported by CLUSTALW2, use NT instead.\n";
		exit;
	}
	#$metver='clustalw2';
	$metver=$msa_program_path;
	my $rc=my_cmd("which $metver 2>&1");
	if ($rc=~/which: no/) {
  	print "ERROR: Could not find $metver, please make sure that $metver is on your path and try again.\n";
	  exit;
	}
	$metver=$metver.' '.$param  if (substr($met,2,1) eq 'W');
	$metver=$metver.' -iteration=alignment' if (substr($met,2,1) eq '2') ;
	$metver=$metver.' -iteration=tree' if (substr($met,2,1) eq '3') ;
	$metver=$metver.($seqtype==0 ? " -type=PROTEIN" : " -type=DNA");
	$metver=$metver." -quiet" if ($debug<2) ;
	log_print(1,0,"metvar: ".$metver);
	return;
}#sub met_init_CW2

################################################
#----------------------------------------------#
################################################

sub make_guide_tree_CW2 {
	my ($infile,$treefile)=@_;
	$cmdstr="($metver -infile=$infile -newtree=$treefile 2>&1) 2>&1";
	my $rc=my_cmd($cmdstr);
	if ($rc=~/err/i) {
    log_print(0,2,"ERROR: $my_name $$ : clustalw2 error:\n$cmdstr\n---\n$rc\n---\n");
	  cleanup(1);
	}
	my @rtime= $rc=~m/\nreal ([0-9\.]+).*\nuser ([0-9\.]+).*\nsys ([0-9\.]+)/;
	print TIMEFILE "$treefile ".join(',',@rtime)."\n";
	return;
}#sub make_guide_tree_CW2

################################################
#----------------------------------------------#
################################################

sub align_seq_CW2 {
	my ($infile,$treefile,$outfile)=@_;
	$cmdstr="($metver -infile=$infile -outfile=$outfile -output=fasta -outorder=input -usetree=$treefile 2>&1) 2>&1";
	my $rc=my_cmd($cmdstr);
	if ($rc=~/err/i) {
    log_print(0,2,"ERROR: $my_name $$ : clustalw2 error:\n$cmdstr\n---\n$rc\n---\n");
	  cleanup(1);
	}
	my @rtime= $rc=~m/\nreal ([0-9\.]+).*\nuser ([0-9\.]+).*\nsys ([0-9\.]+)/;
	print TIMEFILE "$outfile ".join(',',@rtime)."\n";
	return MSA_check2tails($outfile,0);
}#sub align_seq_CW2

################################################
#----------------------------------------------#
################################################

sub align_prof_CW2 {
	my ($pfile1,$pfile2,$tfile1,$tfile2,$tfile3,$ofile)=@_;
	my $cmdstr=($tfile1 eq '') ? '' : "-usetree1=$tfile1";# terminal/internal branch
	$cmdstr="($metver -profile -profile1=$pfile1 -profile2=$pfile2 $cmdstr -usetree2=$tfile2 -output=fasta -outfile=$ofile 2>&1) 2>&1";
	my $rc=my_cmd($cmdstr);
	if ($rc=~/err/i) {
    log_print(0,2,"ERROR: $my_name $$ : clustalw2 error:\n$cmdstr\n---\n$rc\n---\n");
	  cleanup(1);
	}
	my @rtime= $rc=~m/\nreal ([0-9\.]+).*\nuser ([0-9\.]+).*\nsys ([0-9\.]+)/;
	print TIMEFILE "$ofile ".join(',',@rtime)."\n";
	return MSA_check2tails($ofile,0);
}#sub align_prof_CW2


#################################################
##   end clustalw2 functions -------------------#
#################################################

sub MSA_check2tails {
  # Read fasta MSA file and check that sequences match global %hot_seqs
	# Reverse residue order of sequences and return in string
	# MSA_check2tails($input_file,$mode)
	# mode=[1:return tails or shutdown on error,  | 
	#       0:just check, dont shutdown and return boolean, this is for resuming aborted runs]
	#
	my ($file,$mode)=@_;
	if (! -e $file) {
	  return 0 if $mode;
    log_print(0,2,"ERROR: $my_name $$ : File not found: $_[0]");
	  cleanup(1);
	}
	open(INFILE,$file);
	local $/='>';
	$_=<INFILE>;
	if ($_ ne '>') {
		close(INFILE);
	  return 0 if $mode;
    log_print(0,2,"ERROR: $my_name $$ : File not in fasta format: $_[0]");
	  cleanup(1);
	}
	my $sdir= ($file=~m/$extstr[0]$/) ? 0 : 1; #seq dirction from file extention
	my $fastar='';
	while ($_=<INFILE>) {
		my ($name,$seq) = $_=~m/^([^\n\r]*)([^>]*)/s;
		$seq=~s/\s//g;
	  my $seqr=reverse($seq);
    $seqr=~s/.{1,60}/$&\n/g;
	  $fastar=$fastar.">$name\n$seqr";
		my ($sid)= $name=~m/(\d{4})$/;
		$seq=~s/-//g;
		log_print(6,0,"$name $sdir $sid:\n".$seq."\n ref:\n".$hot_seqs{seqs}[$sdir][$sid]."\n-----\n");
		if ($seq ne $hot_seqs{seqs}[$sdir][$sid]) {
  		close(INFILE);
  	  return 0 if $mode;
      log_print(0,2,"ERROR: $my_name $$ : Sequence mismatch in file: $_[0] $name :\n---\n$seq\nshould be:\n$hot_seqs{seqs}[$sdir][$sid]\n---\n");
	    cleanup(1);		  
		}
	}
  close(INFILE);
	return $fastar;
}#sub MSA_check2tails
#################################################
##----------------------------------------------#
#################################################
sub fasta2hotseqs {
	# Read fasta file, write sequence names file, fill sequence structure %hot_seqs
	#
	if (! -e $_[0]) {
    log_print(0,2,"ERROR: $my_name $$ : File not found: $_[0]");
	  cleanup(1);
	}
	open(INFILE,$_[0]);
	local $/=">";
	$_=<INFILE>;
	if ($_ ne '>') {
		close(INFILE);
    log_print(0,2,"ERROR: $my_name $$ : File not in fasta format: $_[0]");
	  cleanup(1);
	}
	my $names_txt='';
	my $i=0;
	while ($_=<INFILE>) {
		# catch names containing '>'
		while (($_!~/\n>$/) && (my $tmp=<INFILE>)) {
			$_=$_.$tmp;
		}
		my $sn=sprintf('seq%04u',$i);
		push @{$hot_seqs{sid}},$sn;
		my ($name,$seq) = $_=~m/^([^\n\r]*)([^>]*)/s;
		push @{$hot_seqs{ids}},$name;
		$names_txt=$names_txt.$sn.' '.$name."\n";
		$seq=~s/[\s-]//g;
		$seq=uc($seq);
		push @{$hot_seqs{seqs}[0]},$seq;
		my $seq1=$seq;
		$seq1=~s/.{1,60}/$&\n/g;
		push @{$hot_seqs{fasta}[0]},">$sn\n$seq1";
		$seq=reverse($seq);
		push @{$hot_seqs{seqs}[1]},$seq;
		$seq=~s/.{1,60}/$&\n/g;
		push @{$hot_seqs{fasta}[1]},">$sn\n$seq";
		$i++;
	}#while INFILE
	close(INFILE);
	$hot_seqs{notu}=$i;
	open(OUTFILE,">seq_names.txt");
	print OUTFILE $names_txt;
	close(OUTFILE);
	return;
}#sub fasta2hotseqs
#################################################
##----------------------------------------------#
#################################################
sub tree2split {
	if (! -e $_[0]) {
    log_print(0,2,"ERROR: $my_name $$ : File not found: $_[0]");
	  cleanup(1);
	}
	my ($nwstr,$nbr,$notu,@otu,@otu2,@bid,@len,@br2,%otus);
	open(INFILE,$_[0]);
#	local $/=';'; v2.05: mafft changed treefile format
	local $/;
	$nwstr=<INFILE>;
	close(INFILE);
	$nwstr=~s/;*$/;/g;
	$nwstr=~s/\s//g;
	
	$subtrees{tree}=$nwstr;
	chop $nwstr;
	if ($hot_seqs{notu}==2) {
	  if ($metcmds{$met} eq "MAF") {
	    $nwstr=~m/.*:([.\d]*),.*:([.\d]*).*/;
	    $subtrees{tree}="1 2 $1 $2\n";
	  }#if MAF
	  $subtrees{nbr}=0;
	  return;
  }#notu==2
	$nwstr=~tr/()/<>/;
	$nbr=0;
	log_print(1,0,"nwstr=\n$nwstr");
	################## terminal branches
	while ($nwstr=~m/[<,]([^,:<>]*):([^,:<>]*)[>,]/) {
		$br2[$nbr]=[-1,-1];
		$otu[$nbr]=$1;
		$otus{$1}=$nbr;
		$len[$nbr]=$2;
		$bid[$nbr][0]=$nbr;
		$nwstr=~s/$1:$2/$nbr/;
		log_print(6,0,"tbrn nbr=$nbr\nnwstr=\n$nwstr\n---\n");
		$nbr++;
	}
	$notu=$nbr;
	############### internal branches
	while ($nwstr=~m/<(\d*),(\d*)>:([^,:<>]*)/) {
		$br2[$nbr]=[$1,$2];
		$otu[$nbr]=sprintf("<%s:%f,%s:%f>",$otu[$1],$len[$1],$otu[$2],$len[$2]);
		#$otu[$nbr]="<".$otu[$1].":".$len[$1].",".$otu[$2].":".$len[$2].">";
		$len[$nbr]=$3;
		$bid[$nbr]=[sort(@{$bid[$1]},@{$bid[$2]})];
		$nwstr=~s/<$1,$2>:$3/$nbr/;
		log_print(6,0,"ibrn nbr=$nbr\nnwstr=\n$nwstr\n---\n");
		$nbr++;
	}
	############## last 3
	my @z=$nwstr=~/,/g; #check if rooted- ',' / unrooted- ',,'
	if ($#z==0) {#change rooted to unrooted
		$nwstr=~m/<(\d*),(\d*)>/;
		$nbr--;
		my $i= ($1==$nbr) ? $2 : $1 ;
		$len[$i]=$len[$1]+$len[$2];
		$nwstr=~s/$nbr/$br2[$nbr][0],$br2[$nbr][1]/;
		log_print(6,0,"unroot nbr=$nbr\nnwstr=\n$nwstr\n---\n");
	}# if rooted
	$nwstr=~m/<(\d*),(\d*),(\d*)>/;
	$otu2[$1]=sprintf("<%s:%f,%s:%f>",$otu[$2],$len[$2],$otu[$3],$len[$3]);
	$otu2[$2]=sprintf("<%s:%f,%s:%f>",$otu[$3],$len[$3],$otu[$1],$len[$1]);
	$otu2[$3]=sprintf("<%s:%f,%s:%f>",$otu[$1],$len[$1],$otu[$2],$len[$2]);
#	$otu2[$1]="<".$otu[$2].":".$len[$2].",".$otu[$3].":".$len[$3].">";
#	$otu2[$2]="<".$otu[$3].":".$len[$3].",".$otu[$1].":".$len[$1].">";
#	$otu2[$3]="<".$otu[$1].":".$len[$1].",".$otu[$2].":".$len[$2].">";
	############## retrace splits
	for my $i (reverse(0..$nbr-1)) {
		if ($br2[$i][0]>-1) {
			for my $j (0..1) {

				$otu2[$br2[$i][$j]]=sprintf("<%s:%f,%s:%f>",$otu2[$i],$len[$i],$otu[$br2[$i][1-$j]],$len[$br2[$i][1-$j]]);
				#$otu2[$br2[$i][$j]]="<".$otu2[$i].":".$len[$i].",".$otu[$br2[$i][1-$j]].":".$len[$br2[$i][1-$j]].">";
			}
		}
	}
	############# recode otus
	$subtrees{notu}=$notu;
	$subtrees{otus}=[sort keys %otus];
	$subtrees{nbr}=$nbr;
	$subtrees{len}=[@len];
	my @i2i;
	my $oid=0;
	for my $oids (@{$subtrees{otus}}) {
		$i2i[$otus{$oids}]=$oid;
		$oid++;
	}

	############ otu ids, complement, and structure fill
	my @splits_txt;
	for my $i (0..$nbr-1) {
		my @a0=@{$bid[$i]};
		for my $j (0..$#a0) {
			$a0[$j]=$i2i[$a0[$j]];
		}
		@a0= sort { $a <=> $b } @a0;
		my @a1=(0..$notu-1);
		for my $j (reverse(@a0)) {
			splice(@a1,$j,1);
		}#for j
		if ($#a0+$#a1+2 != $notu) {
			log_print(0,2,"ERROR: subtrees error:\n  ".($#a0+1)." : ".join(",",@a0)."\n  ".($#a1+1)." : ".join(",",@a1));
	  	cleanup(1);
		}#if
		$subtrees{br}[$i][0]{notu}=$#a0+1;
		$subtrees{br}[$i][1]{notu}=$#a1+1;
		$subtrees{br}[$i][2]{notu}=$notu;
		$subtrees{br}[$i][0]{otu}=[@a0];
		$subtrees{br}[$i][1]{otu}=[@a1];
		$subtrees{br}[$i][2]{otu}=[@a1,@a0];# switched order on purpose, for mafft trees

		my $nwstr0=$otu[$i];
		$nwstr0=~tr/<>/()/;
		$subtrees{br}[$i][0]{tree}=$nwstr0.";";
		$nwstr=$otu2[$i];
		$nwstr=~tr/<>/()/;
		$subtrees{br}[$i][1]{tree}=$nwstr.";";
		if ($metcmds{$met} eq "MAF") {
			$otu2[$3]=sprintf("<%s:%f,%s:%f>",$otu[$1],$len[$1],$otu[$2],$len[$2]);
			$subtrees{br}[$i][2]{tree}=sprintf("(%s:%f,%s:%f);",$nwstr0,($len[$i]/2),$nwstr,($len[$i]/2));
		}
		else {
			$subtrees{br}[$i][2]{tree}=sprintf("(%s:%f,%s;",$nwstr0,$len[$i],substr($nwstr,1));
		}
		$subtrees{br}[$i][0]{name}=sprintf('b%u#%04u',(($i<$notu)?1:0),$i);

		my $split_disp="$subtrees{br}[$i][0]{name} : $subtrees{br}[$i][0]{notu}/$subtrees{br}[$i][1]{notu} : [".join(",",@{$subtrees{br}[$i][0]{otu}})."]/[".join(",",@{$subtrees{br}[$i][1]{otu}})."]\n";
		push @splits_txt,$split_disp;
		if ($debug>1) {
			$split_disp="---- tree2split:\n".$split_disp."  left:    $subtrees{br}[$i][0]{tree}\n".
					"  right:    $subtrees{br}[$i][1]{tree}\n"."  joined:    $subtrees{br}[$i][2]{tree}\n-------\n";
			log_print(6,0,$split_disp);
		}#if
	}# for i
	open(OUTFILE,">splits.txt");
	print OUTFILE join('',@splits_txt);
	close(OUTFILE);
	$subtrees{br}[$subtrees{nbr}][0]{tree}=$subtrees{tree};
	$subtrees{br}[$subtrees{nbr}][0]{otu}=[0..$subtrees{notu}-1];
	$subtrees{br}[$subtrees{nbr}][0]{notu}=$subtrees{notu};
	$subtrees{br}[$subtrees{nbr}][1]{notu}=0;
	$subtrees{br}[$subtrees{nbr}][2]{notu}=0;
	if ($metcmds{$met} eq "MAF") {
		newick2mafft();
	}#if
}#sub tree2split
#################################################
##----------------------------------------------#
#################################################

sub newick2mafft {
	for my $i (0..$subtrees{nbr}) {
		for my $j (0..2) {
			next if ($subtrees{br}[$i][$j]{notu}<2);
			log_print(6,0,"-br $i,$j : $subtrees{br}[$i][$j]{tree}\n");
			log_print(6,0,"-br $i,$j : ".join(",",@{$subtrees{br}[$i][$j]{otu}})."\n");
		  my $tr=$subtrees{br}[$i][$j]{tree};
		  $tr=~tr/()/<>/;
		  my $k;
		  for $k (1..$subtrees{br}[$i][$j]{notu}) {
		  	$tr=~s/$subtrees{'otus'}[$subtrees{br}[$i][$j]{otu}[$k-1]]/$k/;
		  }
		  my @mtr;
		  while ($tr=~/<(\d+):([\d\.]+),(\d+):([\d\.]+)>/) {
		  	if ($1<$3) {
		  		$k=$1;
		  		push @mtr,sprintf("%5d%5d%11.5f%11.5f",$1,$3,$2,$4);
		  	}
		  	else {
		  		$k=$3;
		  		push @mtr,sprintf("%5d%5d%11.5f%11.5f",$3,$1,$4,$2);
		  	}
		  	$tr=~s/<$1:$2,$3:$4>/$k/;
		  }
		  log_print(6,0,"---\n$i $j : $subtrees{br}[$i][$j]{notu} $#mtr\n$subtrees{br}[$i][$j]{tree} \n".join(' ; ',@mtr)."\n---\n");
		 	$subtrees{br}[$i][$j]{tree}=join("\n",@mtr)."\n"; 
		 	log_print(6,0,"---\n$subtrees{br}[$i][$j]{tree}---\n");
		}#j
	}#i
	log_print(3,1,"- $subtrees{tree}\n----\n");
	$subtrees{tree}=$subtrees{br}[$subtrees{nbr}][0]{tree};
	log_print(3,1,"-\n$subtrees{tree}----\n");
#	exit;
	return;
}#sub newick2mafft

################################################
# end of logic functions  ---------------------#
################################################


################################################
# misc. supporting functions    ---------------#
################################################



sub str_print {
  my $file=shift(@_);
  open(OUTFILE,">$file");
  print OUTFILE @_;
  close(OUTFILE);
  return;
}#sub str_print
################################################
#----------------------------------------------#
################################################
sub log_print { # debug_level, whereto 0:2=LOG|OUT|ERR , msgstr
	my ($level,$to,$msgstr)=@_;
	return if ($level>$debug);
	chomp($msgstr);
	$msgstr=$msgstr."\n";
	my ($package, $filename, $line) = caller;
	$msgstr="@ line $line of $filename\n$msgstr\@---\n" 
	    if ($msgstr!~m/^-/);
	print LOGFILE $msgstr if (defined(fileno LOGFILE));
	print STDOUT $msgstr if ($to>0);
	print STDERR $msgstr if ($to>1); 
	return;
}#sub log_print
################################################
#----------------------------------------------#
################################################
sub my_cmd {
	my $cmdstr=$_[0];
	my ($package, $filename, $line) = caller;
	log_print(0,$debug,"---- sh: line $line of $filename\n$cmdstr\n");
	my $rc=`$cmdstr`;
	log_print(1,$debug-1,"---- output:\n$rc\n----\n");
	return $rc;
}#sub my_cmd
################################################
#----------------------------------------------#
################################################
sub cleanup { #1st argument : exit state [0:normal ,1:error] 2nd: msgstr
	$SIG{TERM}='IGNORE';
	my $state=$_[0];
	close(LOGFILE);
 	print TIMEFILE "Total: ".join(',',times,time-$t0)."\n".localtime()."\n";
  close(TIMEFILE);
	if ($state==0) {
  	my_cmd("rm -f prof* tr* in* *$extstr[1] temp* pre") if ($debug<4) ;
  	chdir("..");
		my_cmd("mkdir -p $output_dir");
		if ($tgz) {
			my_cmd("tar -czf $done_tar_file $output_prefix; rm -rf $output_prefix");
			log_print(0,1,"---\n $output_prefix $my_name done : dir saved to $done_tar_file\n".localtime()."\n");
		}
		else {
			my_cmd("rm -rf $output_dir/$output_prefix; mv -f $output_prefix $output_dir") if($basedir!~/$output_dir/);
			log_print(0,1,"---\n $output_prefix $my_name done : dir saved to $output_dir/$output_prefix\n".localtime()."\n");
		}
	  exit;
	} 
	else {
  	chdir("..");
		my_cmd("mkdir -p $output_dir\_err");
		if ($tgz) {
			my_cmd("tar -czf $err_tar_file $output_prefix");
			log_print(0,1,"$output_prefix $my_name error : tmp dir saved to $err_tar_file\n".localtime()."\n");
		}
		else {
			my_cmd("cp -rf $output_prefix $output_dir\_err");
			log_print(0,1,"$output_prefix $my_name error : tmp dir saved to $output_dir\_err/$output_prefix\n".localtime()."\n");
		}
	}
 	my_cmd("rm -rf $output_prefix") if ($debug<4);
	exit;
}#sub cleanup
################################################
#----------------------------------------------#
################################################
sub my_sigtrap {	# 1st argument is signal name
	my ($sig) = @_;
	log_print(0,0,"$my_name $$ caught a SIG$sig -- ".(localtime)."\n");
	wait;
	exit(0);
}# sub my_sigtrap
####################################################
####################################################
###   THE END
####################################################
####################################################
# FFU:
#
#$metcmds{MUS}='(time -p muscle -in input.fasta -out output.fasta 2>&1) 2>&1';
#$metcmds{PRC}='(time -p probcons input.fasta  2>&1 >output.fasta) 2>&1';
#$metcmds{DIA}='(time -p dialign-t '.$ENV{HOME}.'/app/dialign/ input.fasta output.fasta 2>&1) 2>&1';
#$metcmds{PCM}='(time -p pcma input.fasta 2>&1;clustalw -infile=input.aln -convert -output=fasta -outfile=output.fasta 2>&1) 2>&1';
#$metcmds{POA}='(time -p poa -read_fasta input.fasta -toupper -clustal input.aln '.$ENV{HOME}.'/app/poa/blosum80.mat 2>&1) 2>&1;clustalw -infile=input.aln -convert -output=fasta -outfile=output.fasta 2>&1';
#$metcmds{MCO}='(time -p t_coffee input.fasta -special_mode mcoffee -run_name output -output fasta_aln, score_ascii -quiet stdout 2>&1) 2>&1;mv output.fasta_aln output.fasta 2>&1';
#$metcmds{TCO}='(time -p t_coffee input.fasta -run_name output -output fasta_aln, score_ascii -quiet stdout 2>&1) 2>&1;mv output.fasta_aln output.fasta 2>&1';
