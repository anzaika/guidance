#!/usr/bin/perl -w
#************************************************************************************************
# -p : pdb file name
# -w : url to directory with files
# -d : path to address where chimerax files should be created
# -x : name of chimerax file
# -a : alignment file
# -t : tree file
# files to create:
# -c : name of file to hold the color script for the alignment
# -h : name of file to hold the header script for the alignemnt



my @ARGS = @ARGV;  ## This is so later we can re-parse the command line args later if we need to
my $numargv = @ARGS;
my $counter = 0;

my ($pdb_file_name, $www_dir, $files_dir, $chimerax_file, $msa_file, $tree_file, $scf_file, $hdr_file);
my $add_align = "no";
my $add_tree = "no";
my $add_header = "no";

for ($counter = 0; $counter < $numargv; $counter++) {
	if ($ARGS[$counter] =~ /^-p$/i) {                  ## PDB file ##
		$counter++;
		if ($ARGS[$counter] && $ARGS[$counter] !~ /^-/) { $pdb_file_name = $ARGS[$counter]; }
		else { return("WARNING : The argument after -p was not a pdb file!"); $counter--; }
	}
	elsif($ARGS[$counter] =~ /^-w$/i){
		$counter++;
		if ($ARGS[$counter] && $ARGS[$counter] !~ /^-/) { $www_dir = $ARGS[$counter]; }
		else { return("WARNING : The argument after -w was not a url!"); $counter--; }
	}
	elsif($ARGS[$counter] =~ /^-d$/i){
		$counter++;
		if ($ARGS[$counter] && $ARGS[$counter] !~ /^-/) { $files_dir = $ARGS[$counter]; }
		else { return("WARNING : The argument after -d was not a directory!"); $counter--; }
	}
	elsif($ARGS[$counter] =~ /^-x$/i){
		$counter++;
		if ($ARGS[$counter] && $ARGS[$counter] !~ /^-/) { $chimerax_file = $ARGS[$counter]; }
		else { return("WARNING : The argument after -x was not a chimerx file!"); $counter--; }
	}		
	elsif($ARGS[$counter] =~ /^-a$/i){
		$counter++;
		if ($ARGS[$counter] && $ARGS[$counter] !~ /^-/) {
			$msa_file = $ARGS[$counter];
			$add_align = "yes";
		}
		else { return("WARNING : The argument after -a was not MSA!"); $counter--; }
	}
	elsif($ARGS[$counter] =~ /^-t$/i){
		$counter++;
		if ($ARGS[$counter] && $ARGS[$counter] !~ /^-/) {
			$tree_file = $ARGS[$counter];
			$add_tree = "yes";
		}
		else { return("WARNING : The argument after -t was not a tree file!"); $counter--; }
	}
	elsif($ARGS[$counter] =~ /^-c$/i){
		$counter++;
		if ($ARGS[$counter] && $ARGS[$counter] !~ /^-/) { $scf_file = $ARGS[$counter]; }
		else { return("WARNING : The argument after -c was not a scf file!"); $counter--; }
	}
	elsif($ARGS[$counter] =~ /^-h$/i){
		$counter++;
		if ($ARGS[$counter] && $ARGS[$counter] !~ /^-/) {
			$hdr_file = $ARGS[$counter];
			$add_header = "yes";
		}
		else { return("WARNING : The argument after -h was not a hdr file!"); $counter--; }
	}
}
unless (open CHIMERAX, ">".$files_dir.$chimerax_file) {die "cannot open ".$files_dir.$chimerax_file." $!";}
print CHIMERAX <<EndOfFile;
<?xml version="1.0"?>
  <ChimeraPuppet type="std_webdata">
<web_files>
<file  name="$pdb_file_name" format="text" loc="$www_dir$pdb_file_name"/>
EndOfFile
if ($add_align eq "yes"){
	print CHIMERAX <<EndOfFile;
<file  name=\"$msa_file\" format=\"text\" loc=\"$www_dir$msa_file\"/>";
</web_files>
<commands>
  <mid_cmd>preset apply pub 3;repr cpk;show;focus</mid_cmd>

<!-- the following 3 lines locate the Multalign Viewer instance
	that was created by opening the alignment file, and stores a reference
	to the instance as the variable 'mav' -->
  <py_cmd>from MultAlignViewer.MAViewer import MAViewer</py_cmd>
  <py_cmd>from chimera.extension import manager</py_cmd>
  <py_cmd>mav = [inst for inst in manager.instances if isinstance(inst, MAViewer)][-1]</py_cmd>

EndOfFile
}
# if there is no MSA, print the relevant colors in order to color the molecule
else{
	print CHIMERAX <<EndOfColors;
</web_files>
<commands>
  <mid_cmd>colordef CONS10 1.00 1.00 0.59</mid_cmd>
  <mid_cmd>colordef CONS9 0.63 0.15 0.38</mid_cmd>
  <mid_cmd>colordef CONS8 0.94 0.49 0.67</mid_cmd>
  <mid_cmd>colordef CONS7 0.98 0.79 0.87</mid_cmd>
  <mid_cmd>colordef CONS6 0.99 0.93 0.96</mid_cmd>
  <mid_cmd>colordef CONS5 1.00 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS4 0.92 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS3 0.84 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS2 0.55 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS1 0.06 0.78 0.82</mid_cmd>
  <mid_cmd>color CONS10 @/bfactor=0</mid_cmd>
  <mid_cmd>color CONS9 @/bfactor=9 </mid_cmd>
  <mid_cmd>color CONS8 @/bfactor=8</mid_cmd>
  <mid_cmd>color CONS7 @/bfactor=7</mid_cmd>
  <mid_cmd>color CONS6 @/bfactor=6</mid_cmd>
  <mid_cmd>color CONS5 @/bfactor=5</mid_cmd>
  <mid_cmd>color CONS4 @/bfactor=4</mid_cmd>
  <mid_cmd>color CONS3 @/bfactor=3</mid_cmd>
  <mid_cmd>color CONS2 @/bfactor=2</mid_cmd>
  <mid_cmd>color CONS1 @/bfactor=1</mid_cmd>
  <mid_cmd>preset apply pub 3;repr cpk;show;focus</mid_cmd>
EndOfColors
}
# the tree can be shown only if MSA was given
	if ($add_tree eq "yes" and $add_align eq "yes"){
print CHIMERAX <<EndOfTree;
<!-- read in/show Consurf tree -->

<!-- hide initial headers, show phylogeny tree, load the coloring file,
	and make residue letters black.
	This uses two possible sets of calls:  one for the 1.2540 release
	and one for later versions that uses a better API.
	The 'if' condition guarantees that the code will work no
	matter what version the user has -->
<py_cmd>
if hasattr(mav, 'loadScfFile'):
	mav.hideHeaders(mav.headers(shownOnly=True))
	mav.usePhylogenyFile("http://consurf.tau.ac.il/results/1217406203/TheTree.txt", askReorder=False)
	mav.loadScfFile("$www_dir$scf_file")
	mav.useColoringFile(None)
else:
	mav.hideHeaders(mav.seqCanvas.headerDisplayOrder())
	mav.usePhylogenyFile("http://consurf.tau.ac.il/results/1217406203/TheTree.txt")
	mav.regionBrowser.loadScfFile("$www_dir$scf_file")
	from MultAlignViewer.prefs import RC_BLACK
	mav.seqCanvas.setColorFunc(RC_BLACK)
</py_cmd>
EndOfTree
}
# header can be added only if there is an alignment to show
if ($add_header eq "yes" and $add_align eq "yes"){
	print CHIMERAX 
"\n<!-- read in/show Consurf headers -->\n  <py_cmd>mav.readHeaderFile(\"$www_dir$hdr_file\")</py_cmd>";
}
print CHIMERAX "\n  </commands>\n</ChimeraPuppet>";
close CHIMERAX;
