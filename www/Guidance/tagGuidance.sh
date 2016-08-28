#!/bin/csh
if ($#argv == 0) then
	echo "USAGE: tagGuidance.sh <version number>"
else
	set v = $argv[1];
	svn mkdir http://ibis.tau.ac.il/pupkoSVN/tags/guidance.v$v -m "version $v";
#	svn add tags/guidance.v$v;
	svn cp http://ibis.tau.ac.il/pupkoSVN/trunk/Makefile http://ibis.tau.ac.il/pupkoSVN/tags/guidance.v$v/ -m "version $v";
	svn mkdir http://ibis.tau.ac.il/pupkoSVN/tags/guidance.v$v/libs -m "version $v";
#	svn add tags/guidance.v$v/libs;
	svn cp http://ibis.tau.ac.il/pupkoSVN/trunk/libs/Makefile http://ibis.tau.ac.il/pupkoSVN/tags/guidance.v$v/libs/ -m "version $v";
	svn cp http://ibis.tau.ac.il/pupkoSVN/trunk/libs/phylogeny http://ibis.tau.ac.il/pupkoSVN/tags/guidance.v$v/libs/ -m "version $v";
	svn mkdir http://ibis.tau.ac.il/pupkoSVN/tags/guidance.v$v/programs -m "version $v";
#	svn add tags/guidance.v$v/programs;
	svn cp http://ibis.tau.ac.il/pupkoSVN/trunk/programs/Makefile http://ibis.tau.ac.il/pupkoSVN/tags/guidance.v$v/programs/ -m "version $v";
	svn cp http://ibis.tau.ac.il/pupkoSVN/trunk/programs/Makefile.generic http://ibis.tau.ac.il/pupkoSVN/tags/guidance.v$v/programs/ -m "version $v";
	svn cp http://ibis.tau.ac.il/pupkoSVN/trunk/programs/semphy http://ibis.tau.ac.il/pupkoSVN/tags/guidance.v$v/programs/ -m "version $v";
	svn cp http://ibis.tau.ac.il/pupkoSVN/trunk/programs/msa_set_score http://ibis.tau.ac.il/pupkoSVN/tags/guidance.v$v/programs/ -m "version $v";
	svn cp http://ibis.tau.ac.il/pupkoSVN/trunk/programs/isEqualTree http://ibis.tau.ac.il/pupkoSVN/tags/guidance.v$v/programs/ -m "version $v";
	svn cp http://ibis.tau.ac.il/pupkoSVN/trunk/programs/removeTaxa http://ibis.tau.ac.il/pupkoSVN/tags/guidance.v$v/programs/ -m "version $v";
	svn mkdir http://ibis.tau.ac.il/pupkoSVN/tags/guidance.v$v/www -m "version $v";
#	svn add tags/guidance.v$v/www;
	svn cp http://ibis.tau.ac.il/pupkoSVN/trunk/www/bioSequence_scripts_and_constants http://ibis.tau.ac.il/pupkoSVN/tags/guidance.v$v/www/ -m "version $v";
	svn cp http://ibis.tau.ac.il/pupkoSVN/trunk/www/Guidance http://ibis.tau.ac.il/pupkoSVN/tags/guidance.v$v/www/ -m "version $v";
	svn cp http://ibis.tau.ac.il/pupkoSVN/trunk/www/Selecton http://ibis.tau.ac.il/pupkoSVN/tags/guidance.v$v/www/ -m "version $v";
	svn mv http://ibis.tau.ac.il/pupkoSVN/tags/guidance.v$v/www/Guidance/README http://ibis.tau.ac.il/pupkoSVN/tags/guidance.v$v/ -m "version $v";
endif

# Don't forget to remove irrelevant libs and programs from the makefiles 
