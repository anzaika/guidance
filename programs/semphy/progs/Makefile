#! /usr/local/bin/gmake
# makfile for yaep5
Libsources= 
# do not use: fromInstructionFile.cpp

bootstrap_prog:  bootstrap_prog_cmdline.o

.SECONDARY: bootstrap_prog_cmdline.h  bootstrap_prog_cmdline.c

LibCsources= 
# LibCsources += getopt.c getopt1.c 

EXEC = bootstrap_prog
TEST_EXEC_SUB = 
TEST_EXEC = $(addprefix tests/,$(TEST_EXEC_SUB))
LIB = 

#CC=g++
CXX=g++
CC=$(CXX)

#GENGETOPT=/cs/++/phd/ninio/gengetopt-2.7/src/gengetopt
GENGETOPT=gengetopt

libDir=../lib

#LDFLAGS=  
BASICCPPFLAGS = -Wall -Wno-sign-compare -I. -I$(libDir) -I.. -DLOG -ftemplate-depth-25
CPPFLAGS=  -O3 $(BASICCPPFLAGS)
CPPFLAGS+= -I/usr/include/g++-v3
#CPPFLAGS+= -DLOG -DLOGCLS -DMEMCHK

# sources
sources= $(Libsources) $(LibCsources) $(addsuffix .cpp,$(EXEC) $(TEST_EXEC))


.PHONY: tests lib test dat debug %.debug

all: $(EXEC)

test: tests
	+cd tests; make -k test

debug: CPPFLAGS = -g $(BASICCPPFLAGS)
debug: all

%.debug: CPPFLAGS = -g $(BASICCPPFLAGS)
%.debug: %
	@echo "made \""$(*)"\" in debug mode"


lib: $(LIB)

$(LIB): $(Libsources:.cpp=.o) $(LibCsources:.c=.o)
	ar rv $@ $?
	ranlib $@

tags: *.cpp *.h
	etags --members --language=c++ $^ 

$(EXEC) $(TEST_EXEC): ../libSEMPHY.a $(libDir)/libEvolTree.a
tests: $(TEST_EXEC)

include make.dep

install:
	cd ../fast; make -f Makefile.lib install_do


clean:
	-rm -f $(LIB) $(EXEC) $(TEST_EXEC) *.o

depend makedep: _make.dep
	@mv -f _make.dep make.dep 

_make.dep: $(sources)
	@echo making depend
#	$(SHELL) -ec '$(CC) -MM $(CPPFLAGS) $^ | sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@ ; [ -s $@ ] || rm -f $@' 
	@$(SHELL) -ec '$(CC) -MM $(CPPFLAGS) $^ > $@' 
_fast:
	cd ../fast;make -f Makefile.lib -k all

fast.% _fast.%:
	cd ../fast;make -f Makefile.lib -k $(*)


# commandline (gengetopts)
%_cmdline.h %_cmdline.c: %.ggo
	$(GENGETOPT) $(GENGETOPTFLAGS) -i$< -F$(*)_cmdline -f $(*)_cmdline_parser  --arg-struct-name=$(*)_args_info 

#%.h %.c: %.ggo
#	$(GENGETOPT)  --file-name=$(*) -i$<


%.dat.q: %.dat
	awk 'BEGIN{RS="[\n\r]+";};{print "\" "$$0" \"\r"}' $< >  $@
#	cat $@

DAT = cpREV45.dat.q dayhoff.dat.q jones.dat.q mtREV24.dat.q wag.dat.q

dat: $(DAT)

cleandat:
	rm $(DAT)

datMatrixHolder.o: $(DAT)

.PRECIOUS: $(DAT)

# DO NOT DELETE
