#! /usr/local/bin/gmake
# makfile for yaep5
Libsources= StdAfx.cpp cAngles.cpp cCoeffMatrices.cpp cDRotationMatrix.cpp cDataProbModel.cpp cEigenVals.cpp cEigenValsOptimise.cpp  cOccuranceData.cpp cOccurancePair.cpp cOrthonormalMatrix.cpp cProbModel.cpp  cProbs.cpp cRotationMatrix.cpp cRow.cpp cSquareMatrix.cpp cUpdateAngle.cpp readInput.cpp errorMsg.cpp cOptManager.cpp 
#cLineOpt.cpp cProbModelOptimizer.cpp

#ProbModel.cpp StdAfx.cpp cAngles.cpp cCoeffMatrices.cpp cDRotationMatrix.cpp cEigenVals.cpp cLineOpt.cpp cOccuranceData.cpp cOccurancePair.cpp cOrthonormalMatrix.cpp cProbModelOptimizer.cpp cProbs.cpp cRotationMatrix.cpp cRow.cpp cSquareMatrix.cpp readInput.cpp cProbModel.cpp cEigenValsOptimise.cpp



# ProbModel.cpp  StdAfx.cpp cEigenVals.cpp cLineOpt.cpp cOccuranceData.cpp cOccurancePair.cpp cOrthonormalMatrix.cpp cProbModelOptimizer.cpp cProbs.cpp cRow.cpp cSquareMatrix.cpp cRotationMatrix.cpp cAngles.cpp


EXEC = semphyOptManager semphyOptMat
TEST_EXEC_SUB = 
TEST_EXEC= $(addprefix tests/,$(TEST_EXEC_SUB))
LIB = libSEMPHYMatOpt.a

DATA=test3.data pi3.data

all: $(LIB) $(libEvol) $(EXEC) 

semphy: main.o 
generateSeq: generateOptions.o
estimateLikelihood: estimateLikelihoodOptions.o
#STANFORD=yes

# ifdef STANFORD
#   CC=/usr/locallocal/bin/g++
# else
CC=g++
#endif


ifndef libDir
libDir=../lib
endif
ifndef libIncDir 
libIncDir=lib
endif
ifndef libEvol
libEvol=lib/libEvolTree.a
endif

#DWBegin
#ifdef STANFORD
#  CPPFLAGS=  -g -Wall -Wno-sign-compare  -I/user2/danwei/share/gengetopt -I$(libDir) -DLOG
#  CPPFLAGS=  -O3 -Wall -Wno-sign-compare  -I/user2/danwei/share/gengetopt -I. -I$(libDir) -DLOG
#else
#  CPPFLAGS=  -g -Wall -Wno-sign-compare  -I. -I$(libDir) -DLOG
CPPFLAGS=  -O3 -Wall -Wno-sign-compare  -I. -I$(libDir) -DLOG
#endif

#LDFLAGS += -static
#LDFLAGS += -R/usr/locallocal/lib

#CPPFLAGS+= -DLOG -DLOGCLS -DMEMCHK

# sources
sources= $(Libsources) $(addsuffix .cpp,$(EXEC) $(TEST_EXEC))

.PHONY: tests lib test debug %.debug tar

all: $(LIB) $(EXEC) 

test: tests
	cd tests; make -k test

debug: CPPFLAGS = -g -Wall -Wno-sign-compare  -I. -I$(libDir) -DLOG
debug: all

%.debug: CPPFLAGS = -g -Wall -Wno-sign-compare  -I. -I$(libDir) -DLOG
%.debug: %
	@echo "made \""$(*)"\" in debug mode" 




$(libEvol) le:
	cd $(libDir);make -f Makefile.lib all

lib: $(LIB)

$(LIB): $(Libsources:.cpp=.o)
	ar rv $@ $?
	ranlib $@

tags: *.cpp *.h $(libDir)/*.h  $(libDir)/*.cpp
	etags --members --language=c++ $^ 

ifdef STANFORD
  $(EXEC) $(TEST_EXEC): $(LIB) $(libEvol) libGetOpt.a
else
  $(EXEC) $(TEST_EXEC): $(LIB)
endif

tests: $(TEST_EXEC)

include make.dep

install:
	cd fast; make install_do


clean:
	-rm -f $(LIB) $(EXEC) $(TEST_EXEC) *.o

cleanlib:
	cd lib;make clean

depend makedep: _make.dep
	@mv -f _make.dep make.dep 

_make.dep: $(sources)
	@echo making depend
#	$(SHELL) -ec '$(CC) -MM $(CPPFLAGS) $^ | sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@ ; [ -s $@ ] || rm -f $@' 
	@$(SHELL) -ec '$(CC) -MM $(CPPFLAGS) -I. $^ > $@' 
_fast:
	cd fast;make -k all

fast.% _fast.%:
	cd fast;make -k $(*)

tar: code.tar

code.tar: $(sources) *.h Makefile matlab/*.m matlab/*.mat $(DATA)
	tar cvf $@ $^

# DO NOT DELETE
