libDir = ../../../lib
CPPFLAGS= -I/usr/include/g++-v3/i386-linux -I/usr/include/g++-v3
CPPFLAGS+=  -O3 -Wall -Wno-sign-compare -I.. -I$(libDir) -DLOG -ftemplate-depth-25
LDFLAGS= -L.. -L$(libDir)
LDLIBS= -lEvolTree -lSEMPHY

CXX=g++
CC=$(CXX)

.PHONY: test tests run debug
.SECONDARY: semphy.%.params semphy.%.text
all:test
test tests: $(TESTS)

SEMPHY := ../semphy
constests: $(CONSTS)

CTESTS = bootstrap_test
CPROGS = constraints_test

SHCONSTESTS =  cons_test1 cons_test2 cons_test3 cons_test2b cons_test3b
SEMPHYTESTNAMES = lys_test1 lys_test2 lys_test3 lys_test4 lys_test5 lys_test6 lys_test7 lys_test8 lys_test9

SHTESTS = $(SHCONSTESTS)

SEMPHYTESTS = $(addprefix semphy.,$(SEMPHYTESTNAMES))

TESTS= $(CTESTS) $(SEMPHYTESTS) $(SHTESTS) 

TMPOUT= $(addsuffix .out.tmp,$(TESTS)) 

run:  $(TMPOUT)



%.debug: SEMPHY = ../semphy.debug
debug: CPPFLAGS += -g -Wall -Wno-sign-compare 
debug: 
	@echo working in debug mode - with SEMPHY=$(SEMPHY)
%.debug: debug %
	@echo "made \""$(*)"\" in debug mode"



$(CTESTS) $(CPROGS): $(libDir)/libEvolTree.a ../libSEMPHY.a

$(addsuffix .out.tmp,$(SHTESTS)): $(SEMPHY)
test: $(addsuffix .test,$(TESTS))

%.test: %.out.tmp %.out.standard
	diff -I "^#" $^

semphy.%.out.tmp: semphy.%.params semphy.%.text $(SEMPHY)
	cat $(word 2,$^)  > $@ && $(word 3,$^) `cat $<`  >> $@ 

cons_test1.out.tmp: $(CPROGS)

cons_test1.out.tmp: lys.tree lys.tree.consraint2 lys.tree.consraint3

%.out.tmp: %
	$(*) >& $@

tests: $(TESTS)


semphy.%.text: %
	grep "echo" $< | awk 'BEGIN{FS="\""};{print $$2}' > $@

semphy.%.params: %
	grep "../semphy" $< | awk 'BEGIN{FS="../semphy"};{print $$2}' > $@

semphy.%.out.standard: %.out.standard
	cp $< $@

#.INTERMEDIATE: $(TMPOUT)


cleantmp:
	-rm -f  $(TMPOUT)


clean: cleantmp
	-rm -f $(CTEST) 
