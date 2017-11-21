CC     = g++
CFLAGS = -c -Wall -O3
LFLAGS = -g -Wall -O3 

all: msa_set_score

install: msa_set_score

clean : 
	rm -f *.o *.exe

OBJECTS = msa_functions.o count_functions.o set_functions.o 


msa_set_score : $(OBJECTS) defs.h msa_set_score.o
	$(CC) $(LFLAGS) -o $@ $(OBJECTS) $@.o 


%.o : %.cpp defs.h
	$(CC) $(CFLAGS) $<
                