//=======================================================================//
//=======================================================================//
// FILE: defs.h
//=======================================================================//
//=======================================================================//

using namespace std;
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <dirent.h>
#include <string>
#include <vector>
#include <math.h>


//----------------------------------------------------------------------//
//Globals

extern int debug_level;



//---------------------------------------------------------------------//
//Structures

typedef struct msa_struct {	
	string file_name;
	vector<string> seq_names;
	vector<string> alignment;
	vector<string> sequences;
	int row;
	long col;
	int  *col_id;
	long **col2res;
	long **res2col;
} t_msa;

typedef struct counters_struct {
  t_msa *ref;
  string alt_source;
	int row;
	long col;
	int nalt;
	int ***res_pair_hit;
	int *col_hit;
} t_cnt;


//----------------------------------------------------------------------//
//Functions

// in msa_functions.cpp:
t_msa *read_fasta(string filename);
void msa_recode(t_msa *msa,t_msa *ref);
void delete_msa(t_msa *msa);

// in count_functions.cpp
t_cnt *init_counters(t_msa *ref);
void add_msa(t_cnt *counter, t_msa *alt);
void delete_counters(t_cnt *counter);

// in set_functions.cpp
t_cnt *set_count(string alt_source_flag,string alt_msa_source,string ref_msa_file);
void print_scores(t_cnt *counter,string output_prefix) ;


//=======================================================================//
//=======================================================================//
// END FILE: defs.h
//=======================================================================//
//=======================================================================//
