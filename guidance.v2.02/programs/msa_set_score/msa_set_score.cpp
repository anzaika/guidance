//=======================================================================//
//=======================================================================//
// FILE: msa_set_score.cpp
//=======================================================================//
//=======================================================================//

#include <cstdlib>
#include "defs.h"

//----------------------------------------------------------------------//
//Globals
int debug_level=0;

//=======================================================================//
//=======================================================================//
// main: msa_set_score
//=======================================================================//
//=======================================================================//

int main(int argc, char *argv[]) { //msa_set_score
  
  string my_name=argv[0];
  
  char *env_debug_level=getenv("DEBUG_LEVEL");
  (env_debug_level==NULL) ? debug_level=0 : sscanf(env_debug_level,"%d",&debug_level);

  string usage="\nmsa_set_score v2.01\n\nUsage:\n\
  \n    msa_set_score  file_with_MSA_to_score   output_files_prefix   -m file_with_a_single_alternative_MSA                \n  Or\
  \n    msa_set_score  file_with_MSA_to_score   output_files_prefix   -d directory_with_multiple_alternative_MSA_files     \n  Or\
  \n    msa_set_score  file_with_MSA_to_score   output_files_prefix   -f file_with_a_list_of_multiple_alternative_MSA_files\n\
  \nAll MSA files must be in FASTA format\n\n";

/*
 string usage="\nmsa_set_score v2.0\n\nUsage:\
  \n    msa_set_score  -m file_with_a_single_alternative_MSA                  file_with_MSA_to_score  output_files_prefix\n  Or\
  \n    msa_set_score  -d directory_with_multiple_alternative_MSA_files       file_with_MSA_to_score  output_files_prefix\n  Or\
  \n    msa_set_score  -f file_with_a_list_of_multiple_alternative_MSA_files  file_with_MSA_to_score  output_files_prefix\n\
  \nAll MSA files must be in FASTA format\n\n";
*/
  if(argc!=5) {
    cout << usage;
    exit(0);
  }//if


//---------------------------------------

  // v2.0 string alt_source_flag=argv[1];
  string alt_source_flag=argv[3];
  if(alt_source_flag!="-d" && alt_source_flag!="-f" && alt_source_flag!="-m") {
    cout << usage+"\n\nIllegal option: "+alt_source_flag+"\n\n";
    exit(0);
  }//if

//---------------------------------------
/* v2.0:
  string alt_msa_source=argv[2];
  string ref_msa_file=argv[3];    
  string output_prefix=argv[4];
  
  if (debug_level>0) cout <<  my_name+" "+alt_source_flag+" "+alt_msa_source+" "+ref_msa_file+" "+output_prefix+"\n" ;

*/  
  
  string alt_msa_source=argv[4];
  string ref_msa_file=argv[1];    
  string output_prefix=argv[2];
  
  if (debug_level>0) cout <<  my_name+" "+ref_msa_file+" "+output_prefix+" "+alt_source_flag+" "+alt_msa_source+"\n" ;

//---------------------------------------
//do it

  t_cnt *counter=set_count(alt_source_flag,alt_msa_source,ref_msa_file);
  print_scores(counter,output_prefix) ;
  delete_counters(counter);
  exit(0);
}//main msa_set_score
  


//=======================================================================//
//=======================================================================//
// END main: msa_set_score
//=======================================================================//
//=======================================================================//



//=======================================================================//
//=======================================================================//
// END FILE: msa_set_score.cpp
//=======================================================================//
//=======================================================================//

