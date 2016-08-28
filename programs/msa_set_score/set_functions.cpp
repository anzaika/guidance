//=======================================================================//
//=======================================================================//
// FILE: set_functions.cpp
//=======================================================================//
//=======================================================================//


#include <cstdlib>
#include "defs.h"


//=======================================================================//
//=======================================================================//
// function: t_cnt *set_count(string alt_source_flag,string alt_msa_source,string ref_msa_file)
//=======================================================================//
//=======================================================================//

t_cnt *set_count(string alt_source_flag, string alt_msa_source, string ref_msa_file) {
  
	//--------------------------------
	// major variables:

  vector<string> alt_file_list;
  t_msa *ref,*alt;
  t_cnt *counter;


	//--------------------------------
	// create alt file name list

  struct stat list_stat;
  if ( stat(alt_msa_source.c_str(),&list_stat) ) {
    cerr << "Error finding alternative MSAs file/directory: "+alt_msa_source+"\n";
    exit(0);
  }
  if(alt_source_flag=="-m") {	// single fasta file
    if (!S_ISREG(list_stat.st_mode)) {
      cerr << "Error finding alternative MSA file: "+alt_msa_source+"\n";
      exit(0);
    }
  	alt_file_list.push_back(alt_msa_source);
  }// flag -m
  else if (alt_source_flag=="-f") {  // file list
    if (!S_ISREG(list_stat.st_mode)) { 
      cerr << "Error finding alternative MSAs list file: "+alt_msa_source+"\n";
      exit(0);
    }
    ifstream alt_list_file;
    alt_list_file.open(alt_msa_source.c_str());
    if(alt_list_file.fail()){
      cerr << "Error opening alternative MSAs list file: "+alt_msa_source+"\n";
      exit(0);
    }
    string alt_filename;
    while(!alt_list_file.eof()){
      getline(alt_list_file,alt_filename);
        if(!alt_filename.empty()){
          alt_file_list.push_back(alt_filename);
        }//if
    }//while
  }// flag -f
  else if (alt_source_flag=="-d") {		// whole directory
    if (!S_ISDIR(list_stat.st_mode)) {
      cerr << "Error finding alternative MSAs directory: "+alt_msa_source+"\n";
      exit(0);
    }
    DIR *dir;
    struct dirent *dent;
    if(!(dir = opendir(alt_msa_source.c_str()))) {
        cerr << "Error opening alternative MSAs directory: "+alt_msa_source+"\n";
        exit(0) ;
    }    
    if (alt_msa_source[alt_msa_source.length()-1]!='/') {
      alt_msa_source.append("/");
    }
    while(dent=readdir(dir)) {
      stat((alt_msa_source+dent->d_name).c_str(),&list_stat);
      if (S_ISREG(list_stat.st_mode)) {
        alt_file_list.push_back(alt_msa_source+dent->d_name);
      }
    }//while readdir
    closedir(dir);
  }//flag -d
  
  if (debug_level>0) {
    cout << "alt_file_list.size = " <<alt_file_list.size()<<"\n";
    for(unsigned int i=0; i<alt_file_list.size();i++){
      cout << alt_file_list[i] << "\n";
    }//
  } // if debug
	//=======================================================================//


	//--------------------------------
	//Load Reference file
	
  ref=read_fasta(ref_msa_file);
  msa_recode(ref,0);
  counter=init_counters(ref);

	//=======================================================================//
  
  counter->alt_source=alt_source_flag+" "+alt_msa_source;
	
	//----------------------------------------------------------//
	//for each alt file: load and count
  for(unsigned int i=0; i<alt_file_list.size();i++){
  	alt=read_fasta(alt_file_list[i]);
    msa_recode(alt,ref);
    add_msa(counter,alt);
    delete_msa(alt);
  } //for alt
	//=======================================================================//

	//----------------------------------------------------------//
	//done
  // note: counter still points to ref, so don't delete it!!!
  return(counter);

}// function: t_cnt *set_count(string alt_source_flag,string alt_msa_source,string ref_msa_file)

//=======================================================================//
//=======================================================================//
// END function: t_cnt *set_count(string alt_source_flag,string alt_msa_source,string ref_msa_file)
//=======================================================================//
//=======================================================================//





//=======================================================================//
//=======================================================================//
// function: void print_scores(t_cnt *counter,string output_prefix)
//=======================================================================//
//=======================================================================//

void print_scores(t_cnt *counter,string output_prefix) {

  int row,row2;
  long col,count;
  long double score,hit;
  string out_name;
  FILE *out;
  
  //----------------------------------------------------------//
  //_msa.scr
  
  out_name=output_prefix+"_msa.scr";
  out=fopen(out_name.c_str(),"w");
  if(out==NULL){
    cerr<<"Error opening output file: "+out_name+"\n";
    exit(0);
  }//if
  
  fprintf(out,"#REF_FILE %s\n",counter->ref->file_name.c_str());
  fprintf(out,"#ROWS %4d  #COLUMNS %6ld\n",counter->row,counter->col);
  fprintf(out,"#ALT_FILES %s \n",counter->alt_source.c_str());
  fprintf(out,"#N_ALT %6d\n",counter->nalt);
  
  //mean residue pair
  count=0;
  hit=0;
  for(col=0;col<counter->col;col++) {
    for(row=0;row<counter->row-1;row++) {
      for(row2=row+1;row2<counter->row;row2++)  {
        if(counter->res_pair_hit[row][row2][col]>-1) {
          count++;
          hit+=double(counter->res_pair_hit[row][row2][col])/counter->nalt;
        }//if res_pair_hit
      }//for row2    
    }//for row
  }//for col
  score=hit/count;
  fprintf(out,"#MEAN_RES_PAIR_SCORE %8.6Lf  ",score);
  
  //mean column
  hit=0;
  for(col=0;col<counter->col;col++) {
    hit+=double(counter->col_hit[col])/counter->nalt;
  }
  score=hit/counter->col;
  fprintf(out,"#MEAN_COL_SCORE %8.6Lf\n",score);
  fprintf(out,"#END\n");
  fclose(out); 
  //end _msa.scr

  //--------------------------------------------------//
  //_col_col
  
  out_name=output_prefix+"_col_col.scr";
  out=fopen(out_name.c_str(),"w");
  if(out==NULL){
    cerr<<"Error opening output file: "+out_name+"\n";
    exit(0);
  }
  
  fprintf(out,"#COL_NUMBER  #COL_SCORE\n");
  for(col=0;col<counter->col;col++) {
    score=double(counter->col_hit[col])/counter->nalt;
    fprintf(out,"%6ld %9.6Lf\n",col+1,score);
  }
  fprintf(out,"#END\n");
  fclose(out); 

  //-------------------------------------------------//
  //res_pair_col  
  
  out_name=output_prefix+"_res_pair_col.scr";
  out=fopen(out_name.c_str(),"w");
  if(out==NULL){
    cerr<<"Error opening output file: "+out_name+"\n";
    exit(0);
  }
  
  fprintf(out,"#COL_NUMBER  #RES_PAIR_COLUMN_SCORE\n");
  for(col=0;col<counter->col;col++){
    hit=0;count=0;
    for(row=0;row<counter->row-1;row++){
      for(row2=row+1;row2<counter->row;row2++){       
        if(counter->res_pair_hit[row][row2][col]>-1){ 
          count++;
          hit+=double(counter->res_pair_hit[row][row2][col])/counter->nalt;
        }
      } //for row2
    } //for row
    if(count>0) {
      score=hit/count;
      fprintf(out,"%6ld %9.6Lf\n",col+1,score);
    }
  }//for col
  fprintf(out,"#END\n");
  fclose(out); 

  //---------------------------------------------------//
  //residue pair sequence
  
  out_name=output_prefix+"_res_pair_seq.scr";
  out=fopen(out_name.c_str(),"w");
  if(out==NULL){
    cerr<<"Error opening output file: "+out_name+"\n";
    exit(0);
  }
  
  fprintf(out,"#ROW_NUMBER #RES_PAIR_SEQUENCE_SCORE\n");
  for(row=0;row<counter->row;row++){
    hit=0;count=0;
    for(row2=0;row2<counter->row;row2++){
      if(row!=row2){
        for(col=0;col<counter->col;col++){
          if(counter->res_pair_hit[row][row2][col]>-1){
            count++;  
            hit+=double(counter->res_pair_hit[row][row2][col])/counter->nalt;
           }
         }  //for col
       }  //if
    } //for row2
    score=hit/count;
    fprintf(out,"%4d %9.6Lf\n",row+1,score);
  }//for row
  fprintf(out,"#END\n");
  fclose(out); 

  //---------------------------------------------------//
  // residue pair sequence pair
  
  out_name=output_prefix+"_res_pair_seq_pair.scr";
  out=fopen(out_name.c_str(),"w");
  if(out==NULL){
    cerr<<"Error opening output file: "+out_name+"\n";
    exit(0);
  }
  
  fprintf(out,"#ROW_NUMBER_1  #ROW_NUMBER_2  #RES_PAIR_SEQ_PAIR_SCORE\n");
  for(row=0;row<counter->row-1;row++) {
    for(row2=row+1;row2<counter->row;row2++) { 
      hit=0;count=0;
      for(col=0;col<counter->col;col++){   
        if(counter->res_pair_hit[row][row2][col]>-1){
          count++;
          hit+=double(counter->res_pair_hit[row][row2][col])/counter->nalt;
        }
      }//for col
    score=hit/count;
    fprintf(out,"%4d %4d %9.6Lf\n",row+1,row2+1,score);
    }//for row2
  }//for row
  fprintf(out,"#END\n");
  fclose(out); 

  //---------------------------------------------------//
  // residue pair residue
  out_name=output_prefix+"_res_pair_res.scr";
  out=fopen(out_name.c_str(),"w");
  if(out==NULL){
    cerr<<"Error opening output file: "+out_name+"\n";
    exit(0);
  }
  
  fprintf(out,"#COL_NUMBER  #ROW_NUMBER  #RES_PAIR_RESIDUE_SCORE\n");
  for(col=0;col<counter->col;col++){  
    for(row=0;row<counter->row;row++){  
      if (counter->ref->col2res[row][col]>-1){ 
        hit=0;count=0;
        for(row2=0;row2<counter->row;row2++){
          if ( (row!=row2) && (counter->res_pair_hit[row][row2][col]>-1) ){ 
            count++;
            hit+=double(counter->res_pair_hit[row][row2][col])/counter->nalt;
          }
        }//for row2
        score=hit/count;
        fprintf(out,"%6ld %4d %9.6Lf\n",col+1,row+1,score);
      }//if col2res>-1
    }//for row
  }//for col
  fprintf(out,"#END\n");
  fclose(out); 

  //---------------------------------------------------//
  // residue pair 
  
  out_name=output_prefix+"_res_pair.scr";
  out=fopen(out_name.c_str(),"w");
  if(out==NULL){
    cerr<<"Error opening output file: "+out_name+"\n";
    exit(0);
  }
  
  fprintf(out,"#COL_NUMBER  #ROW_NUMBER_1  #ROW_NUMBER_2  #RES_PAIR_SCORE\n");
  for(col=0;col<counter->col;col++){
    for(row=0;row<counter->row-1;row++)  {
      for(row2=row+1;row2<counter->row;row2++)  {
        if(counter->res_pair_hit[row][row2][col]!=-1){
          score=double(counter->res_pair_hit[row][row2][col])/counter->nalt;
          fprintf(out,"%6ld %4d %4d %9.6Lf\n",col+1,row+1,row2+1,score);
        }
      }//for row2
    }//for row
  }//for col    
  fprintf(out,"#END\n");
  fclose(out); 

//----------------------------------------------------------// 
//done

}//function: void print_scores(t_cnt *counter,string output_prefix)



//=======================================================================//
//=======================================================================//
// function: void print_scores(t_cnt *counter,string output_prefix)
//=======================================================================//
//=======================================================================//




//=======================================================================//
//=======================================================================//
// END FILE: set_functions.cpp
//=======================================================================//
//=======================================================================//
