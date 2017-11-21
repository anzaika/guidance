//=======================================================================//
//=======================================================================//
// FILE: msa_functions.cpp
//=======================================================================//
//=======================================================================//
 
#include <cstdlib>
#include "defs.h"
//const int MSA_rows = 3; // used in line 120 // define constatn fro win compile
//=======================================================================//
//=======================================================================//
// function: t_msa *read_fasta(string filename) 
//=======================================================================//
//=======================================================================//

t_msa *read_fasta(string filename) {   
  if (debug_level>3) cout << "read_fasta: "+filename+"\n" ;
//------------------------------------------//
//file check and open
  struct stat file_stat;
  if (stat(filename.c_str(),&file_stat) || !S_ISREG(file_stat.st_mode)){
    cerr << "Error finding file "+filename+"\n"  ;
    exit(0);
  }
  ifstream infile;
  infile.open(filename.c_str(),ifstream::in);
  if(infile.fail()){
    cerr << "Error opening file "+filename+"\n"  ;
    exit(0); 
  }//if
//-----------------------------------------//
  t_msa *msa= new(t_msa);  
  msa->file_name=filename;
  string line,seq;
  size_t found;

//------------------------------------------//
  while(!infile.eof()) {   

    getline(infile,line);

    // erase whitespaces
    found=0;
    while ((found=line.find_first_of("\n\r\t ",found))!=string::npos) {
      line.erase(found,1);
    }

    if (line[0]=='>') {//seq header
      if(msa->seq_names.size()>0){ // not first header
        	msa->alignment.push_back(seq);
        	seq.clear();
        }
      msa->seq_names.push_back(line);
    }//if seq header 
    else { //sequence  
     	seq.append(line);
    } // else sequence
  }//while
  infile.close();
  msa->alignment.push_back(seq); // last sequence
  msa->col=seq.length();
  msa->row=msa->seq_names.size();
  
  // check lengths and degap
  if (debug_level>3) cout << "read_fasta: row: "<<msa->row<<" col: "<<msa->col<<" ; checking sequences\n" ;
  for(int row=0;row<msa->row;row++){  

    seq=msa->alignment[row];

    if (long(seq.length())!=msa->col) {
      cerr << "Error: Not an alignment, sequence length mismatch in file " << filename << "\n"  ;
      exit(0); 
    }//if 

 	  found=0;
    while ((found=seq.find_first_of("-",found))!=string::npos) {
      seq.erase(found,1);
    }// 

 	  msa->sequences.push_back(seq);

  }//for row 
  
  if (debug_level>3) cout << "read_fasta done \n" ;
  return(msa);

}//function: t_msa *read_fasta(string filename)

//=======================================================================//
//=======================================================================//
// END function: t_msa *read_fasta(string filename) 
//=======================================================================//
//=======================================================================//





//=======================================================================//
//=======================================================================//
// function: void msa_recode(t_msa *msa,t_msa *ref) 
//=======================================================================//
//=======================================================================//

/*
The original version of the algorithm coded gaps as -1
and characters according to their position in the non-gapped
sequence (starting from 0). 
The coding included 'col2res' -> the sequnce as numbers, and
'res2col' -> the indicies of the characters in the alignment.
In the current version the C matrix coding was adopted for 'col2res':
characters are represented by odd numbers and gaps by even numbers.
The 'res2col' matrix remains the same.
For more info about the Cmatrix recoding see page 3 on the supplemantary material of:
Satija R, Novak A., Mikls I., Lyngs R., and Hein J. (2009) BigFoot:
Bayesian alignment and phylogenetic footprinting with MCMC,
BMC Evolutionary Biology, 9, 217.
http://www.biomedcentral.com/content/supplementary/1471-2148-9-217-s1.pdf
*/
void msa_recode(t_msa *msa,t_msa *ref) {
  long row,msa_row,col,res;
  int seq_order[msa->row];
//  int seq_order[MSA_rows]; // For win compile use a constant 
//--------------------------------------------//
//check against ref
  if (debug_level>0) cout << "msa_recode: "+msa->file_name+"\n" ;

  if (ref!=0) { // synchronize sequences

    if (msa->row!=ref->row) {
      cerr << "Error: Number of sequences mismatch at file "+msa->file_name+" "<< msa->row<<" should be "<<ref->row << "\n"  ;
      exit(0);
    }// if row mismatch

    for(row=0;row<ref->row;row++){

      msa_row=0;
      while ((msa_row<msa->row) && (msa->seq_names[msa_row]!=ref->seq_names[row])) {
        msa_row++;
      }// while

      if (msa_row==msa->row) {
        cerr << "Error: Sequence name mismatch at file "+msa->file_name+" "+ref->seq_names[row]+" not found\n"  ;
        exit(0);
      }//if 

      if (debug_level>6) cout << "msa_recode: Sequence match at file "+msa->file_name+": "<< msa_row << " " << row << "\n"+msa->sequences[msa_row]+"\n"+ref->sequences[row]+"\n"  ;

      if (msa->sequences[msa_row]!=ref->sequences[row]) {
        cerr << "Error: Sequence mismatch at file "+msa->file_name+": "<< msa_row << " " << row << "\n"+msa->sequences[msa_row]+"\nshould be:\n"+ref->sequences[row]+"\n"  ;
        exit(0);
      }//if 

      seq_order[row]=msa_row;     
    }//for row row      
  }//if ref - sync


//--------------------------------------------//
//ALLOCATING MEMORY
  if (debug_level>3) cout << "msa_recode: allocating mem \n" ;
	
  msa->col2res= new long*[msa->row];
  msa->res2col= new long*[msa->row];

  for(row = 0; row < msa->row; row++){
    msa->col2res[row]= new long[msa->col];
    msa->res2col[row]= new long[msa->col];
  }//for row row
  msa->col_id= new int[msa->col];

//--------------------------------------------//
//--------------------------------------------//
//initializing reciprocal res<->col pointers
  for(row=0;row<msa->row;row++){  

    // reorder
    msa_row = (ref==0) ? row : seq_order[row];
    if (debug_level>9) cout << "msa_recode: seq sync: row=" << row <<" msa_row= " <<msa_row <<"\n" ;

    res=0;
	int LastNonGap=-1; //for gaps Cmatrix
    for(col=0;col<msa->col;col++){    

      // note msa_row vs. row !! here the reordering is taking place

      if((msa->alignment[msa_row][col])=='-') { // gaps 
        //msa->col2res[row][col]=-1; //OLD
		msa->col2res[row][col]=LastNonGap+1; //Cmatrix
		if (debug_level>9) cout<<"col2res["<<row<<"]["<<col<<"]="<<msa->col2res[row][col]<<endl;
      }//if gap
      else { //residue
        msa->res2col[row][res]=col;
        //msa->col2res[row][col]=res;
		msa->col2res[row][col]=2*res+1; //Cmatrix
		LastNonGap=2*res+1; //Cmatrix
        msa->col_id[col]=row;           // any residue can serve as col_id, so the last sequence with a residue will do                       
        if (debug_level>9) cout<<"\tres2col["<<row<<"]["<<res<<"]="<<col<<endl;
		if (debug_level>9) cout<<"col2res["<<row<<"]["<<col<<"]="<<msa->col2res[row][col]<<endl;
		res++;
      }//else residue

    }//for col
    
  }//for row 

  if (debug_level>3) cout << " msa_recode done \n" ;
}//function: void msa_recode(t_msa *msa,t_msa *ref)


//=======================================================================//
//=======================================================================//
// END function: void msa_recode(t_msa *msa,t_msa *ref) 
//=======================================================================//
//=======================================================================//


//=======================================================================//
//=======================================================================//
// function: void delete_msa(t_msa *msa)
//=======================================================================//
//=======================================================================//


void delete_msa(t_msa *msa) {
  if (debug_level>3) cout << " delete_msa\n" ;

  msa->sequences.clear();
  msa->alignment.clear();
  msa->seq_names.clear();
  for(int row=0;row<msa->row;row++) {
    delete[]msa->col2res[row];
    delete[]msa->res2col[row];
  }
  delete[]msa->col2res;
  delete[]msa->res2col;
  delete[]msa->col_id;
  delete msa;

  if (debug_level>3) cout << " delete_msa done\n" ;
}//function: void delete_msa(t_msa *msa)


//=======================================================================//
//=======================================================================//
// END function: void delete_msa(t_msa *msa)
//=======================================================================//
//=======================================================================//



//=======================================================================//
//=======================================================================//
// END FILE: msa_functions.cpp
//=======================================================================//
//=======================================================================//




