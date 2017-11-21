//=======================================================================//
//=======================================================================//
// FILE: count_functions.cpp
//=======================================================================//
//=======================================================================//

#include <cstdlib>
#include "defs.h"

//=======================================================================//
//=======================================================================//
// function: t_cnt *init_counters(t_msa *ref)
//=======================================================================//
//=======================================================================//

t_cnt *init_counters(t_msa *ref) {
  if (debug_level>0) cout << "init_counters: "+ref->file_name+"\n" ;
  int row,row2;
  long col;
  t_cnt *counter= new (t_cnt);
  counter->ref=ref;
  counter->row=ref->row;
  counter->col=ref->col;
  counter->nalt=0;
  counter->alt_source="?";

	//---------------------------------------------------------------------//
	//columns counter  

  counter->col_hit=new int[counter->col];

  for(col=0;col<counter->col;col++) {
    counter->col_hit[col]=0;
  }//for col
  
	//---------------------------------------------------------------------//
	//residue pair counter

  counter->res_pair_hit=new int**[counter->row];

  for(row=0;row<counter->row;row++) {
    counter->res_pair_hit[row]=new int*[counter->row];
  }//for row

  
  for (row=0;row<counter->row-1;row++) {
    for (row2=row+1;row2<counter->row;row2++) {
      counter->res_pair_hit[row][row2]=new (nothrow) int[counter->col];
      if (counter->res_pair_hit[row][row2]==0) {
        printf("ERROR:- Out of memory \n");
        exit(0);
      }
      counter->res_pair_hit[row2][row]=counter->res_pair_hit[row][row2];  // note: just the pointer, only one copy of the actual values.
      for (col=0;col<counter->col;col++) {   
        //if ( (ref->col2res[row][col]!=-1) && (ref->col2res[row2][col]!=-1) ) {
		  if ( (ref->col2res[row][col]%2==1) && (ref->col2res[row2][col]%2==1) ) { // non gap Cmatrix
          counter->res_pair_hit[row][row2][col]=0;
        }//if
        else {
          counter->res_pair_hit[row][row2][col]=-1; // undef residue pair
        }//else
		if (debug_level>9) cout<<row<<","<<row2<<","<<col<<":"<<counter->res_pair_hit[row][row2][col]<<endl;
      }//for row
    }//for row2
  }//for row   
	//---------------------------------------------------------------------//
  return(counter);
}//function t_cnt *init_counters(t_msa *ref)

//=======================================================================//
//=======================================================================//
// END function: t_cnt *init_counters(t_msa *ref)
//=======================================================================//
//=======================================================================//


//=======================================================================//
//=======================================================================//
// function: void add_msa(t_cnt *counter, t_msa *alt)
//=======================================================================//
//=======================================================================//

void add_msa(t_cnt *counter, t_msa *alt) {
  if (debug_level>3) cout << " add_msa "+alt->file_name+"\n" ;
  counter->nalt++;
  int row,row2;
  long col,alt_col,res,res_Cmatrix;
  for(col=0;col<counter->col;col++) {

		//---------------------------------------------------------------------//
    //compare cols:
    
    row=counter->ref->col_id[col];
	//res=counter->ref->col2res[row][col];
    res_Cmatrix=counter->ref->col2res[row][col];
	if (res_Cmatrix%2==1){res=0.5*(res_Cmatrix+1)-1;} 	// from the Cmatrix value to the actual residue number
	else {res=-1;} 										// res=-1 for gap
	
    alt_col=alt->res2col[row][res];				//find alt_col to compare to col  
	
    row=0;
    while( (row<counter->row) && (counter->ref->col2res[row][col]==alt->col2res[row][alt_col]) ) {
		if (debug_level>9) cout<<"row:"<<row<<" col:"<<col<<" alt_col:"<<alt_col<<endl;
		row++;
    }//while
    if (row==counter->row) {
      counter->col_hit[col]++;
    }//if

    //----------------------------------------------------------------//
    //compare residue pairs  
    
    for(row=0;row<counter->row-1;row++) {
      //res=counter->ref->col2res[row][col];
		res_Cmatrix=counter->ref->col2res[row][col];
		if (res_Cmatrix%2==1){res=0.5*(res_Cmatrix+1)-1;} // in the Cmatrix the value for each non gap residue is 2*res+1 -> this take it back to the residue number [starting from 0]
		else {res=-1;} 
		if (res>-1) {   // save some needless comparison when a gap
			alt_col=alt->res2col[row][res];   			//find alt_col to compare to col  
			for(row2=row+1;row2<counter->row;row2++) {
				if( (counter->res_pair_hit[row][row2][col]>-1) && (counter->ref->col2res[row2][col]==alt->col2res[row2][alt_col]) )  {
					counter->res_pair_hit[row][row2][col]++;
				}//if
			}//for col2    
		}// if res >-1
    }//for row

  }//for col
  //---------------------------------------------------------------------//
  if (debug_level>3) cout << " add_msa done \n" ;
}// function: void add_msa(t_cnt *counter, t_msa *alt)

//=======================================================================//
//=======================================================================//
// END function: void add_msa(t_cnt *counter, t_msa *alt)
//=======================================================================//
//=======================================================================//



//=======================================================================//
//=======================================================================//
// function: void free_memory_counters(t_cnt *counter)
//=======================================================================//
//=======================================================================//

void delete_counters(t_cnt *counter) {
  int row,row2;
  delete_msa(counter->ref);
  for(row=0;row<counter->row-1;row++) {
    for(row2=row+1;row2<counter->row;row2++) {
      delete[]counter->res_pair_hit[row][row2];
    }//for row2
  }//for row

  for(row=0;row<counter->row;row++) {
    delete[] counter->res_pair_hit[row];
  }//for row
  delete[] counter->res_pair_hit;
  delete counter->col_hit;
  delete counter;
}// function: void free_memory_counters(t_cnt *counter)

//=======================================================================//
//=======================================================================//
// END function: void free_memory_counters(t_cnt *counter)
//=======================================================================//
//=======================================================================//



//=======================================================================//
//=======================================================================//
// END FILE: count_functions.cpp
//=======================================================================//
//=======================================================================//

