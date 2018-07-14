#include "qdldl.h"

#define UNKNOWN -1
#define USED 1
#define UNUSED 0

/* Compute the elimination tree for a quasidefinite matrix
   in compressed sparse column form.
*/

QDLDL_int QDLDL_etree(const QDLDL_int  n,
                      const QDLDL_int* Ap,
                      const QDLDL_int* Ai,
                      QDLDL_int* work,
                      QDLDL_int* Lnz,
                      QDLDL_int* etree){

  QDLDL_int sumLnz = 0;
  QDLDL_int i,j,p;


  for(i = 0; i < n; i++){
  // zero out Lnz and work.  Set all etree values to unknown
    work[i]  = 0;
    Lnz[i]   = 0;
    etree[i] = UNKNOWN;

    //Abort if A is not triu or does not have
    //antries at every diagonal term.  There should
    //be at least one entry in every column, and the
    //lst one in column should be the diagonal term
    if(Ap[i] == Ap[i+1] || Ai[Ap[i+1]-1] != i ){
      return -1;
    }

  }

  for(j = 0; j < n; j++){
    work[j] = j;
    for(p = Ap[j]; p < Ap[j+1]; p++){
      i = Ai[p];
      while(work[i] != j){
        if(etree[i] == UNKNOWN){
          etree[i] = j;
        }
        Lnz[i]++;         //nonzeros in this column
        work[i] = j;
        i = etree[i];
      }
    }
  }

  //count the total nonzeros
  for(i = 0; i < n; i++) sumLnz += Lnz[i];

  return sumLnz;
}



QDLDL_int QDLDL_factor(const QDLDL_int    n,
                  const QDLDL_int*   Ap,
                  const QDLDL_int*   Ai,
                  const QDLDL_float* Ax,
                  QDLDL_int*   Lp,
                  QDLDL_int*   Li,
                  QDLDL_float* Lx,
                  QDLDL_float* D,
                  QDLDL_float* Dinv,
                  const QDLDL_int* Lnz,
                  const QDLDL_int* etree,
                  QDLDL_bool*  bwork,
                  QDLDL_int*   iwork,
                  QDLDL_float* fwork){

  QDLDL_int i,j,k,nnzY, bidx, cidx, nextIdx, nnzE, tmpIdx;
  QDLDL_int *yIdx, *elimBuffer, *LNextSpaceInCol;
  QDLDL_float *yVals;
  QDLDL_bool  *yMarkers;
  QDLDL_int   positiveValuesInD = 0;

  //partition working memory into pieces
  yMarkers        = bwork;
  yIdx            = iwork;
  elimBuffer      = iwork + n;
  LNextSpaceInCol = iwork + n*2;
  yVals           = fwork;

  //Assign basic structure to L.
  Lp[0] = 0;
  for(i = 0; i < n; i++){
    Lp[i+1] = Lp[i] + Lnz[i];   //cumsum, total at the end
  }

  // set all Yidx to be 'unused' initially
  for(i = 0; i < n; i++){
    yMarkers[i]  = UNUSED;
    yVals[i]     = 0.0;
    LNextSpaceInCol[i] = Lp[i];
  }

  // First element of the diagonal D.
  D[0]    = Ax[0];
  Dinv[0] = 1/D[0];

  //Start from 1 here. The upper LH corner is trivially 0
  //in L b/c we are only computing the subdiagonal elements
  for(k = 1; k < n; k++){

    //Initialize D[k] as the last element
    //of this column of A (assumed triu)
    D[k] = Ax[Ap[k+1]-1];

    //number of nonzeros in this row of L
    nnzY = 0;  //number of elements in this row

    //This loop determines where nonzeros
    //will go in the kth row of L, but doesn't
    //compute the actual values
    tmpIdx = (Ap[k+1]-1);
    for(i = Ap[k]; i < tmpIdx; i++){

      bidx        = Ai[i];   // we are working on this element of b
      yVals[bidx] = Ax[i];   // initialise y(bidx) = b(bidx)

      // use the forward elimination tree to figure
      // out which elements must be eliminated after
      // this element of b
      nextIdx = bidx;

      if(yMarkers[nextIdx] == UNUSED){   //this y term not already visited

        yMarkers[nextIdx] = USED;     //I touched this one
        elimBuffer[0]     = nextIdx;  // It goes at the start of the current list
        nnzE              = 1;         //length of unvisited elimination path from here

        nextIdx = etree[bidx];

        while(nextIdx != UNKNOWN && nextIdx < k){
          if(yMarkers[nextIdx] == USED) break;

          yMarkers[nextIdx] = USED;   //I touched this one
          elimBuffer[nnzE] = nextIdx; //It goes in the current list
          nnzE++;                     //the list is one longer than before
          nextIdx = etree[nextIdx];   //one step further along tree

        } //end while

        // now I put the buffered elimination list into
        // my current ordering in reverse order
        while(nnzE){
          yIdx[nnzY++] = elimBuffer[--nnzE];
        } //end while
      } //end if

    } //end for i

    //This for loop places nonzeros values in the k^th row
    for(i = (nnzY-1); i >=0; i--){

      //which column are we working on?
      cidx = yIdx[i];

      // loop along the elements in this
      // column of L and subtract to solve to y
      tmpIdx = LNextSpaceInCol[cidx];
      for(j = Lp[cidx]; j < tmpIdx; j++){
        yVals[Li[j]] -= Lx[j]*yVals[cidx];
      }

      //Now I have the cidx^th element of y = L\b.
      //so compute the corresponding element of
      //this row of L and put it into the right place
      Li[tmpIdx] = k;
      Lx[tmpIdx] = yVals[cidx]*Dinv[cidx];

      //D[k] -= yVals[cidx]*yVals[cidx]*Dinv[cidx];
      D[k] -= yVals[cidx]*Lx[tmpIdx];
      LNextSpaceInCol[cidx]++;

      //reset the yvalues and indices back to zero and UNUSED
      //once I'm done with them
      yVals[cidx]     = 0.0;
      yMarkers[cidx]  = UNUSED;

    } //end for i

    //Maintain a count of the positive entries
    //in D.  If we hit a zero, we can't factor
    //this matrix, so abort
    if(D[k] == 0.0){return -1;}
    if(D[k]  > 0.0){positiveValuesInD++;}

    //compute the inverse of the diagonal
    Dinv[k]= 1/D[k];

  } //end for k

  return positiveValuesInD;

}

// Solves (L+I)x = b
void QDLDL_Lsolve(const QDLDL_int    n,
                  const QDLDL_int*   Lp,
                  const QDLDL_int*   Li,
                  const QDLDL_float* Lx,
                  QDLDL_float* x){

QDLDL_int i,j;
  for(i = 0; i < n; i++){
      for(j = Lp[i]; j < Lp[i+1]; j++){
          x[Li[j]] -= Lx[j]*x[i];
      }
  }
}

// Solves (L+I)'x = b
void QDLDL_Ltsolve(const QDLDL_int    n,
                   const QDLDL_int*   Lp,
                   const QDLDL_int*   Li,
                   const QDLDL_float* Lx,
                   QDLDL_float* x){

QDLDL_int i,j;
  for(i = n-1; i>=0; i--){
      for(j = Lp[i]; j < Lp[i+1]; j++){
          x[i] -= Lx[j]*x[Li[j]];
      }
  }
}

// Solves Ax = b where A has given LDL factors
void QDLDL_solve(const QDLDL_int       n,
                    const QDLDL_int*   Lp,
                    const QDLDL_int*   Li,
                    const QDLDL_float* Lx,
                    const QDLDL_float* Dinv,
                    QDLDL_float* x){

QDLDL_int i;

QDLDL_Lsolve(n,Lp,Li,Lx,x);
for(i = 0; i < n; i++) x[i] *= Dinv[i];
QDLDL_Ltsolve(n,Lp,Li,Lx,x);

}
