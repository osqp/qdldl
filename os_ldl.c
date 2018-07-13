#include "os_ldl.h"

#define UNKNOWN -1
#define USED 1
#define UNUSED 0

/* Compute the elimination tree for a quasidefinite matrix
   in compressed sparse column form.
*/

os_ldl_int osqp_ldl_etree(CscMatrix *A,
                     os_ldl_int* work,
                     os_ldl_int* Lnz,
                     os_ldl_int* etree){

  const os_ldl_int n = A->n;
  os_ldl_int sumLnz = 0;
  os_ldl_int i,j,p;

  // zero out Lnz and work.  Set all etree values to unknown
  for(i = 0; i < n; i++){
    work[i]  = 0;
    Lnz[i]   = 0;
    etree[i] = UNKNOWN;
  }

  for(j = 0; j < n; j++){
    work[j] = j;
    for(p = A->p[j]; p < A->p[j+1]; p++){
      i = A->i[p];
      while(work[i] != j){
        if(etree[i] == UNKNOWN){
          etree[i] = j;
        }
        Lnz[i]++;
        work[i] = j;
        i = etree[i];
      }
    }
  }

  //count the total nonzeros
  for(i = 0; i < n; i++) sumLnz += Lnz[i];

  return sumLnz;
}



void osqp_ldl_factor(CscMatrix *A,
                     CscMatrix *L,
                     os_ldl_float* D,
                     os_ldl_float* Dinv,
                     os_ldl_int* Lnz,
                     os_ldl_int* etree,
                     os_ldl_bool*,bwork,
                     os_ldl_int* iwork,
                     os_ldl_float* fwork){

  os_ldl_int i,j,k,nnzY, bidx, cidx, nextIdx, nnzE, tmpIdx;
  os_ldl_int *yMarkers, *yIdx, *elimBuffer, *LNextSpaceInCol;
  os_ldl_float *yVals;

  //partition working memory into pieces
  yMarkers        = bwork;
  yIdx            = iwork;
  elimBuffer      = iwork + A->n;
  LNextSpaceInCol = iwork + A->n*2;
  yVals           = fwork;

  // First element of the diagonal D
  D[0]    = A->x[0];
  Dinv[0] = 1/D[0];

  //Assign basic structure to L.
  L->n = A->n;
  L->m = A->n;
  L->p[0] = 0;
  for(i = 0; i < A->n; i++){
    L->p[i+1] = L->p[i] + Lnz[i];   //cumsum, total at the end
  }

  // set all Yidx to be 'unused' initially
  for(i = 0; i < L->n; i++){
    yMarkers[i]  = UNUSED;
    yVals[i]     = 0.0;
    LNextSpaceInCol[i] = L->p[i];
  }

  //Start from 1 here. The upper LH corner is trivially 0
  //in L b/c are only computing the subdiagonal elements
  for(k = 1; k < A->n; k++){

    //Initialize D[k] as the last element
    //of this column of triu(A)
    D[k] = A->x[A->p[k+1]-1];

    //number of nonzeros in this row of L
    nnzY = 0;  //number of elements in this row

    //This loop determines where nonzeros
    //will go in the kth row of L, but doesn't
    //compute the actual values
    tmpIdx = (A->p[k+1]-1);
    for(i = A->p[k]; i < tmpIdx; i++){

      bidx        = A->i[i];   // we are working on this element of b
      yVals[bidx] = A->x[i];   // initialise y(bidx) = b(bidx)

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
      for(j = L->p[cidx]; j < tmpIdx; j++){
        yVals[L->i[j]] -= L->x[j]*yVals[cidx];
      }

      //Now I have the cidx^th element of y = L\b.
      //so compute the corresponding element of
      //this row of L and put it into the right place
      L->i[tmpIdx] = k;
      L->x[tmpIdx] = yVals[cidx]*Dinv[cidx];

      //D[k] -= yVals[cidx]*yVals[cidx]*Dinv[cidx];
      D[k] -= yVals[cidx]*L->x[tmpIdx];
      LNextSpaceInCol[cidx]++;

      //reset the yvalues and indices back to zero and UNUSED
      //once I'm done with them
      yVals[cidx]     = 0.0;
      yMarkers[cidx]  = UNUSED;

    } //end for i

    //compute the inverse of the diagonal
    Dinv[k]= 1/D[k];

  } //end for k
}
