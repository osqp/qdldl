module QDLDL

using amd

const int QDLDL_UNKNOWN = -1;
const bool QDLDL_USED   = true;
const bool QDLDL_UNUSED = false;

struct QDLDLFactorisation

    #contains factors L, D^{-1}
    L::SparseMatrixCSC
    Dinv
    #permutation vector
    p

    #constructor for user-specified permutation, including no permutation
    QDLDLFactorisation(A::SparseMatrixCSC{Tv,Ti},p::Array{Ti} = Void)
        if p isa Void
            #no permutation at all
            L, Dinv = factor(triu(A);
        else
            #factor the matrix A using permutation p
            L, Dinv = factor(triu(A(p,p)));
        end
        return new(L,Dinv,p)
    end

    #constructor for the default (AMD) permutation
    QDLDLFactorisation(A::SparseMatrixCSC)
        #use AMD ordering as a default
        p = amd(A);
        return QDLDLFactorisation(A,p);
    end

end

# Solves Ax = b using LDL factors for A.
# Returns x, preserving b
function solve(QDLDL::QDLDLFactorisation,b)
    return QDLDL_solve(QDLDL.L,QDLDL.DInv,b);
end

# Solves Ax = b using LDL factors for A.
# Solves in place (x replaces b)
function solve!(QDLDL::QDLDLFactorisation,b)
    return QDLDL_solve!(QDLDL.L,QDLDL.DInv,b);
end


function factor(A::SparseMatrixCSC{Tv,Ti})

    etree  = Array{Ti}(A.n);
    Lnz    = Array{Ti}(A.n);
    iwork  = Array{Ti}(A.n*3);
    bwork  = Array{Bool}(A.n);
    fwork  = Array{Tv}(A.n);

    #compute elimination gree using QDLDL converted code
    sumLnz = QDLDL_etree!(A,iwork,Lnz,etree)

    if(sumLnz < 0)
        error("A matrix is not upper triangular or has an empty column")
    end

    #allocate space for the L matrix row indices and data
    Lp = Array{Ti}(A.n + 1);
    Li = Array{Ti}(sumLnz);
    Lx = Array{Tv}(sumLnz);

    #allocate for and D D inverse
    D  = Array{Tv}(A.n);
    Di = Array{Tv}(A.n);

    #factor using QDLDL converted code
    posDCount = QDLDL_factor(A::SparseMatricCSC,Lp,Li,Lx,D,Dinv,Lnz,etree,bwork,iwork,fwork)

    if(posDCount < 0)
        error("Zero entry in D (matrix is not quasidefinite)")
    end

    L = SparseMatrixCSC(A.n,A.n,Lp,Li,Lx);

    return L, Dinv

end



# Compute the elimination tree for a quasidefinite matrix
# in compressed sparse column form.

function QDLDL_etree!(A::SparseMatricCSC,work,Lnz,etree){

    for i = 1:n
      # zero out Lnz and work.  Set all etree values to unknown
        work[i]  = 0;
        Lnz[i]   = 0;
        etree[i] = QDLDL_UNKNOWN;

        #Abort if A doesn't have at least one entry
        #one entry in every column
        if(Ap[i] == Ap[i+1])
          return -1;
        end
    end

    for j = 1:n
        work[j] = j;
        for p = Ap[j]:Ap[j+1]
          i = Ai[p];
            if(i > j)
                return -1;
            end
            while(work[i] != j)
                if(etree[i] == QDLDL_UNKNOWN)
                    etree[i] = j;
                end
                Lnz[i]++;         #nonzeros in this column
                sumLnz += 1;      # Increase the total number of nonzeros
                work[i] = j;
                i = etree[i];
            end
        end #end for p
    end

  return sumLnz;
end




function QDLDL_factor(A::SparseMatricCSC,Lp,Li,Lx,D,Dinv,Lnz,etree,bwork,iwork,fwork)


  QDLDL_int i,j,k,nnzY, bidx, cidx, nextIdx, nnzE, tmpIdx;
  QDLDL_int *yIdx, *elimBuffer, *LNextSpaceInCol;
  QDLDL_float *yVals;
  QDLDL_float yVals_cidx;
  QDLDL_bool  *yMarkers;
  QDLDL_int   positiveValuesInD = 0;

  #partition working memory into pieces
  yMarkers        = bwork;
  yIdx            = iwork;
  elimBuffer      = iwork + n;
  LNextSpaceInCol = iwork + n*2;
  yVals           = fwork;


  Lp[0] = 0; #first column starts at index zero

  for(i = 0; i < n; i++){

    #compute L column indices
    Lp[i+1] = Lp[i] + Lnz[i];   #cumsum, total at the end

    # set all Yidx to be 'unused' initially
    #in each column of L, the next available space
    #to start is just the first space in the column
    yMarkers[i]  = QDLDL_UNUSED;
    yVals[i]     = 0.0;
    D[i]         = 0.0;
    LNextSpaceInCol[i] = Lp[i];
  }

  # First element of the diagonal D.
  D[0]     = Ax[0];
  if(D[0] == 0.0){return -1;}
  if(D[0]  > 0.0){positiveValuesInD++;}
  Dinv[0] = 1/D[0];

  #Start from 1 here. The upper LH corner is trivially 0
  #in L b/c we are only computing the subdiagonal elements
  for(k = 1; k < n; k++){

    #NB : For each k, we compute a solution to
    #y = L(0:(k-1),0:k-1))\b, where b is the kth
    #column of A that sits above the diagonal.
    #The solution y is then the kth row of L,
    #with an implied '1' at the diagonal entry.

    #number of nonzeros in this row of L
    nnzY = 0;  #number of elements in this row

    #This loop determines where nonzeros
    #will go in the kth row of L, but doesn't
    #compute the actual values
    tmpIdx = Ap[k+1];

    for(i = Ap[k]; i < tmpIdx; i++){

      bidx = Ai[i];   # we are working on this element of b

      #Initialize D[k] as the element of this column
      #corresponding to the diagonal place.  Don't use
      #this element as part of the elimination step
      #that computes the k^th row of L
      if(bidx == k){
        D[k] = Ax[i];
        continue;
      }

      yVals[bidx] = Ax[i];   # initialise y(bidx) = b(bidx)

      # use the forward elimination tree to figure
      # out which elements must be eliminated after
      # this element of b
      nextIdx = bidx;

      if(yMarkers[nextIdx] == QDLDL_UNUSED){   #this y term not already visited

        yMarkers[nextIdx] = QDLDL_USED;     #I touched this one
        elimBuffer[0]     = nextIdx;  # It goes at the start of the current list
        nnzE              = 1;         #length of unvisited elimination path from here

        nextIdx = etree[bidx];

        while(nextIdx != QDLDL_UNKNOWN && nextIdx < k){
          if(yMarkers[nextIdx] == QDLDL_USED) break;

          yMarkers[nextIdx] = QDLDL_USED;   #I touched this one
          elimBuffer[nnzE] = nextIdx; #It goes in the current list
          nnzE++;                     #the list is one longer than before
          nextIdx = etree[nextIdx];   #one step further along tree

        } #end while

        # now I put the buffered elimination list into
        # my current ordering in reverse order
        while(nnzE){
          yIdx[nnzY++] = elimBuffer[--nnzE];
        } #end while
      } #end if

    } #end for i

    #This for loop places nonzeros values in the k^th row
    for(i = (nnzY-1); i >=0; i--){

      #which column are we working on?
      cidx = yIdx[i];

      # loop along the elements in this
      # column of L and subtract to solve to y
      tmpIdx = LNextSpaceInCol[cidx];
      yVals_cidx = yVals[cidx];
      for(j = Lp[cidx]; j < tmpIdx; j++){
        yVals[Li[j]] -= Lx[j]*yVals_cidx;
      }

      #Now I have the cidx^th element of y = L\b.
      #so compute the corresponding element of
      #this row of L and put it into the right place
      Li[tmpIdx] = k;
      Lx[tmpIdx] = yVals_cidx *Dinv[cidx];

      #D[k] -= yVals[cidx]*yVals[cidx]*Dinv[cidx];
      D[k] -= yVals_cidx*Lx[tmpIdx];
      LNextSpaceInCol[cidx]++;

      #reset the yvalues and indices back to zero and QDLDL_UNUSED
      #once I'm done with them
      yVals[cidx]     = 0.0;
      yMarkers[cidx]  = QDLDL_UNUSED;

    } #end for i

    #Maintain a count of the positive entries
    #in D.  If we hit a zero, we can't factor
    #this matrix, so abort
    if(D[k] == 0.0){return -1;}
    if(D[k]  > 0.0){positiveValuesInD++;}

    #compute the inverse of the diagonal
    Dinv[k]= 1/D[k];

  } #end for k

  return positiveValuesInD;

}

# Solves (L+I)x = b, with x replacing b
function QDLDL_Lsolve!(L,x)

  for i 1:n
      for j = Lp[i]: Lp[i+1]
          x[Li[j]] -= Lx[j]*x[i];
      end
  end
  return nothing
end


# Solves (L+I)x = b.  Returns x, preserving b
function QDLDL_Lsolve(L,b)
   x = copy(b)
   QDLDL_Lsolve!(L,x)
   return x
end


# Solves (L+I)'x = b, with x replacing b
function QDLDL_Ltsolve!(L,x)

  for(i = n:-1:1){
      for(j = Lp[i]:Lp[i+1])
          x[i] -= Lx[j]*x[Li[j]]
      end
  end
  return nothing
end

# Solves (L+I)'x = b.  Returns x, preserving b
function QDLDL_Ltsolve(L,b)
   x = copy(b)
   QDLDL_Ltsolve!(L,x)
   return x
end
