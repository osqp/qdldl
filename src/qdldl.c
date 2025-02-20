/*
 * This file is part of QDLDL, a library for performing the LDL^T factorization
 * of a symmetric indefinite matrix.
 *
 * QDLDL is part of the OSQP project, and is available at https://github.com/osqp/qdldl.
 *
 * Copyright 2018, Paul Goulart, Bartolomeo Stellato, Goran Banjac, The OSQP developers
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * SPDX-License-Identifier: Apache-2.0
 * SPDX-ExternalRef: PACKAGE_MANAGER purl pkg:github/osqp/qdldl
 */
#include "qdldl.h"

#define QDLDL_UNKNOWN (-1)
#define QDLDL_USED (1)
#define QDLDL_UNUSED (0)

/* Compute the elimination tree for a quasidefinite matrix
 * in compressed sparse column form.
 */
QDLDL_int QDLDL_etree(const QDLDL_int n, const QDLDL_int* Ap, const QDLDL_int* Ai, QDLDL_int* work,
                      QDLDL_int* Lnz, QDLDL_int* etree) {
    QDLDL_int i = 0;
    QDLDL_int j = 0;
    QDLDL_int p = 0;
    QDLDL_int sumLnz = 0;

    for(i = 0; i < n; i++) {
        // Zero out Lnz and work.  Set all etree values to unknown
        work[i]  = 0;
        Lnz[i]   = 0;
        etree[i] = QDLDL_UNKNOWN;

        // Abort if A doesn't have at least one entry
        // one entry in every column
        if(Ap[i] == Ap[i + 1]) {
            return -1;
        }
    }

    for(j = 0; j < n; j++) {
        work[j] = j;

        for(p = Ap[j]; p < Ap[j + 1]; p++) {
            i = Ai[p];

            // Abort if entries on lower triangle
            if(i > j) {
                return -1;
            };

            while(work[i] != j) {
                if(etree[i] == QDLDL_UNKNOWN) {
                    etree[i] = j;
                }
                Lnz[i]++; // Nonzeros in this column
                work[i] = j;
                i = etree[i];
            }
        }
    }

    // Compute the total nonzeros in L.  This much
    // space is required to store Li and Lx.  Return
    // error code -2 if the nonzero count will overflow
    // its unteger type.
    sumLnz = 0;

    for(i = 0; i < n; i++) {
        if(sumLnz > QDLDL_INT_MAX - Lnz[i]) {
            sumLnz = -2;
            break;
        } else {
            sumLnz += Lnz[i];
        }
    }

    return sumLnz;
}


QDLDL_int QDLDL_factor(const QDLDL_int n, const QDLDL_int* Ap, const QDLDL_int* Ai,
                       const QDLDL_float* Ax, QDLDL_int* Lp, QDLDL_int* Li, QDLDL_float* Lx,
                       QDLDL_float* D, QDLDL_float* Dinv, const QDLDL_int* Lnz,
                       const QDLDL_int* etree, QDLDL_bool* bwork, QDLDL_int* iwork,
                       QDLDL_float* fwork) {
    QDLDL_int    i = 0;
    QDLDL_int    j = 0;
    QDLDL_int    k = 0;
    QDLDL_int    nnzY = 0;
    QDLDL_int    bidx = 0;
    QDLDL_int    cidx = 0;
    QDLDL_int    nextIdx = 0;
    QDLDL_int    nnzE = 0;
    QDLDL_int    tmpIdx = 0;
    QDLDL_int*   yIdx;
    QDLDL_int*   elimBuffer;
    QDLDL_int*   LNextSpaceInCol;
    QDLDL_float* yVals;
    QDLDL_float  yVals_cidx = 0.0;
    QDLDL_bool*  yMarkers;
    QDLDL_int    positiveValuesInD = 0;

    // Partition working memory into pieces
    yMarkers = bwork;
    yIdx = iwork;
    elimBuffer = iwork + n;
    LNextSpaceInCol = iwork + n * 2;
    yVals = fwork;


    Lp[0] = 0; // First column starts at index zero

    for(i = 0; i < n; i++) {
        // Compute L column indices
        Lp[i + 1] = Lp[i] + Lnz[i]; // cumsum, total at the end

        // Set all Yidx to be 'unused' initially
        // in each column of L, the next available space
        // to start is just the first space in the column
        yMarkers[i] = QDLDL_UNUSED;
        yVals[i] = 0.0;
        D[i] = 0.0;
        LNextSpaceInCol[i] = Lp[i];
    }

    // First element of the diagonal D.
    D[0] = Ax[0];

    if(D[0] == 0.0) {
        return -1;
    }

    if(D[0] > 0.0) {
        positiveValuesInD++;
    }
    Dinv[0] = 1 / D[0];

    // Start from 1 here. The upper LH corner is trivially 0
    // in L b/c we are only computing the subdiagonal elements
    for(k = 1; k < n; k++) {
        // NB : For each k, we compute a solution to
        // y = L(0:(k-1),0:k-1))\b, where b is the kth
        // column of A that sits above the diagonal.
        // The solution y is then the kth row of L,
        // with an implied '1' at the diagonal entry.

        // Number of nonzeros in this row of L
        nnzY = 0; // Number of elements in this row

        // This loop determines where nonzeros
        // will go in the kth row of L, but doesn't
        // compute the actual values
        tmpIdx = Ap[k + 1];

        for(i = Ap[k]; i < tmpIdx; i++) {
            bidx = Ai[i]; // We are working on this element of b

            // Initialize D[k] as the element of this column
            // corresponding to the diagonal place.  Don't use
            // this element as part of the elimination step
            // that computes the k^th row of L
            if(bidx == k) {
                D[k] = Ax[i];
                continue;
            }

            yVals[bidx] = Ax[i]; // Initialise y(bidx) = b(bidx)

            // Use the forward elimination tree to figure
            // out which elements must be eliminated after
            // this element of b
            nextIdx = bidx;

            if(yMarkers[nextIdx] == QDLDL_UNUSED) { // This y term not already visited

                yMarkers[nextIdx] = QDLDL_USED; // I touched this one
                elimBuffer[0] = nextIdx;        // It goes at the start of the current list
                nnzE = 1;                       // Length of unvisited elimination path from here

                nextIdx = etree[bidx];

                while(nextIdx != QDLDL_UNKNOWN && nextIdx < k) {
                    if(yMarkers[nextIdx] == QDLDL_USED)
                        break;

                    yMarkers[nextIdx] = QDLDL_USED; // I touched this one
                    elimBuffer[nnzE] = nextIdx;     // It goes in the current list
                    nnzE++;                         // The list is one longer than before
                    nextIdx = etree[nextIdx];       // One step further along tree

                }

                // Now I put the buffered elimination list into
                // my current ordering in reverse order
                while(nnzE) {
                    yIdx[nnzY++] = elimBuffer[--nnzE];
                }
            }
        }

        // This for loop places nonzeros values in the k^th row
        for(i = (nnzY - 1); i >= 0; i--) {
            //Which column are we working on?
            cidx = yIdx[i];

            // Loop along the elements in this
            // column of L and subtract to solve to y
            tmpIdx = LNextSpaceInCol[cidx];
            yVals_cidx = yVals[cidx];

            for(j = Lp[cidx]; j < tmpIdx; j++) {
                yVals[Li[j]] -= Lx[j] * yVals_cidx;
            }

            // Now I have the cidx^th element of y = L\b.
            // so compute the corresponding element of
            // this row of L and put it into the right place
            Li[tmpIdx] = k;
            Lx[tmpIdx] = yVals_cidx * Dinv[cidx];

            // D[k] -= yVals[cidx]*yVals[cidx]*Dinv[cidx];
            D[k] -= yVals_cidx * Lx[tmpIdx];
            LNextSpaceInCol[cidx]++;

            // Reset the yvalues and indices back to zero and QDLDL_UNUSED
            // once I'm done with them
            yVals[cidx] = 0.0;
            yMarkers[cidx] = QDLDL_UNUSED;

        }

        // Maintain a count of the positive entries
        // in D.  If we hit a zero, we can't factor
        // this matrix, so abort
        if(D[k] == 0.0) {
            return -1;
        }

        if(D[k] > 0.0) {
            positiveValuesInD++;
        }

        // Compute the inverse of the diagonal
        Dinv[k] = 1 / D[k];

    }

    return positiveValuesInD;
}

// Solves (L+I)x = b
void QDLDL_Lsolve(const QDLDL_int n, const QDLDL_int* Lp, const QDLDL_int* Li,
                  const QDLDL_float* Lx, QDLDL_float* x) {
    QDLDL_int i = 0;
    QDLDL_int j = 0;

    for(i = 0; i < n; i++) {
        QDLDL_float val = x[i];

        for(j = Lp[i]; j < Lp[i + 1]; j++) {
            x[Li[j]] -= Lx[j] * val;
        }
    }
}

// Solves (L+I)'x = b
void QDLDL_Ltsolve(const QDLDL_int n, const QDLDL_int* Lp, const QDLDL_int* Li,
                   const QDLDL_float* Lx, QDLDL_float* x) {
    QDLDL_int i = 0;
    QDLDL_int j = 0;

    for(i = n - 1; i >= 0; i--) {
        QDLDL_float val = x[i];

        for(j = Lp[i]; j < Lp[i + 1]; j++) {
            val -= Lx[j] * x[Li[j]];
        }
        x[i] = val;
    }
}

// Solves Ax = b where A has given LDL factors
void QDLDL_solve(const QDLDL_int n, const QDLDL_int* Lp, const QDLDL_int* Li, const QDLDL_float* Lx,
                 const QDLDL_float* Dinv, QDLDL_float* x) {
    QDLDL_int i = 0;

    QDLDL_Lsolve(n, Lp, Li, Lx, x);

    for(i = 0; i < n; i++) {
        x[i] *= Dinv[i];
    }

    QDLDL_Ltsolve(n, Lp, Li, Lx, x);
}
