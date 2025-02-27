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
 */

/* QDLDL TESTER MODULE */

/* THE CODE FOR MINIMAL UNIT TESTING HAS BEEN TAKEN FROM
   http://www.jera.com/techinfo/jtns/jtn002.html */

#include <stdio.h>
#include <stdlib.h>
#include "minunit.h"
#include "qdldl.h"

// Utility functions for solves
QDLDL_float vec_diff_norm(QDLDL_float* x, QDLDL_float* y, QDLDL_int len);
int ldl_factor_solve(QDLDL_int An, QDLDL_int* Ap, QDLDL_int* Ai, QDLDL_float* Ax, QDLDL_float* b);

// Include tests
#include "test_basic.h"
#include "test_identity.h"
#include "test_rank_deficient.h"
#include "test_singleton.h"
#include "test_sym_structure.h"
#include "test_tril_structure.h"
#include "test_two_by_two.h"
#include "test_zero_on_diag.h"
#include "test_osqp_kkt.h"


int tests_run = 0;


static char* all_tests() {
    mu_run_test(test_basic);
    mu_run_test(test_identity);
    mu_run_test(test_rank_deficient);
    mu_run_test(test_singleton);
    mu_run_test(test_sym_structure);
    mu_run_test(test_tril_structure);
    mu_run_test(test_two_by_two);
    mu_run_test(test_zero_on_diag);
    mu_run_test(test_osqp_kkt);

    return 0;
}


QDLDL_float vec_diff_norm(QDLDL_float* x, QDLDL_float* y, QDLDL_int len) {
    QDLDL_float maxDiff = 0.0;
    QDLDL_float elDiff = 0.0;
    QDLDL_int   i = 0;

    for(i = 0; i < len; i++) {
        elDiff = x[i] - y[i];
        maxDiff = (elDiff > maxDiff) ? elDiff : ((-elDiff > maxDiff) ? -elDiff : maxDiff);
    }
    return maxDiff;
}

int ldl_factor_solve(QDLDL_int An, QDLDL_int* Ap, QDLDL_int* Ai, QDLDL_float* Ax, QDLDL_float* b) {
    // Data for L and D factors
    QDLDL_int    Ln = An;
    QDLDL_int*   Lp = 0;
    QDLDL_int*   Li = 0;
    QDLDL_float* Lx;
    QDLDL_float* D;
    QDLDL_float* Dinv;

    // Data for elim tree calculation
    QDLDL_int* etree;
    QDLDL_int* Lnz;
    QDLDL_int  sumLnz = 0;

    // Data for factorisation
    QDLDL_int*   iwork;
    QDLDL_bool*  bwork;
    QDLDL_float* fwork;
    QDLDL_int    factorStatus = 0;

    /*--------------------------------
     * pre-factorisation memory allocations
     *---------------------------------*/

    // These can happen *before* the etree is calculated
    // since the sizes are not sparsity pattern specific

    // For the elimination tree
    etree = (QDLDL_int*) malloc(sizeof(QDLDL_int) * An);
    Lnz = (QDLDL_int*) malloc(sizeof(QDLDL_int) * An);

    // For the L factors.   Li and Lx are sparsity dependent
    // so must be done after the etree is constructed
    Lp = (QDLDL_int*) malloc(sizeof(QDLDL_int) * (An + 1));
    D = (QDLDL_float*) malloc(sizeof(QDLDL_float) * An);
    Dinv = (QDLDL_float*) malloc(sizeof(QDLDL_float) * An);

    // Working memory.  Note that both the etree and factor
    // calls requires a working vector of QDLDL_int, with
    // the factor function requiring 3*An elements and the
    // etree only An elements.   Just allocate the larger
    // amount here and use it in both places
    iwork = (QDLDL_int*) malloc(sizeof(QDLDL_int) * (3 * An));
    bwork = (QDLDL_bool*) malloc(sizeof(QDLDL_bool) * An);
    fwork = (QDLDL_float*) malloc(sizeof(QDLDL_float) * An);


    /*--------------------------------
     * Elimination tree calculation
     *---------------------------------*/
    sumLnz = QDLDL_etree(An, Ap, Ai, iwork, Lnz, etree);

    // Not perfect triu A = bomb
    if(sumLnz < 0) {
        free(Lp);
        free(D);
        free(Dinv);
        free(etree);
        free(Lnz);
        free(iwork);
        free(bwork);
        free(fwork);
        return sumLnz;
    }

    /*--------------------------------
     * LDL factorisation
     *---------------------------------*/

    // First allocate memory for Li and Lx
    Li = (QDLDL_int*) malloc(sizeof(QDLDL_int) * sumLnz);
    Lx = (QDLDL_float*) malloc(sizeof(QDLDL_float) * sumLnz);

    // Now factor
    factorStatus =
            QDLDL_factor(An, Ap, Ai, Ax, Lp, Li, Lx, D, Dinv, Lnz, etree, bwork, iwork, fwork);

    // Zero on the diagonal = bomb
    if(factorStatus < 0) {
        free(Lp);
        free(D);
        free(Dinv);
        free(etree);
        free(Lnz);
        free(iwork);
        free(bwork);
        free(fwork);
        free(Li);
        free(Lx);
        return factorStatus;
    }

    /*--------------------------------
     * solve
     *---------------------------------*/
    QDLDL_solve(Ln, Lp, Li, Lx, Dinv, b);


    /*--------------------------------
     * clean up
     *---------------------------------*/
    free(Lp);
    free(D);
    free(Dinv);
    free(etree);
    free(Lnz);
    free(iwork);
    free(bwork);
    free(fwork);
    free(Li);
    free(Lx);

    return 0;
}


int main(void) {
    char* result = all_tests();

    if(result != 0) {
        printf("%s\n", result);
    } else {
        printf("ALL TESTS PASSED\n");
    }
    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
