/*
 * This file is part of QDLDL, a library for performing the LDL^T factorization
 * of a symmetric indefinite matrix.
 *
 * QDLDL is part of the OSQP project, and is available at https://github.com/osqp/qdldl.
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

/* Test that QDLDL_symbolic + QDLDL_numeric produce the same result as
 * QDLDL_factor, and that QDLDL_numeric can refactorize with new values. */

static char* test_symbolic() {
    // A matrix data (same as test_basic)
    QDLDL_int   Ap[] = { 0, 1, 2, 4, 5, 6, 8, 10, 12, 14, 17 };
    QDLDL_int   Ai[] = { 0, 1, 1, 2, 3, 4, 1, 5, 0, 6, 3, 7, 6, 8, 1, 2, 9 };
    QDLDL_float Ax[] = { 1.0,        0.460641,  -0.121189, 0.417928,  0.177828,   0.1,
                         -0.0290058, -1.0,      0.350321,  -0.441092, -0.0845395, -0.316228,
                         0.178663,   -0.299077, 0.182452,  -1.56506,  -0.1 };
    QDLDL_int   An = 10;

    // RHS and solution to Ax = b
    QDLDL_float b[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
    QDLDL_float xsol[] = { 10.2171,  3.9416,   -5.69096, 9.28661,  50.0,
                           -6.11433, -26.3104, -27.7809, -45.8099, -3.74178 };

    // Allocate all working memory
    QDLDL_int*   etree = (QDLDL_int*) malloc(sizeof(QDLDL_int) * An);
    QDLDL_int*   Lnz   = (QDLDL_int*) malloc(sizeof(QDLDL_int) * An);
    QDLDL_int*   Lp    = (QDLDL_int*) malloc(sizeof(QDLDL_int) * (An + 1));
    QDLDL_float* D     = (QDLDL_float*) malloc(sizeof(QDLDL_float) * An);
    QDLDL_float* Dinv  = (QDLDL_float*) malloc(sizeof(QDLDL_float) * An);
    QDLDL_int*   iwork = (QDLDL_int*) malloc(sizeof(QDLDL_int) * (3 * An));
    QDLDL_bool*  bwork = (QDLDL_bool*) malloc(sizeof(QDLDL_bool) * An);
    QDLDL_float* fwork = (QDLDL_float*) malloc(sizeof(QDLDL_float) * An);

    QDLDL_int sumLnz = QDLDL_etree(An, Ap, Ai, iwork, Lnz, etree);
    mu_assert("etree failed", sumLnz >= 0);

    QDLDL_int*   Li = (QDLDL_int*) malloc(sizeof(QDLDL_int) * sumLnz);
    QDLDL_float* Lx = (QDLDL_float*) malloc(sizeof(QDLDL_float) * sumLnz);

    // Step 1: symbolic factorization
    QDLDL_symbolic(An, Ap, Ai, Lp, Li, Lnz, etree, bwork, iwork);

    // Step 2: numeric factorization using pre-computed pattern
    QDLDL_int factorStatus = QDLDL_numeric(An, Ap, Ai, Ax, Lp, Li, Lx, D, Dinv,
                                           Lnz, etree, bwork, iwork, fwork);
    mu_assert("Numeric factorisation failed", factorStatus >= 0);

    // Solve and check
    QDLDL_int i;
    for(i = 0; i < An; i++) b[i] = i + 1;
    QDLDL_solve(An, Lp, Li, Lx, Dinv, b);
    mu_assert("Solve accuracy failed after symbolic+numeric",
              vec_diff_norm(b, xsol, An) < QDLDL_TESTS_TOL);

    // Step 3: refactorize with same values to verify reuse works
    factorStatus = QDLDL_numeric(An, Ap, Ai, Ax, Lp, Li, Lx, D, Dinv,
                                 Lnz, etree, bwork, iwork, fwork);
    mu_assert("Refactorisation failed", factorStatus >= 0);

    for(i = 0; i < An; i++) b[i] = i + 1;
    QDLDL_solve(An, Lp, Li, Lx, Dinv, b);
    mu_assert("Solve accuracy failed after refactorisation",
              vec_diff_norm(b, xsol, An) < QDLDL_TESTS_TOL);

    free(etree); free(Lnz); free(Lp); free(D); free(Dinv);
    free(iwork); free(bwork); free(fwork); free(Li); free(Lx);

    return 0;
}

/* Test symbolic + numeric on a 2x2 system, then refactorize with different values. */
static char* test_symbolic_refactor() {
    // A = [1, 1; 1, -1]
    QDLDL_int   Ap[] = { 0, 1, 3 };
    QDLDL_int   Ai[] = { 0, 0, 1 };
    QDLDL_float Ax1[] = { 1.0, 1.0, -1.0 };
    QDLDL_int   An = 2;

    QDLDL_float b1[] = { 2, 4 };
    QDLDL_float xsol1[] = { 3, -1 };

    // Second set of values with same sparsity: A = [2, 3; 3, -2]
    QDLDL_float Ax2[] = { 2.0, 3.0, -2.0 };
    QDLDL_float b2[] = { 1, 1 };
    // A = [2 3; 3 -2], b = [1;1]
    // det = -4-9 = -13, A^{-1} = (1/-13)*[-2 -3; -3 2]
    // x = [5/13; 1/13]
    QDLDL_float xsol2[] = { 5.0 / 13.0, 1.0 / 13.0 };

    QDLDL_int*   etree = (QDLDL_int*) malloc(sizeof(QDLDL_int) * An);
    QDLDL_int*   Lnz   = (QDLDL_int*) malloc(sizeof(QDLDL_int) * An);
    QDLDL_int*   Lp    = (QDLDL_int*) malloc(sizeof(QDLDL_int) * (An + 1));
    QDLDL_float* D     = (QDLDL_float*) malloc(sizeof(QDLDL_float) * An);
    QDLDL_float* Dinv  = (QDLDL_float*) malloc(sizeof(QDLDL_float) * An);
    QDLDL_int*   iwork = (QDLDL_int*) malloc(sizeof(QDLDL_int) * (3 * An));
    QDLDL_bool*  bwork = (QDLDL_bool*) malloc(sizeof(QDLDL_bool) * An);
    QDLDL_float* fwork = (QDLDL_float*) malloc(sizeof(QDLDL_float) * An);

    QDLDL_int sumLnz = QDLDL_etree(An, Ap, Ai, iwork, Lnz, etree);
    mu_assert("etree failed", sumLnz >= 0);

    QDLDL_int*   Li = (QDLDL_int*) malloc(sizeof(QDLDL_int) * sumLnz);
    QDLDL_float* Lx = (QDLDL_float*) malloc(sizeof(QDLDL_float) * sumLnz);

    // Symbolic once
    QDLDL_symbolic(An, Ap, Ai, Lp, Li, Lnz, etree, bwork, iwork);

    // Numeric with first values
    QDLDL_int factorStatus = QDLDL_numeric(An, Ap, Ai, Ax1, Lp, Li, Lx, D, Dinv,
                                           Lnz, etree, bwork, iwork, fwork);
    mu_assert("First factorisation failed", factorStatus >= 0);

    QDLDL_solve(An, Lp, Li, Lx, Dinv, b1);
    mu_assert("First solve accuracy failed",
              vec_diff_norm(b1, xsol1, An) < QDLDL_TESTS_TOL);

    // Numeric with second values (same pattern, different numbers)
    factorStatus = QDLDL_numeric(An, Ap, Ai, Ax2, Lp, Li, Lx, D, Dinv,
                                 Lnz, etree, bwork, iwork, fwork);
    mu_assert("Second factorisation failed", factorStatus >= 0);

    QDLDL_solve(An, Lp, Li, Lx, Dinv, b2);
    mu_assert("Second solve accuracy failed",
              vec_diff_norm(b2, xsol2, An) < QDLDL_TESTS_TOL);

    free(etree); free(Lnz); free(Lp); free(D); free(Dinv);
    free(iwork); free(bwork); free(fwork); free(Li); free(Lx);

    return 0;
}
