/*
 * This file is part of QDLDL, a library for performing the LDL^T factorization
 * of a symmetric indefinite matrix.
 *
 * QDLDL is part of the OSQP project, and is available at https://github.com/osqp/qdldl.
 *
 * Copyright 2018, Paul Goulart, The OSQP developers
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

static char* test_zero_on_diag()
{
  //A matrix data
  QDLDL_int Ap[]  = {0, 1, 2, 5};
  QDLDL_int Ai[]  = {0 ,0, 0, 1, 2};
  QDLDL_float Ax[] = {4,1,2,1,-3};
  QDLDL_int An = 3;

  // RHS and solution to Ax = b
  QDLDL_float b[]    = {6,9,12};
  QDLDL_float xsol[] = {17,-46,-8};

  //x replaces b during solve (should fill due to zero in middle)
  //NB : this system is solvable, but not by LDL
  int status = ldl_factor_solve(An,Ap,Ai,Ax,b);

  mu_assert("Factorisation failed", status >= 0);
  mu_assert("Solve accuracy failed", vec_diff_norm(b,xsol,An) < QDLDL_TESTS_TOL);

  return 0;
}
