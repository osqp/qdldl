/*
 * This file is part of QDLDL, a library for performing the LDL^T factorization
 * of a symmetric indefinite matrix.
 *
 * QDLDL is part of the OSQP project, and is available at https://github.com/osqp/qdldl.
 *
 * Copyright 2018, Paul Goulart, Bartolomeo Stellato, The OSQP developers
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

static char* test_basic()
{
  //A matrix data
  QDLDL_int Ap[]  = {0, 1, 2, 4, 5, 6, 8, 10, 12, 14, 17};
  QDLDL_int Ai[]  = {0, 1, 1, 2, 3, 4, 1, 5, 0, 6, 3, 7, 6, 8, 1, 2, 9};
  QDLDL_float Ax[] = {1.0, 0.460641, -0.121189, 0.417928, 0.177828,
                      0.1, -0.0290058, -1.0, 0.350321, -0.441092, -0.0845395,
                      -0.316228, 0.178663, -0.299077, 0.182452, -1.56506, -0.1};
  QDLDL_int An = 10;

  // RHS and solution to Ax = b
  QDLDL_float b[]    = {1,2,3,4,5,6,7,8,9,10};
  QDLDL_float xsol[] = {10.2171, 3.9416, -5.69096, 9.28661, 50.0, -6.11433,
                     -26.3104, -27.7809, -45.8099, -3.74178};

  //x replaces b during solve
  int status = ldl_factor_solve(An,Ap,Ai,Ax,b);

  mu_assert("Factorisation failed", status >= 0);
  mu_assert("Solve accuracy failed", vec_diff_norm(b,xsol,An) < QDLDL_TESTS_TOL);

  return 0;
}
