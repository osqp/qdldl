/*
 * This file is part of QDLDL, a library for performing the LDL^T factorization
 * of a symmetric indefinite matrix.
 *
 * QDLDL is part of the OSQP project, and is available at https://github.com/osqp/qdldl.
 *
 * Copyright 2018, Bartolomeo Stellato, The OSQP developers
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

static char* test_osqp_kkt()
{
  // Unordered A
  QDLDL_int Ap[]  = {0, 1, 2, 5, 6, 7, 8, 12};
  QDLDL_int Ai[]  = {0, 1, 2, 1, 0, 3, 4, 5, 5, 6, 4, 3};
  QDLDL_float Ax[] = {-0.25,  -0.25,   1.0,   0.513578,   0.529142,  -0.25,  -0.25,   1.10274,   0.15538,   1.25883,   0.13458,   0.621134};

  // Ordered A (this works)
  // QDLDL_int Ap[]  = {0, 1, 2, 5, 6, 7, 8, 12};
  // QDLDL_int Ai[]  = {0, 1, 0, 1, 2, 3, 4, 5, 3, 4, 5, 6};
  // QDLDL_float Ax[] = {-0.25, -0.25, 0.529142, 0.513578, 1.0, -0.25, -0.25, 1.10274, 0.621134, 0.13458, 0.15538, 1.25883};
  QDLDL_int An = 7;

  // RHS and solution to Ax = b
  QDLDL_float b[]    = {-0.595598, -0.0193715, -0.576156, -0.168746, 0.61543, 0.419073, 1.31087};
  QDLDL_float xsol[] = {1.13141, -1.1367, -0.591044, 1.68867, -2.24209, 0.32254, 0.407998};

  //x replaces b during solve
  int status = ldl_factor_solve(An,Ap,Ai,Ax,b);

  mu_assert("Factorisation failed", status >= 0);
  mu_assert("Solve accuracy failed", vec_diff_norm(b,xsol,An) < QDLDL_TESTS_TOL);

  return 0;
}

