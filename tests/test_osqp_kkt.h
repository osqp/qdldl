static char* test_osqp_kkt()
{
  //A matrix data
  QDLDL_int Ap[]  = {0, 1, 2, 5, 6, 7, 8, 12};
  QDLDL_int Ai[]  = {0, 1, 2, 1, 0, 3, 4, 5, 5, 6, 4, 3};
  QDLDL_float Ax[] = {-0.25000000,  -0.25000000,   1.00000000,   0.51357812,   0.52914209,  -0.25000000,  -0.25000000,   1.10274361,   0.15537975,   1.25882928,   0.13457995,   0.62113383};
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

