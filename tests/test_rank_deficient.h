static char* test_rank_deficient()
{
  //A matrix data
  QDLDL_int Ap[]  = {0, 1, 3};
  QDLDL_int Ai[]  = {0, 0, 1};
  QDLDL_float Ax[] = {1.0, 1.0, 1.0};
  QDLDL_int An = 2;

  // RHS and solution to Ax = b (should be impossible)
  QDLDL_float b[]    = {1,1};
  QDLDL_float xsol[] = {1,1};

  //x replaces b during solve
  int status = ldl_factor_solve(An,Ap,Ai,Ax,b);

  mu_assert("Rank deficiency not detected", status < 0);

  return 0;
}
