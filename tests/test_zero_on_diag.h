static char* test_zero_on_diag()
{
  //A matrix data
  QDLDL_int Ap[]  = {0, 1, 2, 5};
  QDLDL_int Ai[]  = {0 ,0, 0, 1, 2};
  QDLDL_float Ax[] = {4,1,2,1,-3};
  QDLDL_int An = 3;

  // RHS and solution to Ax = b
  QDLDL_float b[]    = {1,2,3};
  QDLDL_float xsol[] = {4,-11,-2};

  //x replaces b during solve (should fill due to zero in middle)
  //NB : this system is solvable, but not by LDL
  int status = ldl_factor_solve(An,Ap,Ai,Ax,b);

  mu_assert("Factorisation failed", status < 0);

  return 0;
}
