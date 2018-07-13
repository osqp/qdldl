#ifndef OS_LDL_H
# define OS_LDL_H

//Define external bool, int and float types if they
//are not already defined externally.  If you wish
//to have your own types defined here, then define
//OS_LDL_TYPES_DEFINED elsewhere and manually define
//the types os_ldl_bool, os_ldl_int and os_ldl_float

#ifndef OS_LDL_TYPES_DEFINED
  #include <stdbool.h>
  #define OS_LDL_TYPES_DEFINED
  typedef os_ldl_bool   bool;
  typedef os_ldl_int    long long;
  typedef os_ldl_float  double;
# endif


# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

/**
 *  Matrix in compressed-column or triplet form.  The same structure
 *  is used for csc, csr and sparse triplet form
 */
struct SparseMatrix_ {
  os_ldl_int    nzmax; ///< maximum number of entries.
  os_ldl_int    m;     ///< number of rows
  os_ldl_int    n;     ///< number of columns
  os_ldl_int   *p;     ///< column or row pointers (size n+1) (col indices (size nzmax)
                  // start from 0 when using triplet format (direct KKT matrix
                  // formation))
  os_ldl_int   *i;     ///< row indices, size nzmax starting from 0
  os_ldl_float *x;     ///< numerical values, size nzmax
};

typedef struct SparseMatrix_ CscMatrix; // Compressed sparse column matrix


/**
  * Compute the elimination tree for a quasidefinite matrix
  * in compressed sparse column form, where the input matrix is
  * assumed to contain data for the upper triangular part of A only,
  * and there are no duplicate indices.
  *
  * Returns an elimination tree for the factorization A = LDL^T and a
  * count of the nonzeros in each column of L that are strictly below the
  * diagonal.
  *
  * Does not use MALLOC.  It is assumed that the arrays work, Lnz, and
  * etree will be allocated with a number of elements equal to A->n.
  *
  * @param  A      CscMatrix (upper triangular part only)
  * @param  work   work vector (no meaning on return)
  * @param  Lnz    count of nonzeros in each column of L below diagonal
  * @param  etree  elimination tree
  * @return total  sum of Lnz (i.e. total nonzeros in L below diagonal)
  *
  *
*/

 os_ldl_int os_ldl_etree(const CscMatrix *A,
                os_ldl_int* work,
                os_ldl_int* Lnz,
                os_ldl_int* etree);

/**
  * Compute an LDL decomposition for a quasidefinite matrix
  * in compressed sparse column form, where the input matrix is
  * assumed to contain data for the upper triangular part of A only,
  * and there are no duplicate indices.
  *
  * Returns factors L, D and Dinv = 1./D.
  *
  * Does not use MALLOC.  It is assumed that L will be a CscMatrix
  * with sufficient space allocated, with a number of nonzeros equal
  * to the count given as a return value by osqp_ldl_etree
  *
  * @param  A      CscMatrix (upper triangular part only) (not modified)
  * @param  L      CscMatrix (size determined by osqp_ldl_etree
  * @param  D      vectorized factor D.  Length is A->n
  * @param  Dinv   reciprocal of D.  Length is A->n
  * @param  Lnz    count of nonzeros in each column of L below diagonal,
  *                as given by osqp_ldl_etree (not modified)
  * @param  etree  elimination tree as as given by osqp_ldl_etree (not modified)
  * @param  bwork  working array of bools. Length is A->n
  * @param  iwork  working array of integers. Length is 3*A->n
  * @return fwork  working array of floats. Length is A->n
  *
  *
*/


void os_ldl_factor(CscMatrix *A,
                     CscMatrix *L,
                     os_ldl_float* D,
                     os_ldl_float* Dinv,
                     os_ldl_int* Lnz,
                     os_ldl_int* etree,
                     os_ldl_bool* bwork,
                     os_ldl_int* iwork,
                     os_ldl_float* fwork);


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef OSQP_LDL_H
