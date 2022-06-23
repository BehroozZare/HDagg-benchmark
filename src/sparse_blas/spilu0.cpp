//
// Created by behrooz on 10/3/21.
//

#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>
#include <sparse_blas_lib.h>
#include <omp.h>

namespace sym_lib
{

    inline double sparse_dot_product(int l1, int u1, int l2, int u2, const int* indices, const double* data)
    {
        double result = 0.0;
        while (l1 < u1 && l2 < u2) {
            if (indices[l1] == indices[l2])     // matching column?
                result += data[l1++] * data[l2++];
            else if (indices[l1] < indices[l2])       // else proceed until we find matching columns
                l1++;
            else
                l2++;
        }
        return result;
    }

/*
  Purpose:

    REARRANGE_CR sorts a sparse compressed row matrix.

  Discussion:

    This routine guarantees that the entries in the CR matrix
    are properly sorted.

    After the sorting, the entries of the matrix are rearranged in such
    a way that the entries of each column are listed in ascending order
    of their column values.

    The matrix A is assumed to be stored in compressed row format.  Only
    the nonzero entries of A are stored.  The vector JA stores the
    column index of the nonzero value.  The nonzero values are sorted
    by row, and the compressed row vector IA then has the property that
    the entries in A and JA that correspond to row I occur in indices
    IA[I] through IA[I+1]-1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 July 2007

  Author:

    Lili Ju

  Parameters:

    Input, int N, the order of the system.

    Input, int NZ_NUM, the number of nonzeros.

    Input, int IA[N+1], the compressed row index vector.

    Input/output, int JA[NZ_NUM], the column indices of the matrix values.
    On output, the order of the entries of JA may have changed because of
    the sorting.

    Input/output, double A[NZ_NUM], the matrix values.  On output, the
    order of the entries may have changed because of the sorting.
*/
    void rearrange_cr(int n, int nz_num, int ia[], int ja[], double a[])
    {
        double dtemp;
        int i;
        int is;
        int itemp;
        int j;
        int j1;
        int j2;
        int k;

        for (i = 0; i < n; i++) {
            j1 = ia[i];
            j2 = ia[i + 1];
            is = j2 - j1;

            for (k = 1; k < is; k++) {
                for (j = j1; j < j2 - k; j++) {
                    if (ja[j + 1] < ja[j]) {
                        itemp = ja[j + 1];
                        ja[j + 1] = ja[j];
                        ja[j] = itemp;

                        dtemp = a[j + 1];
                        a[j + 1] = a[j];
                        a[j] = dtemp;
                    }
                }
            }
        }
        return;
    }


/******************************************************************************/
/*
  Purpose:

    ILU_CR computes the incomplete LU factorization of a matrix.

  Discussion:

    The matrix A is assumed to be stored in compressed row format.  Only
    the nonzero entries of A are stored.  The vector JA stores the
    column index of the nonzero value.  The nonzero values are sorted
    by row, and the compressed row vector IA then has the property that
    the entries in A and JA that correspond to row I occur in indices
    IA[I] through IA[I+1]-1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 July 2007

  Author:

    Lili Ju
*/
  //Parameters:
  //=============================== Up looking ==============================
    void spilu0_ul_csr_serial(int n, int nnz, const int *Ap, const int *Ai,
                              const double *Ax, const int *A_diag, double *l)
    {
        std::vector<int> iw(n);
        for (int j = 0; j < n; j++ ){
            iw[j] = EMPTY;
        }
        int jrow=0, jw;
        double tl;
        // copying A
        for (int k = 0; k < nnz; k++ ){
            l[k] = Ax[k];
        }
        for (int i = 0; i < n; i++ ){
            // iw points to the nonzero entries in row i.
            for (int k = Ap[i]; k < Ap[i+1]; k++ ){
                iw[Ai[k]] = k;
            }
            int j;
            for (j = Ap[i]; j < Ap[i+1]; ++j) {
                jrow = Ai[j];
                if ( i <= jrow ){
                    break;
                }
                tl = l[j] / l[A_diag[jrow]];
                l[j] = tl;
                for (int jj = A_diag[jrow] + 1; jj < Ap[jrow+1]; jj++ ){
                    jw = iw[Ai[jj]];
                    if ( jw != EMPTY ){
                        l[jw] = l[jw] - tl * l[jj];
                    }
                }
            }
            assert(A_diag[i] == j);
            assert(jrow == i);//if not, it mean it does not have a diagonal value
            assert(l[j] != 0.0); //zero diagonal check
            for (int k = Ap[i]; k < Ap[i+1]; k++ ){
                iw[Ai[k]] = EMPTY;
            }
        }
    }



    void spilu0_ul_csr_levelset(int n, int nnz, const int *Ap, const int *Ai,
                                const double *Ax, const int *A_diag, double *l,
                                int nLevels, const int* LevelPtr, const int* LevelSet, int* temp){

        /*
            Copy A.
        */
        #pragma omp parallel
        {
//            std::vector<int> iw(n);
//            for (int j = 0; j < n; j++ ){
//                iw[j] = EMPTY;
//            }
            int id = omp_get_thread_num();
            int *iw = temp+(id*n);
            int jrow=0, jw;
            double tl;
            // copying A
            // copying A
            #pragma omp for
            for (int k = 0; k < nnz; k++ ){
                l[k] = Ax[k];
            }

            for (int i1 = 0; i1 < nLevels; ++i1) {
                #pragma omp  for schedule(static)
                for (int j1 = LevelPtr[i1]; j1 < LevelPtr[i1 + 1]; ++j1) {
                    int i = LevelSet[j1];
                    // iw points to the nonzero entries in row i.
                    for (int k = Ap[i]; k < Ap[i+1]; k++ ){
                        iw[Ai[k]] = k;
                    }
                    int j;
                    for (j = Ap[i]; j < Ap[i+1]; ++j) {
                        jrow = Ai[j];
                        if ( i <= jrow ){
                            break;
                        }
                        tl = l[j] / l[A_diag[jrow]];
                        l[j] = tl;
                        for (int jj = A_diag[jrow] + 1; jj < Ap[jrow+1]; jj++ ){
                            jw = iw[Ai[jj]];
                            if ( jw != EMPTY ){
                                l[jw] = l[jw] - tl * l[jj];
                            }
                        }
                    }
                    assert(A_diag[i] == j);
                    assert(jrow == i);//if not, it mean it does not have a diagonal value
                    assert(l[j] != 0.0); //zero diagonal check
                    for (int k = Ap[i]; k < Ap[i+1]; k++ ){
                        iw[Ai[k]] = EMPTY;
                    }
                }
            }
        }
    }


    void spilu0_ul_csr_group_levelset(int n, int nnz, const int *Ap, const int *Ai,
                                const double *Ax, const int *A_diag, double *l,
                                int nLevels, const int* LevelPtr, const int* LevelSet,
                                const int* group_ptr, const int* group_set, int* temp){

        /*
            Copy A.
        */
        #pragma omp parallel
        {
            //            std::vector<int> iw(n);
            //            for (int j = 0; j < n; j++ ){
            //                iw[j] = EMPTY;
            //            }
            int id = omp_get_thread_num();
            int *iw = temp+(id*n);
            int jrow=0, jw;
            double tl;
            // copying A
            // copying A
            #pragma omp for
            for (int k = 0; k < nnz; k++ ){
                l[k] = Ax[k];
            }

            for (int i1 = 0; i1 < nLevels; ++i1) {
                #pragma omp  for schedule(static)
                for (int j1 = LevelPtr[i1]; j1 < LevelPtr[i1 + 1]; ++j1) {
                    int g_ptr = LevelSet[j1];
                    for(int g = group_ptr[g_ptr]; g < group_ptr[g_ptr + 1]; g++){
                        int i = group_set[g];
                        // iw points to the nonzero entries in row i.
                        for (int k = Ap[i]; k < Ap[i+1]; k++ ){
                            iw[Ai[k]] = k;
                        }
                        int j;
                        for (j = Ap[i]; j < Ap[i+1]; ++j) {
                            jrow = Ai[j];
                            if ( i <= jrow ){
                                break;
                            }
                            tl = l[j] / l[A_diag[jrow]];
                            l[j] = tl;
                            for (int jj = A_diag[jrow] + 1; jj < Ap[jrow+1]; jj++ ){
                                jw = iw[Ai[jj]];
                                if ( jw != EMPTY ){
                                    l[jw] = l[jw] - tl * l[jj];
                                }
                            }
                        }
                        assert(A_diag[i] == j);
                        assert(jrow == i);//if not, it mean it does not have a diagonal value
                        assert(l[j] != 0.0); //zero diagonal check
                        for (int k = Ap[i]; k < Ap[i+1]; k++ ){
                            iw[Ai[k]] = EMPTY;
                        }
                    }
                }
            }
        }
    }




    void spilu0_ul_csr_lbc(int n, int nnz, const int *Ap, const int *Ai,
                           const double *Ax, const int *A_diag, double *l,
                           int level_no, const int* level_ptr, const int* par_ptr, const int* partition, int* temp){

        /*
            Copy A.
        */
        #pragma omp parallel
        {
            int id = omp_get_thread_num();
            int *iw = temp+(id*n);
//            std::vector<int> iw(n);
//            for (int j = 0; j < n; j++ ){
//                iw[j] = EMPTY;
//            }
            int jrow=0, jw;
            double tl;
            // copying A
            // copying A
            #pragma omp for
            for (int k = 0; k < nnz; k++ ){
                l[k] = Ax[k];
            }

            for (int i1 = 0; i1 < level_no; ++i1) {
                #pragma omp  for schedule(static)
                for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1) {
                    for (int k1 = par_ptr[j1]; k1 < par_ptr[j1 + 1]; ++k1) {
                        int i = partition[k1];
                        // iw points to the nonzero entries in row i.
                        for (int k = Ap[i]; k < Ap[i + 1]; k++) {
                            iw[Ai[k]] = k;
                        }
                        int j;
                        for (j = Ap[i]; j < Ap[i + 1]; ++j) {
                            jrow = Ai[j];
                            if (i <= jrow) {
                                break;
                            }
                            tl = l[j] / l[A_diag[jrow]];
                            l[j] = tl;
                            for (int jj = A_diag[jrow] + 1; jj < Ap[jrow + 1]; jj++) {
                                jw = iw[Ai[jj]];
                                if (jw != EMPTY) {
                                    l[jw] = l[jw] - tl * l[jj];
                                }
                            }
                        }
                        assert(A_diag[i] == j);
                        assert(jrow == i);//if not, it mean it does not have a diagonal value
                        assert(l[j] != 0.0); //zero diagonal check
                        for (int k = Ap[i]; k < Ap[i + 1]; k++) {
                            iw[Ai[k]] = EMPTY;
                        }
                    }
                }
            }
        }
    }


} // namespace sym_lib