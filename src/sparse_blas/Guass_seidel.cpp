//
// Created by george on 2019-10-09.
//
#include <omp.h>
#include <vector>
#include <algorithm>
#include "sparse_blas_lib.h"

namespace sym_lib
{

    //=========================== Right Looking SpTrSv ==========================
    /*
     * It is a csr version of guass-seidel
     * @param n Number of iterations or node
     * @param Lp the pointer array in CSR version
     * @param Li the index array in CSR version
     * @param Lx the value array in CSR version
     * @return x the output
     */
    void gs_csr(int n, int *Lp, int *Li,int *idiag, double *values, double *y, const double *b)
    {
        int i,j;
        double sum;

        for (i = 0; i < n; i++) {
            sum = b[i];
            for (j = Lp[i]; j < Lp[i + 1]; j++) {
                sum -= values[j]*y[Li[j]];
            }
            y[i] = sum*idiag[i];
        } // for each row
    }
} // namespace sym_lib