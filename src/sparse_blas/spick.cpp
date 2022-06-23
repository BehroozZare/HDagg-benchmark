//
// Created by behrooz on 2021-11-04.
//

#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>

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


    //=============================== up Looking ==============================
    bool spick_csr_UL_serial(int n, const int *Lp, const int *Li, double *Lx,
                             double *new_Lp, double* new_Li, double* new_Lx, double& new_nnz){
        for (int i = 0; i < n; ++i) {
            for (int k = Lp[i]; k < Lp[i+1]; ++k) {
                const int j = Li[k];     // column
                double dp = sparse_dot_product(
                        Lp[i], Lp[i+1] - 1,       // i-th row minus diagonal
                        Lp[j], Lp[j+1] - 1,       // j-th row minus diagonal
                        Li, &new_Lx[0]);

                const double A_ij = Lx[k];

                if (j < i) {        // below diagonal?
                    const double L_jj = new_Lx[Lp[j+1] - 1];    // diagonal is last entry of j-th row
                    if(L_jj != 0){
                        new_Lx[k] = (A_ij - dp) / L_jj;
                    }
                } else if (j == i){    // on the diagonal?
                    if(A_ij - dp >= 0){
                        new_Lx[k] = std::sqrt(A_ij - dp);
                    }
                } else                // above diagonal -- input should be triangular!
                throw std::logic_error("Matrix passed in the wrong format - should be triangular");
            }
        }
        //        #ifndef NDEBUG
        //        for (int l = 0; l < n; ++l) {
        //            assert(Lx[Lp[l]-1] != 0);
        //            if(Lx[Lp[l]-1] == 0){
        //                finish_flag = false;
        //                break;
        //            }
        //        }
        //        #endif
        return true;
    }



} // namespace sym_lib