//
// Created by george on 2019-10-09.
//
#include <omp.h>
#include <vector>
#include <algorithm>
#include "sparse_blas_lib.h"
#include "../dense_blas/includes/BLAS.h"

namespace sym_lib
{
    //=========================== Left Looking SpTrSv ==========================
    void sptrsv_csr(int n, int *Lp, int *Li, double *Lx, double *x)
      {
//            int i, j;
            for (int i = 0; i < n; i++){
                  for (int j = Lp[i]; j < Lp[i + 1] - 1; j++){
                        x[i] -= Lx[j] * x[Li[j]]; //S2
                  }
                  x[i] /= Lx[Lp[i + 1] - 1]; //S1
            }
      }


    void sptrsv_csr_levelset(int n, const int *Lp, const int *Li, const double *Lx,
                               int levels, const int *levelPtr, const int *levelSet,
                             double *x)
      {
            #pragma omp parallel
            {
                  for (int l = 0; l < levels; l++)
                  {
                        #pragma omp for schedule(auto)
                        for (int k = levelPtr[l]; k < levelPtr[l + 1]; ++k)
                        {
                              int i = levelSet[k];
                              for (int j = Lp[i]; j < Lp[i + 1] - 1; j++)
                              {
                                    x[i] -= Lx[j] * x[Li[j]];//S1
                              }
                              x[i] /= Lx[Lp[i + 1] - 1]; //S2
                        }
                  }
            }
      }


   /*
    * It is left looking Sparse Triangular Solve
    * @param n Number of iterations or node
    * @param Lp the pointer array in CSC version
    * @param Li the index array in CSC version
    * @param Lx the value array in CSC version
    * @return x the output
    * @param levels number of levels in the DAg
    * @param levelPtr the pointer array in CSC format
    * that point to starting and ending point of nodes in a level
    * @param LevelSet the array that store nodes sorted based on their level
    * @param groupPtr the array pointer for groups
    * @param groupSet the array set for groups. Nodes are sorted based on their group
    */
    void sptrsv_csr_group_levelset(int *Lp, int *Li, double *Lx, double *x,
                                   int level_no, int *level_ptr, int *level_set,
                                   int *groupPtr, int *groupSet)
    {
#pragma omp parallel
        {
            for (int i1 = 0; i1 < level_no; ++i1)
            {
#pragma omp for schedule(auto)
                for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1)
                {
                    int group_idx = level_set[j1];
                    for (int k = groupPtr[group_idx]; k < groupPtr[group_idx + 1]; ++k)
                    {
                        int i = groupSet[k];
                        for (int j = Lp[i]; j < Lp[i + 1] - 1; j++)
                        {
                            x[i] -= Lx[j] * x[Li[j]];
                        }
                        x[i] /= Lx[Lp[i + 1] - 1];
                    }
                }
            }
        }
    }

    /*
    * It is left looking Sparse Triangular Solve
    * @param n Number of iterations or node
    * @param Lp the pointer array in CSC version
    * @param Li the index array in CSC version
    * @param Lx the value array in CSC version
    * @return x the output
    * @param level_no number of levels in the DAg
    * @param level_ptr the pointer array in CSC format
    * that point to starting and ending point of nodes in a level
    * @param par_ptr the array that point to the beginning and ending of a w-partition
    * @param partition the array that store nodes sorted based on their w-partition
    */
    void sptrsv_csr_lbc(int n, int *Lp, int *Li, double *Lx, double *x,
                          int level_no, int *level_ptr, int *par_ptr, int *partition){
            //======================================================
            #pragma omp parallel
            {
                  // iterate over l-partitions
                  for (int i1 = 0; i1 < level_no; ++i1)
                  {
                      // Iterate over all the w-partitions of a l-partition
                        #pragma omp for schedule(auto)
                        for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1)
                        {
                              // Iterate over all the node of a w-partition
                              for (int k1 = par_ptr[j1]; k1 < par_ptr[j1 + 1]; ++k1)
                              {
                                    //Detect the node
                                    int i = partition[k1];
                                    //Do the computation
                                    for (int j = Lp[i]; j < Lp[i + 1] - 1; j++)
                                    {
                                          x[i] -= Lx[j] * x[Li[j]];
                                    }
                                    x[i] /= Lx[Lp[i + 1] - 1];
                              }
                        }
                  }
            }
      }

    /*
    * It is left looking Sparse Triangular Solve
    * @param n Number of iterations or node
    * @param Lp the pointer array in CSC version
    * @param Li the index array in CSC version
    * @param Lx the value array in CSC version
    * @return x the output
    * @param level_no number of levels in the DAg
    * @param level_ptr the pointer array in CSC format
    * that point to starting and ending point of nodes in a level
    * @param par_ptr the array that point to the beginning and ending of a w-partition
    * @param partition the array that store nodes sorted based on their w-partition
    * @param groupPtr the array pointer for groups
    * @param groupSet the array set for groups. Nodes are sorted based on their group
    */
    void sptrsv_csr_group_lbc(int n, int *Lp, int *Li, double *Lx, double *x,
                              int level_no, int *level_ptr,
                              int *par_ptr, int *partition, int *groupPtr, int *groupSet){
        #pragma omp parallel
        {
            for (int i1 = 0; i1 < level_no; ++i1){
                #pragma omp for schedule(auto)
                for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1){
                    for (int k1 = par_ptr[j1]; k1 < par_ptr[j1 + 1]; ++k1){
                        int p = partition[k1];
                        for (int k = groupPtr[p]; k < groupPtr[p + 1]; ++k){
                            int i = groupSet[k];
                            for (int j = Lp[i]; j < Lp[i + 1] - 1; j++){
                                x[i] -= Lx[j] * x[Li[j]];
                            }
                            x[i] /= Lx[Lp[i + 1] - 1];
                        }
                    }
                }
            }
        }
    }







    void sptrsv_csr_lbc_buffer(int n, int *Lp, int *Li, double *Lx, double *x,
                   int level_no, int *level_ptr,
                   int *par_ptr, int *partition){
        //======================================================
        #pragma omp parallel
        {
            for (int i1 = 0; i1 < level_no; ++i1) // l-partition
            {
                #pragma omp for schedule(auto)
                for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1)
                {
                    for (int k1 = par_ptr[j1]; k1 < par_ptr[j1 + 1]; ++k1)
                    {
                        int i = partition[k1];
                        for (int j = Lp[k1]; j < Lp[k1 + 1] - 1; j++)
                        {
                            x[i] -= Lx[j] * x[Li[j]];
                        }
                        x[i] /= Lx[Lp[k1 + 1] - 1];
                    }
                }
            }
        }
    }

    void sptrsv_csr_lbc_double_buffer(int n, int *Lp, int *Li, double *Lx, double *x, double *x_copy,
                          int level_no, int *level_ptr,
                          int *par_ptr, int *partition){
    //======================================================
        #pragma omp parallel
        {
            for (int i1 = 0; i1 < level_no; ++i1) // l-partition
            {
                #pragma omp for schedule(auto)
                for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1)
                {
                    for (int k1 = par_ptr[j1]; k1 < par_ptr[j1 + 1]; ++k1)
                    {
                        for (int j = Lp[k1]; j < Lp[k1 + 1] - 1; j++)
                        {
                            x_copy[k1] -= Lx[j] * x_copy[Li[j]];
                        }
                        x_copy[k1] /= Lx[Lp[k1 + 1] - 1];
                    }
                }
            }
            #pragma omp for schedule(auto)
            for(int i = 0; i < n; i++){
                x[partition[i]] = x_copy[i];
            }
        }
    }

    void sptrsv_csr_w_sort_Hlevel_lbc(int n, int *Lp, int *Li, double *Lx, double *x,
                      int level_no, int *level_ptr,
                      int *par_ptr, int *partition){
        #pragma omp parallel
        {
              for (int i1 = 0; i1 < level_no; ++i1) // l-partition
              {
                    #pragma omp for schedule(auto)
                    for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1)
                    {
                          for (int k1 = par_ptr[j1]; k1 < par_ptr[j1 + 1]; ++k1)
                          {
                                int i = partition[k1];
                                for (int j = Lp[i]; j < Lp[i + 1] - 1; j++)
                                {
                                      x[i] -= Lx[j] * x[Li[j]];
                                }
                                x[i] /= Lx[Lp[i + 1] - 1];
                          }
                    }
              }
        };
    }



    //=========================== Right Looking SpTrSv ==========================
    /*
     * It is right looking Sparse Triangular Solve
     * @param n Number of iterations or node
     * @param Lp the pointer array in CSC version
     * @param Li the index array in CSC version
     * @param Lx the value array in CSC version
     * @return x the output
     */
    void sptrsv_csc(int n, int *Lp, int *Li, double *Lx, double *x)
    {
        for (int i = 0; i < n; i++)
        {
            x[i] /= Lx[Lp[i]]; //S1
            for (int j = Lp[i] + 1; j < Lp[i + 1]; j++)
            {
                x[Li[j]] -= Lx[j] * x[i]; //S2
            }
        }
    }

    /*
     * It is right looking Sparse Triangular Solve
     * @param n Number of iterations or node
     * @param Lp the pointer array in CSC version
     * @param Li the index array in CSC version
     * @param Lx the value array in CSC version
     * @return x the output
     * @param levels number of levels in the DAg
     * @param levelPtr the pointer array in CSC format
     * that point to starting and ending point of nodes in a level
     * @param LevelSet the array that store nodes sorted based on their level
     */
    void sptrsv_csc_levelset(int n, int *Lp, int *Li, double *Lx, double *x,
                             int levels, const int *levelPtr, const int *levelSet)
    {
        #pragma omp parallel
        {
            for (int l = 0; l < levels; l++)
            {
                #pragma omp for schedule(auto)
                for (int k = levelPtr[l]; k < levelPtr[l + 1]; ++k)
                {
                    int i = levelSet[k];
                    x[i] /= Lx[Lp[i]]; //S1
                    for (int j = Lp[i] + 1; j < Lp[i + 1]; j++)
                    {
                        #pragma omp atomic
                        x[Li[j]] -= Lx[j] * x[i]; //S2
                    }
                }
            }
        }
    }


    /*
     * It is left looking Sparse Triangular Solve
     * @param n Number of iterations or node
     * @param Lp the pointer array in CSC version
     * @param Li the index array in CSC version
     * @param Lx the value array in CSC version
     * @return x the output
     * @param levels number of levels in the DAg
     * @param levelPtr the pointer array in CSC format
     * that point to starting and ending point of nodes in a level
     * @param LevelSet the array that store nodes sorted based on their level
     * @param groupPtr the array pointer for groups
     * @param groupSet the array set for groups. Nodes are sorted based on their group
     */
    void sptrsv_csc_group_levelset(int *Lp, int *Li, double *Lx, double *x,
                                   int level_no, int *level_ptr, int *level_set,
                                   int *groupPtr, int *groupSet)
    {
        #pragma omp parallel
        {
            for (int i1 = 0; i1 < level_no; ++i1)
            {
                #pragma omp for schedule(auto)
                for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1)
                {
                    int group_idx = level_set[j1];
                    for (int k = groupPtr[group_idx]; k < groupPtr[group_idx + 1]; ++k)
                    {
                        int i = groupSet[k];
                        x[i] /= Lx[Lp[i]]; //S1
                        for (int j = Lp[i] + 1; j < Lp[i + 1]; j++)
                        {
                            #pragma omp atomic
                            x[Li[j]] -= Lx[j] * x[i]; //S2
                        }
                    }
                }
            }
        };
    }

    /*
     * It is right looking Sparse Triangular Solve
     * @param n Number of iterations or node
     * @param Lp the pointer array in CSC version
     * @param Li the index array in CSC version
     * @param Lx the value array in CSC version
     * @return x the output
     * @param level_no number of levels in the DAg
     * @param level_ptr the pointer array in CSC format
     * that point to starting and ending point of nodes in a level
     * @param par_ptr the array that point to the beginning and ending of a w-partition
     * @param partition the array that store nodes sorted based on their w-partition
     */
    void sptrsv_csc_lbc(int n, int *Lp, int *Li, double *Lx, double *x,
                        int level_no, int *level_ptr, int *par_ptr, int *partition)
    {

        #pragma omp parallel
        {
            // iterate over l-partitions
            for (int i1 = 0; i1 < level_no; ++i1)
            {
                // Iterate over all the w-partitions of a l-partition
                #pragma omp for schedule(auto)
                for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1)
                {
                    // Iterate over all the node of a w-partition
                    for (int k1 = par_ptr[j1]; k1 < par_ptr[j1 + 1]; ++k1)
                    {
                        //Detect the node
                        int i = partition[k1];
                        //Do the computation
                        x[i] /= Lx[Lp[i]]; //S1
                        for (int j = Lp[i] + 1; j < Lp[i + 1]; j++)
                        {
                            #pragma omp atomic
                            x[Li[j]] -= Lx[j] * x[i]; //S2
                        }
                    }
                }
            }
        }
    }


    /*
    * It is Right looking Sparse Triangular Solve
    * @param n Number of iterations or node
    * @param Lp the pointer array in CSC version
    * @param Li the index array in CSC version
    * @param Lx the value array in CSC version
    * @return x the output
    * @param level_no number of levels in the DAg
    * @param level_ptr the pointer array in CSC format
    * that point to starting and ending point of nodes in a level
    * @param par_ptr the array that point to the beginning and ending of a w-partition
    * @param partition the array that store nodes sorted based on their w-partition
    * @param groupPtr the array pointer for groups
    * @param groupSet the array set for groups. Nodes are sorted based on their group
    */
    void sptrsv_csc_group_lbc(int n, int *Lp, int *Li, double *Lx, double *x,
                              int level_no, int *level_ptr,
                              int *par_ptr, int *partition, int *groupPtr, int *groupSet)
    {
        #pragma omp parallel
        {
            for (int i1 = 0; i1 < level_no; ++i1)
            {
                #pragma omp for schedule(auto)
                for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1)
                {
                    for (int k1 = par_ptr[j1]; k1 < par_ptr[j1 + 1]; ++k1)
                    {
                        int p = partition[k1];

                        for (int k = groupPtr[p]; k < groupPtr[p + 1]; ++k)
                        {
                            int i = groupSet[k];
                            //Do the computation
                            x[i] /= Lx[Lp[i]]; //S1
                            for (int j = Lp[i] + 1; j < Lp[i + 1]; j++)
                            {
                                #pragma omp atomic
                                x[Li[j]] -= Lx[j] * x[i]; //S2
                            }
                        }
                    }
                }
            }
        };
    }


    void sptrsv_csc_block(int n, const int *Lp, const int *Li, const int *nrows,
                          double *Lx, double *x, int num_nodes,
                          const int *supernodes) {
        int i, p, k;
        auto *tempvec = new double[n]();

        for (i = 0; i < num_nodes; i++) {
            int super = supernodes[i];

            int width = supernodes[i + 1] - super;
            int nrow = nrows[i];

            custom_blas::lsolve_BLAS(nrow, width, &Lx[Lp[i]], &x[super]);
            custom_blas::matvec_BLAS(nrow, nrow - width, width, &Lx[Lp[i] + width], &x[super], tempvec);

            for (p = Lp[i] + width, k = 0; p < Lp[i] + nrow; p++, k++) {
                int idx = Li[p];
                x[idx] -= tempvec[k];
                tempvec[k] = 0;
            }
        }

        delete[]tempvec;
    }


    void sptrsv_csc_levelset_block(int n, const int *Lp, const int *Li,
                                   const int *nrows, double *Lx, double *x,
                                   const int *levelPtr,
                                   const int *levels, int n_lev,
                                   const int *supernodes, const int *sup2node,
                                   double **tempvecs) {
        int i, j, p, k;

        for (i = 0; i < n_lev; i++) {
            #pragma omp parallel default(shared) private(j, k, p)
            {
                int id = omp_get_thread_num();
                double *tempvec = tempvecs[id];

                #pragma omp for schedule(auto)
                for (j = levelPtr[i]; j < levelPtr[i + 1]; j++) {
                    int super = levels[j];

                    // sup2node[j] stores the index of current supernode in the block set
                    int index = sup2node[j];
                    int width = supernodes[index + 1] - super;
                    int nrow = nrows[index];

                    custom_blas::lsolve_BLAS(nrow, width, &Lx[Lp[index]], &x[super]);
                    custom_blas::matvec_BLAS(nrow, nrow - width, width, &Lx[Lp[index] + width], &x[super],
                                tempvec);

                    for (p = Lp[index] + width, k = 0; p < Lp[index] + nrow; p++, k++) {
                        int idx = Li[p];
                        #pragma omp atomic
                        x[idx] -= tempvec[k];
                        tempvec[k] = 0;
                    }
                }
            }
        }
    }


//    int blockedLsolve_mrhs(int n, size_t *Lp, int *Li, double *Lx,
//                           size_t *Li_ptr, int *col2sup, int *sup2col,
//                           int supNo, double *x, int n_rhs, int max_col){
//        int i, p, k;
//        double one[2], zero[2];
//        one[0] = 1.0;
//        one[1] = 0.;
//        zero[0] = 0.;
//        zero[1] = 0.;
//        int ione = 1;
//        int off_set = max_col < 0 ? n : max_col;
//        auto *tempvec = new double[n_rhs*off_set]();
//
//        for (i = 0; i < supNo; i++) {
//            int curCol = sup2col[i];
//            int nxtCol = sup2col[i + 1];
//            int supWdt = nxtCol - curCol;
//            int nSupR = Li_ptr[nxtCol] - Li_ptr[curCol];//row size of supernode
//            double *Ltrng = &Lx[Lp[curCol]];//first nnz of current supernode
//
//            SYM_DTRSM("L", "L", "N", "N", &supWdt,&n_rhs,one,Ltrng,
//                      &nSupR,&x[curCol],&n);
//
//            Ltrng = &Lx[Lp[curCol] + supWdt];//first nnz of below diagonal
//            int tmpRow = nSupR - supWdt;
//            SYM_DGEMM("N", "N", &tmpRow, &n_rhs, &supWdt, one, Ltrng,
//                      &nSupR,
//                      &x[curCol], &n,
//                      one, tempvec, &off_set);
//
//
//            for (int l = Li_ptr[curCol] + supWdt, k = 0; l < Li_ptr[nxtCol]; l++, k++) {
//                int idx = Li[l];
//                for (int j = 0; j < n_rhs; ++j) {
//                    auto tmp = k + j*off_set;
//                    x[idx+(j*n)] -= tempvec[tmp];
//                    tempvec[tmp] = 0;
//                }
//            }
//        }
//        delete[]tempvec;
//        return 1;
//    }
} // namespace sym_lib