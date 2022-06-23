//
// Created by labuser (Bangtian Liu) on 9/15/20.
//

#ifndef LBC_LIB_EXECUTOR_H
#define LBC_LIB_EXECUTOR_H
#include <algorithm>
#include <omp.h>
#include <cstring>
#include <functional>

namespace sym_lib
{
    int fs_csr_executor_sgroup(int n, int *Lp, int *Li, double *Lx, double *b, double *x, int *groupPtr, int *groupSet,  int ngroup,
                               int levels, int *levelPtr, int *levelSet);

    /*
     * @brief It is the executor that use levelset and grouping
     * @param n: number of nodes rows)
     * @param
     * @param
     */
    int CSRGeneralGroupingExecuter(int n, int *Lp, int *Li, double *Lx, double *x,
                                   int* group_range, int levels, int *levelPtr, int *levelSet);


    /*
     * @brief This is the grouping algorithm with levelset and serial and parallel region
     * @param Lp: pointer to the columns in CSR format
     * @param Li: columns indices in CSR format
     * @param Lx: Values in CSR format
     * @param x: the output x which is initialize by b
     * @param group_range: each node is a group of consecutive nodes which is shown by a range
     * @param levels: is the number of levels in levelset
     * @param levelPtr: is the pointer that shows the nodes belong to a level in levelset
     * @param Levelset: the nodes in each level divided by levelPtr value
     */
    int CSRGeneralGroupingExecutorWRegion(const int *Lp, const int *Li, const double *Lx,
                                          double *x,
                                          const int* group_range,
                                          int levels,
                                          const int *levelPtr, const int *levelSet);


    /*
     * @brief
     * @param
     * @param
     * @param
     */
    int CSRLevelExecuter(int n, int *Lp, int *Li, double *Lx, double *x,
                         int levels, int *levelPtr, int *levelSet);

    /**
     * @brief profile serial implementation of triangualr solver based on CSR format
     * @param n number of rows
     * @param Lp pointers to starting address of one row
     * @param Li index array of non-zeros
     * @param flops number of floating operations for serial code
     * @param access_nnz  the total memory access operations on sparse matrix and x array
     * @param reuse_nnz   the reused memory access across two consecutive iterations
     */
    void fs_csr_stat(int n, int *Lp, int *Li,  int &flops, int &access_nnz, int &reuse_nnz);

    /**
     * @brief profile serial implementation of incompelete choleksy factorization based on CSC format
     * @param n
     * @param val
     * @param colPtr
     * @param rowIdx
     * @param flops
     * @param access_nnz
     * @param reuse_nnz
     */
    void ic0_csc_stat(int n, int * colPtr, int *rowIdx, int &flops, int &access_nnz, int &reuse_nnz);

    /**
     * brief  profile parallel implementation of triangular solver based on CSR
     * @param Lp row pointer in the CSR format
     * @param Li index array in the CSR format
     * @param groupPtr Pointer to the starting location of one group
     * @param groupSet Pointer to the column indices in one group
     * @param levels number of levels
     * @param levelPtr Pointer to the starting location of one level
     * @param levelSet Pointer to index array of one level
     * @param lcost store the maximum difference between nodes for each level
     */
    void fs_csr_levelset_stat(int *Lp, int *Li, int *groupPtr, int *groupSet,
                              int levels, int *levelPtr, int *levelSet, int *lcost);
    /**
     * @brief similar to previous function, but works for coarsening level method
     * @param n
     * @param Lp
     * @param Li
     * @param level_no
     * @param level_ptr
     * @param par_ptr
     * @param partition
     * @param lcost
     */
    void sptrsv_csr_lbc_stat(int n, int *Lp, int *Li,
                             int level_no, int *level_ptr,
                             int *par_ptr, int *partition, int *lcost);

    /**
     * @brief profiling executor code based on coarsening and grouping method
     * @param n
     * @param Lp
     * @param Li
     * @param level_no
     * @param level_ptr
     * @param par_ptr
     * @param partition
     * @param groupPtr
     * @param groupSet
     * @param lcost
     */
    void sptrsv_csr_group_lbc_stat(int n, int *Lp, int *Li,
                                   int level_no, int *level_ptr,
                                   int *par_ptr, int *partition, int *groupPtr, int *groupSet, int *lcost);


}

#endif //LBC_LIB_EXECUTOR_H
