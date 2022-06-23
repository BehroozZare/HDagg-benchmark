//
// Created by kazem on 10/9/19.
//

#ifndef FUSION_SPARSEBLASLIB_H
#define FUSION_SPARSEBLASLIB_H

#include "def.h"

namespace sym_lib {

    //============================= SPTRSV ============================
    ///\description It is left looking Sparse Triangular Solve
    ///\input n Number of iterations or node
    ///\input Lp the pointer array in CSC version
    ///\input Li the index array in CSC version
    ///\input Lx the value array in CSC version
    ///\inout x the output
    void sptrsv_csr(int n, int *Lp, int *Li, double *Lx, double *x);


    ///\Description It is a parallel left looking Sparse Triangular Solve which uses wavefronts schedule
    ///\input n Number of iterations or node
    ///\input Lp the pointer array in CSC version
    ///\input Li the index array in CSC version
    ///\input Lx the value array in CSC version
    ///\input levels number of levels in the DAg
    ///\input levelPtr the pointer array in CSC format
    /// that point to starting and ending point of nodes in a level
    ///\input LevelSet the array that store nodes sorted based on their level
    ///\inout x the output
    void sptrsv_csr_levelset(int n, const int *Lp, const int *Li, const double *Lx,
                        int levels, const int *levelPtr, const int *levelSet,
                        double *x);

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
                                   int level_no, int *level_ptr,
                                   int *levelSet, int *groupPtr, int *groupSet);

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
    void
    sptrsv_csr_lbc(int n, int *Lp, int *Li, double *Lx, double *x,
                   int level_no, int *level_ptr,
                   int *par_ptr, int *partition);


    /**
    * @brief perform SpTrsv_CSR based on grouping and LBC
    * @param n number of rows
    * @param Lp Row pointer in the CSR format
    * @param Li Index array in the CSR format
    * @param Lx Val array in the CSR format
    * @param x  right hand side
    * @param level_no number of levels
    * @param level_ptr  pointer to one coarsen level
    * @param par_ptr  pointer to one partition
    * @param partition pointer to array storing node Index
    * @param groupPtr Pointer to the starting location of one group
    * @param groupSet Pointer to the column indices in one group
    */
    void sptrsv_csr_group_lbc(int n, int *Lp, int *Li, double *Lx, double *x,
                              int level_no, int *level_ptr,
                              int *par_ptr, int *partition, int *groupPtr, int *groupSet);

    void
    sptrsv_csr_lbc_buffer(int n, int *Lp, int *Li, double *Lx, double *x,
                          int level_no, int *level_ptr,
                          int *par_ptr, int *partition);

    void
    sptrsv_csr_lbc_double_buffer(int n, int *Lp, int *Li, double *Lx, double *x, double *x_copy,
                                 int level_no, int *level_ptr,
                                 int *par_ptr, int *partition);


    void sptrsv_csr_lbc_seq(int n, int *Lp, int *Li, double *Lx, double *x,
                            int level_no, int *level_ptr,
                            int *par_ptr, int *partition);


    /**
     * @brief perform SpTrsv_CSR based on greedymerging and LBC
     * @param n number of rows
     * @param Lp Row pointer in the CSR format
     * @param Li Index array in the CSR format
     * @param Lx Val array in the CSR format
     * @param x  right hand side
     * @param level_no number of levels
     * @param level_ptr  pointer to one coarsen level
     * @param par_ptr  pointer to one partition
     * @param partition pointer to array storing node Index
    */
    void sptrsv_csr_w_sort_Hlevel_lbc(int n, int *Lp, int *Li, double *Lx, double *x,
                                      int level_no, int *level_ptr,
                                      int *par_ptr, int *partition);

    //==================================== SpTrSv Right Looking ====================================
    /*
     * @brief Right Looking Sparse Triangular Solver using CSC format
     * @param n number of iterations or nodes of the dependency DAGs
     * @param Lp the pointer array in CSC format
     * @param Li the index array in CSC format
     * @param Lx the value array in CSC format (it is the output in here)
     * @return x The result is stored in vector x
     */
    void sptrsv_csc(int n, int *Lp, int *Li, double *Lx, double *x);

    /*
     * @brief Right Looking Sparse Triangular Solver using CSC format
     * @param n number of iterations or nodes of the dependency DAGs
     * @param Lp the pointer array in CSC format
     * @param Li the index array in CSC format
     * @param Lx the value array in CSC format (it is the output in here)
     * @param levels: The number of levels in the levelset schedule
     * @param LevelPtr The pointer that shows the starting and ending point of nodes in each level
     * @param LevelSet The nodes sorted based on the level in this array
     * @return x The result is stored in vector x
     */
    void sptrsv_csc_levelset(int n, int *Lp, int *Li, double *Lx, double *x,
                             int levels, const int *levelPtr, const int *levelSet);

    /*
     * @brief Right Looking Sparse Triangular Solver using CSC format
     * @param n: Number of nodes
     * @param Lp: The pointer in CSC format
     * @param Li: The index array in CSC format
     * @param Lx: The value array in CSC format
     * @param nLevels: number of levels in LBC scheduler
     * @param levelPtr: the level ptr (like CSC format) (shows the start and end of the set of w-partitions in a lavel)
     * @param ParPtr: The pointer that shows the nodes inside each w-partition (think of it as two CSC on top of each other)
     * @param Partitions: the nodes sorted based on the w-partitions that they are in
     * @return x: the output vector which is initialized by b first
     */
    void sptrsv_csc_lbc(int n, int *Lp, int *Li, double *Lx, double *x,
                        int level_no, int *level_ptr, int *par_ptr, int *partition);

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
                                   int *groupPtr, int *groupSet);

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
                              int *par_ptr, int *partition, int *groupPtr, int *groupSet);


    /// \param n number of columns/rows
    /// \param Lp pointer array in CSC format
    /// \param Li index array in CSC format
    /// \param nrows
    /// \param Lx value array in CSC format
    /// \param x the vector to solve for initialized with b
    /// \param num_nodes
    /// \param supernodes
    void sptrsv_csc_block(int n, const int *Lp, const int *Li, const int *nrows,
                          double *Lx, double *x, int num_nodes,
                          const int *supernodes);

    /// \param n number of columns/rows
    /// \param Lp pointer array in CSC format
    /// \param Li index array in CSC format
    /// \param nrows
    /// \param Lx value array in CSC format
    /// \param x the vector to solve for initialized with b
    /// \param num_nodes
    /// \param supernodes
    void sptrsv_csc_levelset_block(int n, const int *Lp, const int *Li,
                                   const int *nrows, double *Lx, double *x,
                                   const int *levelPtr,
                                   const int *levels, int n_lev,
                                   const int *supernodes, const int *sup2node,
                                   double **tempvecs);


    /// \param n number of columns/rows
    /// \param Lp pointer array in CSC format
    /// \param Li index array in CSC format
    /// \param nrows
    /// \param Lx value array in CSC format
    /// \param x the vector to solve for initialized with b
    /// \param num_nodes
    /// \param supernodes
    int blockedLsolve_mrhs(int n, size_t *Lp, int *Li, double *Lx, size_t *Li_ptr, int *col2sup, int *sup2col,
                           int supNo, double *x, int n_rhs, int max_col);


    //============================================== ICh0 ===========================================
    //=================================== Right Looking Incomplete Cholesky ===============================
    /*
     * @brief Right looking incomplete cholesky
     * @param n number of iterations or nodes of the dependency DAGs
     * @param Lp the pointer array in CSC format
     * @param Li the index array in CSC format
     * @return Lx the value array in CSC format (it is the output in here)
     */
    bool spic0_csc_RL_serial(int n, int *Lp, int *Li, double *Lx);

    /*
     * @brief Right looking incomplete cholesky
     * @param n number of iterations or nodes of the dependency DAGs
     * @param Lp the pointer array in CSC format
     * @param Li the index array in CSC format
     * @param Lx the value array in CSC format (it is the output in here)
     * @param nLevels is the number of wavefront using levelset scheduling
     * @param LevelPtr The pointer that shows the starting and ending point of nodes in each level
     * @param LevelSet The nodes sorted based on the level in this array
     */
    void spico_csc_RL_levelset(int n, const int *Lp, const int *Li, double *Lx,
                               int levels, const int *levelPtr, const int *levelSet);

    /*
    * @brief Right looking incomplete cholesky
    * @param n number of iterations or nodes of the dependency DAGs
    * @param Lp the pointer array in CSC format
    * @param Li the index array in CSC format
    * @param Lx the value array in CSC format (it is the output in here)
    * @param nLevels is the number of wavefront using levelset scheduling
    * @param LevelPtr The pointer that shows the starting and ending point of nodes in each level
    * @param LevelSet The nodes sorted based on the level in this array
    * @param groupPtr the array pointer for groups
    * @param groupSet the array set for groups. Nodes are sorted based on their group
    */
    void spico_csc_RL_group_levelset(int n, const int *Lp, const int *Li, double *Lx,
                                     int levels, const int *levelPtr,
                                     const int *levelSet, int *groupPtr, int *groupSet);

    /*
     * @brief Right looking incomplete cholesky
     * @param n: Number of nodes
     * @param Lp: The pointer in CSC format
     * @param Li: The index array in CSC format
     * @param Lx: The value array in CSC format
     * @param nLevels: number of levels in LBC scheduler
     * @param level_ptr: the level ptr (like CSC format) (shows the start and end of the set of w-partitions in a lavel)
     * @param par_ptr: The pointer that shows the nodes inside each w-partition (think of it as two CSC on top of each other)
     * @param partition: the nodes sorted based on the w-partitions that they are in
     */
    void
    spico_csc_RL_lbc(int n, int *Lp, int *Li, double *Lx, int level_no, int *level_ptr, int *par_ptr, int *partition);

    /*
     * @brief Right looking incomplete cholesky
     * @param n: Number of nodes
     * @param Lp: The pointer in CSC format
     * @param Li: The index array in CSC format
     * @param Lx: The value array in CSC format
     * @param nLevels: number of levels in LBC scheduler
     * @param level_ptr: the level ptr (like CSC format) (shows the start and end of the set of w-partitions in a lavel)
     * @param par_ptr: The pointer that shows the nodes inside each w-partition (think of it as two CSC on top of each other)
     * @param partition: the nodes sorted based on the w-partitions that they are in
     * @param groupPtr the array pointer for groups
     * @param groupSet the array set for groups. Nodes are sorted based on their group
     */
    void spico_csc_RL_group_lbc(int n, int *Lp, int *Li, double *Lx,
                                int level_no, int *level_ptr,
                                int *par_ptr, int *partition, int *groupPtr, int *groupSet);


    void spico_csc_lbc_reordered(int n, double *val, int *rowPtr, int *colIdx,
                                 int level_no, int *level_ptr,
                                 int *par_ptr, int *partition);


    void spico_csc_group_levelset_v2(int n, const int *Lp, const int *Li, double *Lx,
                                     int levels, const int *levelPtr,
                                     const int *levelSet, int *groupPtr, int *groupSet);


    void spico_csc_lbc(int n, double *Lx, int *Lp, int *Li,
                       int level_no, int *level_ptr,
                       int *par_ptr, int *partition);




    //=================================== Left Looking Incomplete Cholesky ===============================
    /*
     * @brief left looking incomplete cholesky
     * @param n number of iterations or nodes of the dependency DAGs
     * @param Lp the pointer array in CSC format
     * @param Li the index array in CSC format
     * @param Lx the value array in CSC format (it is the output in here)
     * @param PrunePtr In incomplete cholesky it is the
     * pointer in CSR format of row j of the input matrix
     * @param PrunePtr In incomplete cholesky it is the
     * index set in CSR format of row j of the input matrix
     */
    bool spic0_csc_LL_serial(int n, const int *Lp, const int *Li, double *Lx,
                             const int *PrunePtr, const int *PruneSet);


    /*
     * @brief left looking incomplete cholesky
     * @param n number of iterations or nodes of the dependency DAGs
     * @param Lp the pointer array in CSC format
     * @param Li the index array in CSC format
     * @param Lx the value array in CSC format (it is the output in here)
     * @param PrunePtr In incomplete cholesky it is the
     * pointer in CSR format of row j of the input matrix
     * @param PrunePtr In incomplete cholesky it is the
     * index set in CSR format of row j of the input matrix
     */
    bool spic0_csc_LL_serial_test(int n, const int *Lp, const int *Li, double *Lx,
                                  const int *PrunePtr, const int *PruneSet);

    /*
     * @brief left looking incomplete cholesky
     * @param n number of iterations or nodes of the dependency DAGs
     * @param Lp the pointer array in CSC format
     * @param Li the index array in CSC format
     * @param Lx the value array in CSC format (it is the output in here)
     * @param PrunePtr In incomplete cholesky it is the
     * pointer in CSR format of row j of the input matrix
     * @param PrunePtr In incomplete cholesky it is the
     * index set in CSR format of row j of the input matrix
     * @param nLevels is the number of wavefront using levelset scheduling
     * @param LevelPtr The pointer that shows the starting and ending point of nodes in each level
     * @param LevelSet The nodes sorted based on the level in this array
     * @param groupPtr the array pointer for groups
     * @param groupSet the array set for groups. Nodes are sorted based on their group
     */
    bool spic0_csc_LL_group_Levelset(int n, const int *Lp, const int *Li, double *Lx,
                                     const int *PrunePtr, const int *PruneSet,
                                     int levels, const int *levelPtr, const int *levelSet,
                                     int *groupPtr, int *groupSet);

    /*
     * @brief left looking incomplete cholesky
     * @param n number of iterations or nodes of the dependency DAGs
     * @param Lp the pointer array in CSC format
     * @param Li the index array in CSC format
     * @param Lx the value array in CSC format (it is the output in here)
     * @param PrunePtr In incomplete cholesky it is the
     * pointer in CSR format of row j of the input matrix
     * @param PrunePtr In incomplete cholesky it is the
     * index set in CSR format of row j of the input matrix
     * @param nLevels is the number of wavefront using levelset scheduling
     * @param LevelPtr The pointer that shows the starting and ending point of nodes in each level
     * @param LevelSet The nodes sorted based on the level in this array
     */
    bool spic0_csc_LL_Levelset(int n, const int *Lp, const int *Li, double *Lx,
                               const int *PrunePtr, const int *PruneSet,
                               int nLevels, const int *LevelPtr, const int *LevelSet);


    /*
     * @brief left looking incomplete cholesky
     * @param n: Number of nodes
     * @param Ap: The pointer in CSC format
     * @param Ai: The index array in CSC format
     * @param Ax: The value array in CSC format
     * @param prunePtr: the prune pointer (like CSC format)
     * @param pruneSet: the prune set
     * @param nLevels: number of levels in LBC scheduler
     * @param levelPtr: the level ptr (like CSC format) (shows the start and end of the set of w-partitions in a lavel)
     * @param ParPtr: The pointer that shows the nodes inside each w-partition (think of it as two CSC on top of each other)
     * @param Partitions: the nodes sorted based on the w-partitions that they are in
     */
    bool spic0_csc_LL_LBC(int n, const int *Lp, const int *Li, double *Lx,
                          const int *PrunePtr, const int *PruneSet,
                          int nLevels, const int *LevelPtr, const int *ParPtr, const int *Partitions);

    /*
     * @brief left looking incomplete cholesky
     * @param n: Number of nodes
     * @param Ap: The pointer in CSC format
     * @param Ai: The index array in CSC format
     * @param Ax: The value array in CSC format
     * @param prunePtr: the prune pointer (like CSC format)
     * @param pruneSet: the prune set
     * @param nLevels: number of levels in LBC scheduler
     * @param levelPtr: the level ptr (like CSC format) (shows the start and end of the set of w-partitions in a lavel)
     * @param ParPtr: The pointer that shows the nodes inside each w-partition (think of it as two CSC on top of each other)
     * @param Partitions: the nodes sorted based on the w-partitions that they are in
     * @param groupPtr the array pointer for groups
     * @param groupSet the array set for groups. Nodes are sorted based on their group
     */
    bool spic0_csc_LL_group_LBC(int n, const int *Lp, const int *Li, double *Lx,
                                const int *PrunePtr, const int *PruneSet,
                                int nLevels, const int *LevelPtr, const int *ParPtr, const int *Partitions,
                                int *groupPtr, int *groupSet);


    ///\description: This function do a dot product on two rows in CSR format
    ///\input l1: Is the beginning index of row 1
    ///\input u1: Is the end index of row 1
    ///\input l2: Is the beginning index of row 2
    ///\input u2: Is the end index of row 2
    ///\input indices: Is the index array (stores column indices of row 1 and row 2)
    ///\input data: Is the current partial factorized matrix
    ///\output The dot product of the two rows EXCLUDING diagonal index of the smaller row (the one with smaller index)
    inline double sparse_dot_product(int l1, int u1, int l2, int u2, const int *indices, const double *data);


    ///\description: This function perform an incomplete up-looking cholesky factorization
    ///\input n: number of columns/rows
    ///\input Lp: pointer array in CSR format
    ///\input Li: Index array in CSR format
    ///\input Lx: Value array in CSR format
    ///\input new_data: The output buffer filled with zero
    ///\output a flag that show whether something is wrong or not
    bool spic0_csr_uL_serial(int n, const int *Lp, const int *Li, double *Lx, double *new_data);


    ///\description: This function perform an incomplete up-looking cholesky factorization with levelset scheduling
    ///\input n: number of columns/rows
    ///\input Lp: pointer array in CSR format
    ///\input Li: Index array in CSR format
    ///\input Lx: Value array in CSR format
    ///\input new_data: The output buffer filled with zero
    ///\output is the pointer same as new_data
    bool spic0_csr_UL_Levelset(int n, const int *Lp, const int *Li, double *Lx, double *new_data,
                               int nLevels, const int *LevelPtr, const int *LevelSet);


    ///\description: This function perform an incomplete up-looking cholesky factorization with levelset scheduling
    ///\input n: number of columns/rows
    ///\input Lp: pointer array in CSR format
    ///\input Li: Index array in CSR format
    ///\input Lx: Value array in CSR format
    ///\input new_data: The output buffer filled with zero
    ///\output is the pointer same as new_data
    ///\input nLevels: Index array in CSR format
    ///\input LevelPtr: Value array in CSR format
    ///\input LevelSet: The output buffer filled with zero
    ///\input group_ptr: The group_ptr array that define groups in group_set array
    ///\input group_set: Array with sorted nodes based on their group
    bool spic0_csr_UL_group_Levelset(int n, const int *Lp, const int *Li, double *Lx, double *new_data,
                                     int nLevels, const int *LevelPtr, const int *LevelSet,
                                     const int *group_ptr, const int *group_set);

    ///\description: This function perform an incomplete up-looking cholesky factorization with LC scheduling
    ///\input n: number of columns/rows
    ///\input Lp: pointer array in CSR format
    ///\input Li: Index array in CSR format
    ///\input Lx: Value array in CSR format
    ///\input new_data: The output buffer filled with zero
    ///\output is the pointer same as new_data
    bool spic0_csr_UL_LBC(int n, const int *Lp, const int *Li, double *Lx, double *new_data,
                          int nLevels, const int *LevelPtr, const int *ParPtr, const int *Partitions);


    //===========================================================================================
    //                                 Sparse Incomplete LU0
    //===========================================================================================


    ///\Description: This function reqrrange the data so that column indices are in ascending order
    ///\Input int N, the order of the system.
    ///\Input int NZ_NUM, the number of nonzeros.
    ///\Input int IA[N+1], the compressed row index vector.
    ///\Input_output int JA[NZ_NUM], the column indices of the matrix values.
    /// the order of the entries of JA may have changed because of
    /// the sorting.
    ///\Input_output double A[NZ_NUM], the matrix values.  On output, the
    /// order of the entries may have changed because of the sorting.
    void rearrange_cr(int n, int nz_num, int ia[], int ja[], double a[]);

    ///\Description: The Serial algorithm for incomplete LU decomposition with zero fill-in
    ///\Input int N, the order of the system.
    ///\Input int NZ_NUM, the number of nonzeros.
    ///\Input int Ap[N+1], Ai[NZ_NUM], the row and column indices
    /// of the matrix values.  The row vector has been compressed.
    ///\Input double Ax[NZ_NUM], the matrix values.
    ///\Input int A_diag[N], the index of the diagonal element of each row.
    ///\Output double l[NZ_NUM], the ILU factorization of A.
    void
    spilu0_ul_csr_serial(int n, int nnz, const int *Ap, const int *Ai, const double *Ax, const int *A_diag, double *l);


    ///\Description: The Levelset algorithm for incomplete LU decomposition with zero fill-in
    ///\Input int N, the order of the system.
    ///\Input int NZ_NUM, the number of nonzeros.
    ///\Input int Ap[N+1], Ji[NZ_NUM], the row and column indices
    /// of the matrix values.  The row vector has been compressed.
    ///\Input double Ax[NZ_NUM], the matrix values.
    ///\Input int A_diag[N], the index of the diagonal element of each row.
    ///\Output double l[NZ_NUM], the ILU factorization of A.
    ///\Input nlevels Number of levels in levelset schedule
    ///\Input LevelPtr Pointer to Levelset array that shows nodes inside each level
    ///\Input Levelset Nodes' id sorted based on their level
    ///\Input temp: A variable used in the executor. Should have the same size of the full value array
    void spilu0_ul_csr_levelset(int n, int nz_num, const int *Ap, const int *Ai,
                                const double *Ax, const int *A_diag, double *l,
                                int nLevels, const int *LevelPtr, const int *LevelSet, int *temp);


    ///\Description: The Levelset algorithm for incomplete LU decomposition with zero fill-in
    ///\Input int N, the order of the system.
    ///\Input int NZ_NUM, the number of nonzeros.
    ///\Input int Ap[N+1], Ji[NZ_NUM], the row and column indices
    /// of the matrix values.  The row vector has been compressed.
    ///\Input double Ax[NZ_NUM], the matrix values.
    ///\Input int A_diag[N], the index of the diagonal element of each row.
    ///\Output double l[NZ_NUM], the ILU factorization of A.
    ///\Input nlevels Number of levels in levelset schedule
    ///\Input LevelPtr Pointer to Levelset array that shows nodes inside each level
    ///\Input Levelset Nodes' id sorted based on their level
    ///\Input LevelPtr Pointer to Levelset array that shows nodes inside each level
    ///\Input Levelset Nodes' id sorted based on their level
    ///\Input group_ptr pointer to group_set array
    ///\Input group_set vertices
    ///\Input temp: A variable used in the executor. Should have the same size of the full value array
    void spilu0_ul_csr_group_levelset(int n, int nnz, const int *Ap, const int *Ai,
                                      const double *Ax, const int *A_diag, double *l,
                                      int nLevels, const int *LevelPtr, const int *LevelSet,
                                      const int *group_ptr, const int *group_set, int *temp);

    ///\Description: The Serial algorithm for incomplete LU decomposition with zero fill-in
    ///\Input int N, the order of the system.
    ///\Input int NZ_NUM, the number of nonzeros.
    ///\Input int Ap[N+1], Ai[NZ_NUM], the row and column indices
    /// of the matrix values.  The row vector has been compressed.
    ///\Input double Ax[NZ_NUM], the matrix values.
    ///\Input int A_diag[N], the index of the diagonal element of each row.
    ///\Output double l[NZ_NUM], the ILU factorization of A.
    ///\Input level_no Number of levels in LBC schedule
    ///\Input level_ptr Pointer to par_ptr array where w-partitions are defined
    ///\Input par_ptr pointer to node_ptr array where nodes for each part are defined
    ///\Input partition Nodes' id sorted based on their partition
    void spilu0_ul_csr_lbc(int n, int nnz, const int *Ap, const int *Ai,
                           const double *Ax, const int *A_diag, double *l,
                           int level_no, const int *level_ptr, const int *par_ptr, const int *partition, int *temp);


    //===========================================================================================
    //                                 Sparse Incomplete IC with K fill-ins
    //===========================================================================================
    bool spick_csr_UL_serial(int n, const int *Lp, const int *Li, double *Lx,
                             double *new_Lp, double *new_Li, double *new_Lx, double &new_nnz);


} // namespace sym_lib

#endif //FUSION_SPARSEBLASLIB_H
