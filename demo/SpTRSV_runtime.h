//
// Created by behrooz on 1/17/21.
//

#ifndef SPTRSV_FINAL_H
#define SPTRSV_FINAL_H


#define ANALYSE_STATISTIC false

#include <sparse_inspector.h>
#include <lbc.h>
#include <glc.h>
#include <Group.h>
#include <Utils.h>
#include <executor.h>
#include <StatSpMat.h>

#include <algorithm>
#include <vector>
#include <cmath>
#include <numeric>
#include <list>
#include <mkl.h>

#include "FusionDemo.h"
#include "sparse_blas_lib.h"
#include "SuperNodalTools.h"
#include "BCSCMatrix.h"
#include "BLAS.h"

#ifdef DAGP

#include "dagp_utils.h"
#include <dgraphReader.h>
#include <info.h>
#include <rvcycle.h>

#endif

#ifdef SPMP

#include "SpMP_Interface.h"

#endif


namespace sym_lib
{
    /*
     * @brief Computing Data Dependency Graph of the kernel and also apply grouping to the DAG
     * @param A Matrix A in CSR version
     * @return DAG_ptr The pointer array in CSC format
     * @return DAG_set The index array in CSC format. It stores the child of each node
     * @return groupPtr The pointer array in CSC format
     * @return groupSet The index array in CSC format. It stores the group's nodes
     * @return groupInv For each node, we can determine in which group it resides
     * @return ngroup number of groups
    */
    timing_measurement computingGroupedDAG_CSR(const CSR* A, std::vector<int>& DAG_ptr, std::vector<int>& DAG_set,
                                               std::vector<int>& groupPtr, std::vector<int>& groupSet, int& ngroup){
        timing_measurement  DAG_time;
        DAG_time.start_timer();

        int n = A->n;
        int nnz = A->nnz;
        groupPtr.resize(n + 1, 0);
        groupSet.resize(n, 0);
        std::vector<int> groupInv(n, 0);

        group g(A->n, A->p, A->i);
        g.inspection_sptrsvcsr_v1(groupPtr.data(), groupSet.data(), ngroup, groupInv.data());

        std::vector<std::vector<int>> DAG;
        DAG.resize(ngroup);

        fs_csr_inspector_dep(ngroup, groupPtr.data(), groupSet.data(), groupInv.data(), A->p, A->i, DAG);

        size_t count=0;
        for(int j = 0; j < DAG.size(); ++j) {
            DAG[j].erase(std::unique(DAG[j].begin(), DAG[j].end()), DAG[j].end());
            count+=DAG[j].size();
        }
        DAG_ptr.resize(n + 1, 0);
        DAG_set.resize(count + n, 0);

        long int cti,edges=0;
        for(cti = 0, edges = 0; cti < ngroup; cti++){
            DAG_ptr[cti] = edges;
            DAG_set[edges++] = cti;
            for(int ctj = 0; ctj < DAG[cti].size(); ctj++) {
                DAG_set[edges++] = DAG[cti][ctj];
            }
        }
        DAG_ptr[cti] = edges;
        DAG_time.measure_elapsed_time();

        std::vector<bool> check(n, false);
        return DAG_time;
    }

    /*
     * @brief Computing Data Dependency Graph of the kernel and also apply grouping to the DAG
     * @param A Matrix A in CSR version
     * @return DAG_ptr The pointer array in CSC format
     * @return DAG_set The index array in CSC format. It stores the child of each node
     * @return groupPtr The pointer array in CSC format
     * @return groupSet The index array in CSC format. It stores the group's nodes
     * @return groupInv For each node, we can determine in which group it resides
     * @return ngroup number of groups
    */
    timing_measurement computingGroupedDAG_CSR(int n, std::vector<int>& CSR_DAG_ptr, std::vector<int>& CSR_DAG_set,
                                               std::vector<int>& DAG_ptr, std::vector<int>& DAG_set,
                                               std::vector<int>& groupPtr, std::vector<int>& groupSet, int& ngroup){
        timing_measurement  DAG_time;
        DAG_time.start_timer();

        groupPtr.resize(n + 1, 0);
        groupSet.resize(n, 0);
        std::vector<int> groupInv(n, 0);

        group g(n, CSR_DAG_ptr.data(), CSR_DAG_set.data());
        g.inspection_sptrsvcsr_v1(groupPtr.data(), groupSet.data(), ngroup, groupInv.data());

        std::vector<std::vector<int>> DAG;
        DAG.resize(ngroup);

        fs_csr_inspector_dep(ngroup, groupPtr.data(), groupSet.data(),
                             groupInv.data(), CSR_DAG_ptr.data(), CSR_DAG_set.data(), DAG);

        size_t count=0;
        for(int j = 0; j < DAG.size(); ++j) {
            DAG[j].erase(std::unique(DAG[j].begin(), DAG[j].end()), DAG[j].end());
            count+=DAG[j].size();
        }
        DAG_ptr.resize(n + 1, 0);
        DAG_set.resize(count + n, 0);

        long int cti,edges=0;
        for(cti = 0, edges = 0; cti < ngroup; cti++){
            DAG_ptr[cti] = edges;
            DAG_set[edges++] = cti;
            for(int ctj = 0; ctj < DAG[cti].size(); ctj++) {
                DAG_set[edges++] = DAG[cti][ctj];
            }
        }
        DAG_ptr[cti] = edges;
        DAG_time.measure_elapsed_time();

        std::vector<bool> check(n, false);
        return DAG_time;
    }

    /*
     * @brief Computing Data Dependency Graph of the kernel and store it in CSC format
     * @param A Matrix A in CSR version
     * @return DAG_ptr The pointer array in CSC format
     * @return DAG_set The index array in CSC format. It stores the child of each node
    */
    timing_measurement computingDAG_CSR(const CSR* A, std::vector<int>& DAG_ptr, std::vector<int>& DAG_set){
        timing_measurement  DAG_time;
        int n_ = A->n;
        DAG_time.start_timer();

        std::vector<std::vector<int>> DAG;
        DAG.resize(A->n);
        fs_csr_inspector_dep(A->n, A->p, A->i, DAG);
        size_t count = 0;
        for (auto & j : DAG){
            j.erase(std::unique(j.begin(), j.end()), j.end());
            count += j.size();
        }

        DAG_ptr.reserve(n_ + 1);
        DAG_set.resize(count + n_);

        long int cti, edges = 0;
        for (cti = 0, edges = 0; cti < n_; cti++)
        {
            // for each edge (cti, ctj), save the starting point for ctjs of cti
            DAG_ptr.push_back(edges);
            // Add the (cti,cti) edge because it is not inside the DAG
            DAG_set[edges] = cti;
            edges++;
            // Iterate over all the dependent ctjs that is cti -> ctj
            for(int ctj = 0; ctj < DAG[cti].size(); ctj++)
            {
                // Add the edge into the graph_edge_idx
                DAG_set[edges] = DAG[cti][ctj];
                edges++;
            }
        }
        // This is the number of all the edges, just like CSR or CSC format
        DAG_ptr.push_back(edges);
        DAG_time.measure_elapsed_time();
        return DAG_time;
    }

    /*
     * @brief Computing Data Dependency Graph of the kernel and store it in CSC format
     * @param n Number of Nodes in the DAG
     * @param DAG_ptr The pointer array in CSC format
     * @param DAG_set The index array in CSC format. It stores the child of each node
     * @return LevelPtr the pointer array in CSC format that store pointers to the nodes inside a level
     * @return LevelSet the nodes are sorted based on level in this array
    */
    timing_measurement computingLevelSet_CSC(int n, const int* DAG_ptr, const int* DAG_set,
                                             std::vector<int>& LevelPtr, std::vector<int>& LevelSet, int& nlevels){
        timing_measurement LevelSet_time;
        LevelSet_time.start_timer();
        LevelPtr.resize(n + 1);
        LevelSet.resize(n);
        nlevels = build_levelSet_CSC_V2(n, DAG_ptr, DAG_set,
                                        LevelPtr.data(), LevelSet.data());
        LevelSet_time.measure_elapsed_time();
        return LevelSet_time;
    }

    /*
     * @brief assign the nodes' level to each node and store it in Node2Level
     * @param LevelPtr The pointer array in CSC format
     * @param LevelSet The index array in CSC format
     * @param nlevels Number of levels
     * @param num_nodes
     * @return Node2Level array that stores the data
     */
    timing_measurement computingNode2Level(const std::vector<int>& LevelPtr, const std::vector<int>& LevelSet, int nlevels,
                                           int num_nodes, std::vector<int>& Node2Level){
        timing_measurement Node2Level_time;
        Node2Level_time.start_timer();
        Node2Level.resize(num_nodes);
        for (int lvl = 0; lvl < nlevels; ++lvl) {
            for (int j = LevelPtr[lvl]; j < LevelPtr[lvl + 1]; ++j) {
                int node = LevelSet[j];
                Node2Level[node] = lvl;
            }
        }

        Node2Level_time.measure_elapsed_time();
        return Node2Level_time;
    }

    enum Algorithm { SERIAL, LEVELSET, GROUP_LEVELSET, BLOCK_LEVELSET, BLOCK_SERIAL};

    /*
    * Serial Code to have a base for correctness and runtime
    */
    //================================ SpTrSv Serial LL ================================
    class Sptrsv_LL_FLOPS_MEMORY : public FusionDemo {
    protected:
        double flops = 0;
        double double_element = 0;
        double int_element = 0;
        timing_measurement fused_code() override {
            timing_measurement t1;
            auto n = n_;
            auto Lp = L1_csr_->p;
            auto Li = L1_csr_->i;
            auto Lx = L1_csr_->x;
            auto x = x_in_;
            flops = 0;
            double_element = 0;
            int_element = 0;
            t1.start_timer();
            for (int i = 0; i < n; i++)
            {
                int_element+=2;
                double_element++;
                for (int j = Lp[i]; j < Lp[i + 1] - 1; j++)
                {
                    x[i] -= Lx[j] * x[Li[j]]; //S2
                    int_element++;
                    double_element+=2;
                    flops+=2;
                }
                x[i] /= Lx[Lp[i + 1] - 1]; //S1
                double_element++;
                flops++;
            }
            t1.measure_elapsed_time();
            copy_vector(0, n_, x_in_, x_);
            return t1;
        }

    public:
        Sptrsv_LL_FLOPS_MEMORY(CSR *L, CSC *L_csc,
                         double *correct_x, std::string name) : FusionDemo(L->n, L->nnz, name) {
            L1_csr_ = L;
            L1_csc_ = L_csc;
            correct_x_ = correct_x;
        };
        double getFLOPS() const{
            return flops;
        }

        double getBytes() const{
            return double_element * sizeof(double) + int_element * sizeof(int);
        }
        ~Sptrsv_LL_FLOPS_MEMORY() override = default;
    };
    //================================ SpTrSv Serial LL ================================
    class Sptrsv_LL_Serial : public FusionDemo {
    protected:
        timing_measurement fused_code() override {
            timing_measurement t1;
            t1.start_timer();
            sptrsv_csr(n_, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_);
            t1.measure_elapsed_time();
            copy_vector(0, n_, x_in_, x_);
            return t1;
        }

    public:
        Sptrsv_LL_Serial(CSR *L, CSC *L_csc,
                         double *correct_x, std::string name) : FusionDemo(L->n, L->nnz, name) {
            L1_csr_ = L;
            L1_csc_ = L_csc;
            correct_x_ = correct_x;
        };

        ~Sptrsv_LL_Serial() = default;
    };
    //================================ SpTrSv Blocked Vectorized Serial LL ================================
    class Sptrsv_LL_Blocked_Serial : public FusionDemo {
    protected:
        timing_measurement schedule_time;
        supernodaltools::SuperNodal* super_node_obj;
        void build_set() override {
            schedule_time.start_timer();
            super_node_obj = new supernodaltools::SuperNodal(supernodaltools::CSR_TYPE,
                                                             L1_csc_, L1_csr_, 0);
            schedule_time.measure_elapsed_time();
        }

        timing_measurement fused_code() override {
            timing_measurement t1;
            auto num_nodes = super_node_obj->getNumSuperNodes();
            auto supernode_ptr = super_node_obj->getSuperNodeGroups();
            auto Lp = L1_csr_->p;
            auto Li = L1_csr_->i;
            auto Lx = L1_csr_->x;
            auto n = L1_csr_->n;
            auto x = x_in_;
            std::vector<double> tempvec(n, 0);
            t1.start_timer();
            for (int super_ptr = 0; super_ptr < num_nodes; super_ptr++) {
                int start_row = supernode_ptr[super_ptr];
                int end_row = supernode_ptr[super_ptr + 1];
                //ncol is the number of columns in the off-diagonal block
                int num_off_diag_col = Lp[start_row + 1] - Lp[start_row] - 1;
                assert(num_off_diag_col >= 0);
                int nrows = end_row - start_row;
                assert(nrows > 0);
                //Solving the independent part
                //Copy x[Li[col]] into a continues buffer
                for (int col_ptr = Lp[supernode_ptr[super_ptr]], k = 0;
                     col_ptr < Lp[supernode_ptr[super_ptr]] + num_off_diag_col; col_ptr++, k++) {
                    tempvec[k] = x[Li[col_ptr]];
                }
                custom_blas::SpTrSv_MatVecCSR_BLAS(nrows, num_off_diag_col, &Lx[Lp[start_row]], tempvec.data(), &x[start_row]);
                custom_blas::SpTrSv_LSolveCSR_BLAS(nrows, num_off_diag_col, &Lx[Lp[start_row] + num_off_diag_col], &x[start_row]);
            }

            t1.measure_elapsed_time();
            copy_vector(0, n_, x_in_, x_);
            return t1;
        }

    public:
        Sptrsv_LL_Blocked_Serial(CSR *L, CSC *L_csc,
                                 double *correct_x, std::string name) : FusionDemo(L->n, L->nnz, name) {
            L1_csr_ = L;
            L1_csc_ = L_csc;
            correct_x_ = correct_x;
            super_node_obj = nullptr;
        };

        double getSchedulingTime(){return schedule_time.elapsed_time;}
        ~Sptrsv_LL_Blocked_Serial() {
            delete super_node_obj;
        };
    };
    //================================ SpTrSv Blocked Vectorized LevelSet LL ================================
    class SpTrsv_LL_Blocked_LevelSet : public Sptrsv_LL_Serial {
    protected:
        int nlevels;
        std::vector<int> level_ptr;
        std::vector<int> level_set;
        std::vector<int> DAG_ptr;
        std::vector<int> DAG_set;
        int num_supernode, num_dep_edges;
        timing_measurement scheduling_time;
        supernodaltools::SuperNodal* super_node_obj;

        void build_set() override {
            scheduling_time.start_timer();
            //Create the DAG
            assert(L1_csc_ != nullptr);
            assert(L1_csr_ != nullptr);
            super_node_obj = new supernodaltools::SuperNodal(supernodaltools::CSR_TYPE,
                                                             L1_csc_, L1_csr_,
                                                             L1_csc_, 0);
            super_node_obj->getDependencyDAGCSCformat(DAG_ptr, DAG_set, num_supernode, num_dep_edges);
            //Create the Levelset
            auto levelset_time =  computingLevelSet_CSC(num_supernode, DAG_ptr.data(), DAG_set.data(),
                                                        level_ptr, level_set, nlevels);
            scheduling_time.measure_elapsed_time();
        }

        timing_measurement fused_code() override {
            //sym_lib::rhs_init(L1_csc_->n, L1_csc_->p, L1_csc_->i, L1_csc_->x, x_); // x is b
            timing_measurement t1;
            auto supernode_ptr = super_node_obj->getSuperNodeGroups();
            auto Lp = L1_csr_->p;
            auto Li = L1_csr_->i;
            auto Lx = L1_csr_->x;
            auto n = L1_csr_->n;
            auto x = x_in_;
            t1.start_timer();
            #pragma omp parallel
            {
                std::vector<double> tempvec(n, 0);
                for (int l = 0; l < nlevels; l++)
                {
                    #pragma omp for schedule(auto)
                    for (int super_ptr = level_ptr[l]; super_ptr < level_ptr[l + 1]; ++super_ptr)
                    {
                        int super = level_set[super_ptr];
                        int start_row = supernode_ptr[super];
                        int end_row = supernode_ptr[super + 1];
                        //ncol is the number of columns in the off-diagonal block
                        int num_off_diag_col = Lp[start_row + 1] - Lp[start_row] - 1;
                        assert(num_off_diag_col >= 0);
                        int nrows = end_row - start_row;
                        assert(nrows > 0);
                        //Solving the independent part
                        //Copy x[Li[col]] into a continues buffer
                        for (int col_ptr = Lp[supernode_ptr[super]], k = 0;
                             col_ptr < Lp[supernode_ptr[super]] + num_off_diag_col; col_ptr++, k++) {
                            tempvec[k] = x[Li[col_ptr]];
                        }
                        custom_blas::SpTrSv_MatVecCSR_BLAS(nrows, num_off_diag_col, &Lx[Lp[start_row]], tempvec.data(), &x[start_row]);
                        custom_blas::SpTrSv_LSolveCSR_BLAS(nrows, num_off_diag_col, &Lx[Lp[start_row] + num_off_diag_col], &x[start_row]);
                    }
                }
            }
            t1.measure_elapsed_time();
            sym_lib::copy_vector(0, n_, x_in_, x_);
            return t1;
        }

    public:
        /*
         * @brief Class constructor
         * @param L Sparse Matrix with CSR format
         * @param L_csc Sparse Matrix with CSC format
         * @param correct_x The right answer for x
         * @param name The name of the algorithm
         * @param nt number of threads
         */
        SpTrsv_LL_Blocked_LevelSet(CSR *L, CSC *L_csc,
                                   double *correct_x, std::string name,
                                   int nt)
                : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
            L1_csr_ = L;
            L1_csc_ = L_csc;
            correct_x_ = correct_x;
        };
        double getSchedulingTime(){return scheduling_time.elapsed_time;}
        ~SpTrsv_LL_Blocked_LevelSet() = default;
    };
    //================================ SpTrSv LevelSet LL ================================
    class SpTrsv_LL_LevelSet : public Sptrsv_LL_Serial {
    protected:
        int nlevels, nthreads;
        std::vector<int> level_ptr;
        std::vector<int> level_set;
        std::vector<int> DAG_ptr;
        std::vector<int> DAG_set;
        timing_measurement scheduling_time;
        void build_set() override {
            //Create the DAG for automation stuff
            //auto DAG_time = computingDAG_CSR(L1_csr_, DAG_ptr, DAG_set);
            scheduling_time.start_timer();
            //Create the Levelset
            auto levelset_time =  computingLevelSet_CSC(n_, L1_csc_->p, L1_csc_->i,
                                                        level_ptr, level_set, nlevels);
            scheduling_time.measure_elapsed_time();
        }

        timing_measurement fused_code() override {
            timing_measurement t1;
            t1.start_timer();
            sptrsv_csr_levelset(L1_csr_->n, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_,
                                nlevels, level_ptr.data(), level_set.data());
            t1.measure_elapsed_time();
            sym_lib::copy_vector(0, n_, x_in_, x_);
            return t1;
        }

    public:
        /*
         * @brief Class constructor
         * @param L Sparse Matrix with CSR format
         * @param L_csc Sparse Matrix with CSC format
         * @param correct_x The right answer for x
         * @param name The name of the algorithm
         * @param nt number of threads
         */
        SpTrsv_LL_LevelSet(CSR *L, CSC *L_csc,
                           double *correct_x, std::string name,
                           int nt)
                           : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
            L1_csr_ = L;
            L1_csc_ = L_csc;
            correct_x_ = correct_x;
            this->nthreads = nt;
        };
        double getSchedulingTime(){return scheduling_time.elapsed_time;};
        void getWaveStatistic(int& nlevels, int& part_no){
            part_no = n_;
            nlevels = this->nlevels;
        }
        ~SpTrsv_LL_LevelSet() = default;
    };
    //================================ SpTrSv Parallel LevelSet LL ================================
    class SpTrsv_LL_Parallel_Levelset : public Sptrsv_LL_Serial {
    protected:
        int nlevels, nthreads;
        std::vector<int> DAG_ptr;
        std::vector<int> DAG_set;
        std::vector<int> level_ptr, level_set;
        timing_measurement scheduling_time;
        CSR* LBC_A;
        void build_set() override {
            //Convert lbc CSR version to SpMP CSR version
            SpMP::CSR* A = new SpMP::CSR();
            Convert_LBCCSR_to_SpMP(LBC_A, A);
            scheduling_time.start_timer();
            nlevels = GLC::levelsetCSRParallel_SpMP(A, nthreads, level_ptr, level_set);

            scheduling_time.measure_elapsed_time();
            delete A;
        }

        timing_measurement fused_code() override {
            timing_measurement t1;
            t1.start_timer();
            sptrsv_csr_levelset(L1_csr_->n, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_,
                                nlevels, level_ptr.data(), level_set.data());
            t1.measure_elapsed_time();
            sym_lib::copy_vector(0, n_, x_in_, x_);
            return t1;
        }

    public:
        /*
         * @brief Class constructor
         * @param L Sparse Matrix with CSR format
         * @param L_csc Sparse Matrix with CSC format
         * @param correct_x The right answer for x
         * @param name The name of the algorithm
         * @param nt number of threads
         */
        SpTrsv_LL_Parallel_Levelset(CSR* A, CSR *L, CSC *L_csc,
                       double *correct_x, std::string name,
                       int nt)
                : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
            this -> LBC_A = A;
            L1_csr_ = L;
            L1_csc_ = L_csc;
            correct_x_ = correct_x;
            this->nthreads = nt;
            omp_set_num_threads(nt);
        };

        double getSchedulingTime(){return scheduling_time.elapsed_time;};
        void getWaveStatistic(int& nlevels, int& part_no){
            part_no = n_;
            nlevels = this->nlevels;
        }

        ~SpTrsv_LL_Parallel_Levelset() {
        };
    };
    //================================ SpICh0 Levelset + Tree UL ================================
    class SpTrSv_LL_Tree_Levelset_No_unPACK : public Sptrsv_LL_Serial {
    protected:
        int nthreads;
        timing_measurement Scheduling_time;
        int ngroups;
        std::vector<int> DAG_ptr;
        std::vector<int> DAG_set;
        std::vector<int> group_set, group_ptr;
        std::vector<int> level_ptr, level_set;
        int nlevels;
        void build_set() override {
            Scheduling_time.start_timer();
            //Create the DAG
            GLC::partialSparsification(n_, nnz_, L1_csc_->p, L1_csc_->i, DAG_ptr, DAG_set);

            GLC::treeBasedGrouping(n_, DAG_ptr, DAG_set,
                                   ngroups, group_ptr, group_set, false);
            std::vector<int> group_DAG_ptr, group_DAG_set;
            GLC::buildGroupDAG(n_, ngroups, group_ptr.data(), group_set.data(),
                               DAG_ptr.data(), DAG_set.data(), group_DAG_ptr, group_DAG_set);

            computingLevelSet_CSC(ngroups, group_DAG_ptr.data(), group_DAG_set.data(),
                                  level_ptr, level_set, nlevels);

            for(int g = 0; g < ngroups; g++){
                std::sort(group_set.data() + group_ptr[g], group_set.data() + group_ptr[g + 1]);
            }
            Scheduling_time.measure_elapsed_time();
        }

        timing_measurement fused_code() override {
            timing_measurement t1;
            t1.start_timer();
            sptrsv_csr_group_levelset(L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_,
                                      nlevels, level_ptr.data(),
                                      level_set.data(), group_ptr.data(), group_set.data());
            t1.measure_elapsed_time();
            sym_lib::copy_vector(0, n_, x_in_, x_);
            return t1;
        }

    public:
        /*
         * @brief Class constructor
         * @param L Sparse Matrix with CSR format
         * @param L_csc Sparse Matrix with CSC format
         * @param correct_x The right answer for x
         * @param name The name of the algorithm
         * @param nt number of threads
         */
        SpTrSv_LL_Tree_Levelset_No_unPACK(CSR *L, CSC *L_csc,
                                          double *correct_x, std::string name,
                                          int nt)
                                          : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
            L1_csr_ = L;
            L1_csc_ = L_csc;
            correct_x_ = correct_x;
            nthreads = nt;
        };

        int getNlevels(){
            return nlevels;
        }
        void getWaveStatistic(int& nlevels, int& part_no){
            part_no = ngroups;
            nlevels = this->nlevels;
        }
        double getSchedulingTime() { return Scheduling_time.elapsed_time;}
        ~SpTrSv_LL_Tree_Levelset_No_unPACK() = default;
    };
    //================================ SpTrSv Block Vectorized LBC LL ================================
    class SpTrSv_LL_Blocked_LBC : public Sptrsv_LL_Serial {
    protected:
        int nthreads;
        timing_measurement Scheduling_time1, Scheduling_time2;
        int num_supernodes, num_dep_edges;
        std::vector<int> supernode_ptr;
        int num_w_part;
        int lp_, cp_, ic_;
        std::vector<int> DAG_ptr;
        std::vector<int> DAG_set;
        std::vector<int> DAG_inv_ptr;
        std::vector<int> DAG_inv_set;
        std::vector<double> cost;
        supernodaltools::SuperNodal* supernode_obj;
        int final_level_no;
        int *final_level_ptr, *final_part_ptr, *final_node_ptr;
        int part_no;
        bool first_time;
        void build_set() override {
            Scheduling_time2.start_timer();
            assert(L1_csc_ != nullptr);
            assert(L1_csr_ != nullptr);
            if(!first_time){
                delete final_level_ptr;
                delete final_node_ptr;
                delete final_part_ptr;
                first_time = false;
            }
            get_coarse_levelSet_DAG_CSC_tree(num_supernodes, DAG_inv_ptr.data(), DAG_inv_set.data(), -1,
                                             final_level_no,
                                             final_level_ptr,part_no,
                                             final_part_ptr,final_node_ptr,
                                             lp_, ic_, cp_, cost.data());
            Scheduling_time2.measure_elapsed_time();
        }

        timing_measurement fused_code() override {
            //sym_lib::rhs_init(L1_csc_->n, L1_csc_->p, L1_csc_->i, L1_csc_->x, x_); // x is b
            timing_measurement t1;
            auto Lp = L1_csr_->p;
            auto Li = L1_csr_->i;
            auto Lx = L1_csr_->x;
            auto x = x_in_;

            t1.start_timer();
            #pragma omp parallel
            {
                std::vector<double> tempvec(n_, 0);
                for (int lvl = 0; lvl < final_level_no; ++lvl) {
                    #pragma omp for
                    for (int w_ptr = final_level_ptr[lvl]; w_ptr < final_level_ptr[lvl + 1]; ++w_ptr) {
                        for(int super_ptr = final_part_ptr[w_ptr]; super_ptr < final_part_ptr[w_ptr + 1]; super_ptr++){
                            int super = final_node_ptr[super_ptr];
                            int start_row = supernode_ptr[super];
                            int end_row = supernode_ptr[super + 1];
                            //ncol is the number of columns in the off-diagonal block
                            int num_off_diag_col = Lp[start_row + 1] - Lp[start_row] - 1;
                            assert(num_off_diag_col >= 0);
                            int nrows = end_row - start_row;
                            assert(nrows > 0);
                            //Solving the independent part
                            //Copy x[Li[col]] into a continues buffer
                            for (int col_ptr = Lp[supernode_ptr[super]], k = 0;
                                 col_ptr < Lp[supernode_ptr[super]] + num_off_diag_col; col_ptr++, k++) {
                                tempvec[k] = x[Li[col_ptr]];
                            }
                            custom_blas::SpTrSv_MatVecCSR_BLAS(nrows, num_off_diag_col, &Lx[Lp[start_row]], tempvec.data(), &x[start_row]);
                            custom_blas::SpTrSv_LSolveCSR_BLAS(nrows, num_off_diag_col, &Lx[Lp[start_row] + num_off_diag_col], &x[start_row]);
                        }
                    }
                }
            }
            t1.measure_elapsed_time();
            sym_lib::copy_vector(0, n_, x_in_, x_);
            return t1;
        }

    public:
        /*
         * @brief Class constructor
         * @param L Sparse Matrix with CSR format
         * @param L_csc Sparse Matrix with CSC format
         * @param correct_x The right answer for x
         * @param name The name of the algorithm
         * @param nt number of threads
         */
        SpTrSv_LL_Blocked_LBC(CSR *L, CSC *L_csc,
                              double *correct_x, std::string name,
                              int lp)
                : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
            L1_csr_ = L;
            L1_csc_ = L_csc;
            correct_x_ = correct_x;
            lp_=lp;
            nthreads = lp;


            //Calculating The Dependency DAG and the levelset
            timing_measurement block_LBC_time;
            Scheduling_time1.start_timer();
            supernodaltools::SuperNodal super_node_obj(supernodaltools::CSR_TYPE,
                                                                  L1_csc_, L1_csr_, L1_csc_, 0);

            super_node_obj.getDependencyDAGCSCformat(DAG_ptr, DAG_set,num_supernodes, num_dep_edges);

            CSC DAG(num_supernodes,num_supernodes, num_dep_edges);
            std::copy(DAG_ptr.begin(), DAG_ptr.end(), DAG.p);
            std::copy(DAG_set.begin(), DAG_set.end(), DAG.i);
            auto DAG_inv = csc_to_csr(&DAG);

            DAG_inv_ptr.resize(num_supernodes + 1);
            DAG_inv_set.resize(num_dep_edges);
            std::copy(DAG_inv->p, DAG_inv->p + num_supernodes + 1, DAG_inv_ptr.begin());
            std::copy(DAG_inv->i, DAG_inv->i + num_dep_edges, DAG_inv_set.begin());

            cost.resize(num_supernodes, 0);
            super_node_obj.getSuperNodeGroups(supernode_ptr);
            for (int super = 0; super < num_supernodes; ++super) {
                for(int i = supernode_ptr[super]; i < supernode_ptr[super + 1]; i++){
                    cost[super] += L1_csr_->p[i + 1] - L1_csr_->p[i];
                }
            }
            first_time = true;
            Scheduling_time1.start_timer();
        };

        double getSchedulingTime() { return Scheduling_time1.elapsed_time + Scheduling_time2.elapsed_time; }

        void setP2_P3(int p2, int p3) {
            this->cp_ = p2;
            this->ic_ = p3;
        }
        ~SpTrSv_LL_Blocked_LBC() = default;
    };
    //================================ SpTrSv Blocked LC LL ================================
    class SpTrSv_LL_Blocked_LC : public Sptrsv_LL_Serial {
    protected:
        int nthreads;
        int max_level;
        timing_measurement Scheduling_time1, Scheduling_time2;
        std::vector<int> blocked_level_set, blocked_level_ptr, blocked_node_to_level;
        int blocked_levelNo, num_supernodes, num_dep_edges;
        std::vector<int> supernode_ptr;
        std::vector<int> LBC_level_ptr, LBC_w_ptr, LBC_node_ptr;
        int LBC_level_no;
        int num_w_part;
        std::vector<int> WM;
        std::vector<int> DAG_ptr;
        std::vector<int> DAG_set;
        std::vector<double> cost;


        //First partition parameters
        double max_cc_nodes;
        double min_cc_nodes;
        double whole_nodes;
        double max_out_degree;
        double max_in_degree;

        void build_set() override {
            Scheduling_time2.start_timer();
            assert(L1_csc_ != nullptr);
            assert(L1_csr_ != nullptr);

            LBC_level_ptr.clear();
            LBC_w_ptr.clear();
            LBC_node_ptr.clear();

            if(WM.empty()){
                std::cerr << "The WM is not defined" << std::endl;
            }

            skeletonCoarsening(
                    num_supernodes, DAG_ptr.data(), DAG_set.data(),
                    LBC_level_no, LBC_level_ptr, num_w_part,
                    LBC_w_ptr, LBC_node_ptr,
                    WM, nthreads,
                    blocked_levelNo, blocked_level_set.data(), blocked_level_ptr.data(),
                    blocked_node_to_level.data());

            Scheduling_time2.measure_elapsed_time();
        }

        timing_measurement fused_code() override {
            timing_measurement t1;
            auto Lp = L1_csr_->p;
            auto Li = L1_csr_->i;
            auto Lx = L1_csr_->x;
            auto n = L1_csr_->n;
            auto x = x_in_;
            t1.start_timer();
            #pragma omp parallel
            {
                std::vector<double> tempvec(n);
                for (int lvl = 0; lvl < LBC_level_no; ++lvl) {
                    #pragma omp for schedule(static)
                    for (int w_ptr = LBC_level_ptr[lvl]; w_ptr < LBC_level_ptr[lvl + 1]; ++w_ptr) {
                        for(int super_ptr = LBC_w_ptr[w_ptr]; super_ptr < LBC_w_ptr[w_ptr + 1]; super_ptr++){
                            int super = LBC_node_ptr[super_ptr];
                            int start_row = supernode_ptr[super];
                            int end_row = supernode_ptr[super + 1];
                            //ncol is the number of columns in the off-diagonal block
                            int num_off_diag_col = Lp[start_row + 1] - Lp[start_row] - 1;
                            assert(num_off_diag_col >= 0);
                            int nrows = end_row - start_row;
                            assert(nrows > 0);
                            //Solving the independent part
                            //Copy x[Li[col]] into a continues buffer
                            for (int col_ptr = Lp[supernode_ptr[super]], k = 0;
                            col_ptr < Lp[supernode_ptr[super]] + num_off_diag_col; col_ptr++, k++) {
                                tempvec[k] = x[Li[col_ptr]];
                            }
                            custom_blas::SpTrSv_MatVecCSR_BLAS(nrows, num_off_diag_col, &Lx[Lp[start_row]], tempvec.data(), &x[start_row]);
                            custom_blas::SpTrSv_LSolveCSR_BLAS(nrows, num_off_diag_col, &Lx[Lp[start_row] + num_off_diag_col], &x[start_row]);
                        }
                    }
                }
            }
            t1.measure_elapsed_time();
            sym_lib::copy_vector(0, n_, x_in_, x_);
            return t1;
        }

    public:
        /*
         * @brief Class constructor
         * @param L Sparse Matrix with CSR format
         * @param L_csc Sparse Matrix with CSC format
         * @param correct_x The right answer for x
         * @param name The name of the algorithm
         * @param nt number of threads
         */
        SpTrSv_LL_Blocked_LC(CSR *L, CSC *L_csc, double *correct_x, std::string name, int lp)
        : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
            L1_csr_ = L;
            L1_csc_ = L_csc;
            correct_x_ = correct_x;
            nthreads = lp;

            //Calculating The Dependency DAG and the levelset
            timing_measurement block_LBC_time;
            Scheduling_time1.start_timer();
            supernodaltools::SuperNodal supernode_obj(supernodaltools::CSR_TYPE,
                                                      L1_csc_, L1_csr_, L1_csc_, 0);

            supernode_obj.getDependencyDAGCSCformat(DAG_ptr, DAG_set,num_supernodes, num_dep_edges);

            computingLevelSet_CSC(num_supernodes, DAG_ptr.data(), DAG_set.data(),
                                  blocked_level_ptr, blocked_level_set, blocked_levelNo);

            computingNode2Level(blocked_level_ptr, blocked_level_set,
                                blocked_levelNo, num_supernodes, blocked_node_to_level);

            cost.resize(num_supernodes, 0);
            supernode_obj.getSuperNodeGroups(supernode_ptr);
            for (int super = 0; super < num_supernodes; ++super) {
                for(int i = supernode_ptr[super]; i < supernode_ptr[super + 1]; i++){
                    cost[super] += L1_csc_->p[i + 1] - L1_csc_->p[i];
                }
            }
            Scheduling_time1.start_timer();
        };

        double getSchedulingTime() { return Scheduling_time1.elapsed_time + Scheduling_time2.elapsed_time; }

        double getMaxCC(){
            return max_cc_nodes;
        }

        double getminCC(){
            return min_cc_nodes;
        }

        double getMaxToWholeRatio(){
            return max_cc_nodes / whole_nodes;
        }

        double getMaxOutDegree(){
            return max_out_degree;
        }

        double getMaxInDegree(){
            return max_in_degree;
        }

        double getNumLevels(){
            return blocked_levelNo;
        }

        void setMaxLevel(int max_level){
            this->max_level = max_level;
        }

        void setWM(int wm){
            WM.clear();
            for(int j = 0; j < blocked_levelNo + 1; j+=wm){
                WM.push_back(j);
            }
            if(WM.back() != blocked_levelNo + 1){
                WM.push_back(blocked_levelNo + 1);
            }
        }


        void computeCharacteristics(double& first, double& second, double& third, double& avg){

            std::vector<double> max_array(LBC_level_no, 0);
            std::vector<double> whole_array(LBC_level_no, 0);
            for(int lvl = 0; lvl < LBC_level_no; lvl++){
                double total_nodes = 0;
                double max_size = 0;
                for (int w_ptr = LBC_level_ptr[lvl]; w_ptr < LBC_level_ptr[lvl + 1]; ++w_ptr) {
                    int cc_size = LBC_w_ptr[w_ptr + 1] - LBC_w_ptr[w_ptr];
                    if(cc_size > max_size){
                        max_size = cc_size;
                    }
                    total_nodes += cc_size;
                }
                if(lvl == 0){
                    first = max_size / total_nodes;
                }
                if(lvl == 1){
                    second = max_size / total_nodes;
                }
                if(lvl == 2){
                    third = max_size / total_nodes;
                }
                max_array[lvl] = max_size;
                whole_array[lvl] = total_nodes;
            }

            double max_cc_sizes = 0;
            for(auto& iter: max_array){
                max_cc_sizes += iter;
            }
            avg = max_cc_sizes / n_;

            #ifndef NDEBUG
            double sum = 0;
            for(auto &iter: whole_array){
                sum += iter;
            }
            assert(sum == n_);
            #endif
        }

        ~SpTrSv_LL_Blocked_LC()=default;
    };
    //================================ SpTrSv LC LL ================================
    class SpTrSv_LL_LC : public Sptrsv_LL_Serial {
    protected:
        int nthreads;
        timing_measurement Scheduling_time1;
        timing_measurement Scheduling_time2;
        std::vector<int> DAG_ptr;
        std::vector<int> DAG_set;
        std::vector<int> level_set, level_ptr, node_to_level;
        int levelNo;
        std::vector<int> final_level_ptr, final_part_ptr, final_node_ptr;
        int final_level_no;
        int part_no;
        std::vector<int> WM;
        void build_set() override {
            Scheduling_time2.start_timer();
            std::vector<double> cost(n_, 0);
            for (int i = 0; i < n_; ++i) {
                cost[i] = L1_csr_->p[i + 1] - L1_csr_->p[i];
            }
            final_level_ptr.clear();
            final_part_ptr.clear();
            final_node_ptr.clear();
            if(WM.empty()){
                std::cerr << "The WM is not defined" << std::endl;
            }
            skeletonCoarsening(
                    n_, L1_csc_->p, L1_csc_->i,
                    final_level_no, final_level_ptr, part_no,
                    final_part_ptr, final_node_ptr,
                    WM, nthreads,
                    levelNo, level_set.data(), level_ptr.data(),
                    node_to_level.data());

            Scheduling_time2.measure_elapsed_time();
        }

        timing_measurement fused_code() override {
            timing_measurement t1;
            t1.start_timer();
            sptrsv_csr_lbc(n_, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_,
                           final_level_no, final_level_ptr.data(),
                           final_part_ptr.data(), final_node_ptr.data());
            t1.measure_elapsed_time();
            sym_lib::copy_vector(0, n_, x_in_, x_);
            return t1;
        }

    public:
        /*
         * @brief Class constructor
         * @param L Sparse Matrix with CSR format
         * @param L_csc Sparse Matrix with CSC format
         * @param correct_x The right answer for x
         * @param name The name of the algorithm
         * @param nt number of threads
         */
        SpTrSv_LL_LC(CSR *L, CSC *L_csc, double *correct_x, std::string name, int lp)
        : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
            L1_csr_ = L;
            L1_csc_ = L_csc;
            correct_x_ = correct_x;
            nthreads = lp;

            Scheduling_time1.start_timer();
//            //Create the DAG
//            computingDAG_CSR(L1_csr_, DAG_ptr, DAG_set);
            //Create the Levelset
            computingLevelSet_CSC(n_, L1_csc_->p, L1_csc_->i, level_ptr, level_set, levelNo);

            computingNode2Level(level_ptr, level_set, levelNo, n_, node_to_level);
            Scheduling_time1.measure_elapsed_time();
        };

        double getSchedulingTime() { return Scheduling_time1.elapsed_time + Scheduling_time2.elapsed_time; }

        void setWM(int wm){
            WM.clear();
            for(int j = 0; j < levelNo + 1; j+=wm){
                WM.push_back(j);
            }
            if(WM.back() != levelNo + 1){
                WM.push_back(levelNo + 1);
            }
        }

        void computeCharacteristics(double& first, double& second, double& third, double& avg){

            std::vector<double> max_array(final_level_no, 0);
            std::vector<double> whole_array(final_level_no, 0);
            for(int lvl = 0; lvl < final_level_no; lvl++){
                double total_nodes = 0;
                double max_size = 0;
                for (int w_ptr = final_level_ptr[lvl]; w_ptr < final_level_ptr[lvl + 1]; ++w_ptr) {
                    int cc_size = final_part_ptr[w_ptr + 1] - final_part_ptr[w_ptr];
                    if(cc_size > max_size){
                        max_size = cc_size;
                    }
                    total_nodes += cc_size;
                }
                if(lvl == 0){
                    first = max_size / total_nodes;
                }
                if(lvl == 1){
                    second = max_size / total_nodes;
                }
                if(lvl == 2){
                    third = max_size / total_nodes;
                }
                max_array[lvl] = max_size;
                whole_array[lvl] = total_nodes;
            }

            double max_cc_sizes = 0;
            for(auto& iter: max_array){
                max_cc_sizes += iter;
            }
            avg = max_cc_sizes / n_;

            #ifndef NDEBUG
            double sum = 0;
            for(auto &iter: whole_array){
                sum += iter;
            }
            assert(sum == n_);
            #endif
        }
        ~SpTrSv_LL_LC() = default;
    };
    //================================ SpTrSv LC LL ================================
    class SpTrSv_LL_DLC : public Sptrsv_LL_Serial {
    protected:
        int nthreads;
        timing_measurement Scheduling_time1;
        timing_measurement Scheduling_time2;
        std::vector<int> DAG_ptr;
        std::vector<int> DAG_set;
        std::vector<int> level_set, level_ptr, node_to_level;
        int levelNo;
        std::vector<int> final_level_ptr, final_part_ptr, final_node_ptr;
        int final_level_no;
        int part_no;
        std::vector<int> WM;
        void build_set() override {
            Scheduling_time2.start_timer();
            std::vector<double> cost(n_, 0);
            for (int i = 0; i < n_; ++i) {
                cost[i] = L1_csr_->p[i + 1] - L1_csr_->p[i];
            }
            final_level_ptr.clear();
            final_part_ptr.clear();
            final_node_ptr.clear();
            computeWM();
            if(WM.empty()){
                std::cerr << "The WM is not defined" << std::endl;
            }
            skeletonCoarsening(
                    n_, DAG_ptr.data(), DAG_set.data(),
                    final_level_no, final_level_ptr, part_no,
                    final_part_ptr, final_node_ptr,
                    WM, nthreads,
                    levelNo, level_set.data(), level_ptr.data(),
                    node_to_level.data());

            Scheduling_time2.measure_elapsed_time();
        }

        timing_measurement fused_code() override {
            timing_measurement t1;
            t1.start_timer();
            sptrsv_csr_lbc(n_, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_,
                           final_level_no, final_level_ptr.data(),
                           final_part_ptr.data(), final_node_ptr.data());
            t1.measure_elapsed_time();
            sym_lib::copy_vector(0, n_, x_in_, x_);
            return t1;
        }

    public:
        /*
         * @brief Class constructor
         * @param L Sparse Matrix with CSR format
         * @param L_csc Sparse Matrix with CSC format
         * @param correct_x The right answer for x
         * @param name The name of the algorithm
         * @param nt number of threads
         */
        SpTrSv_LL_DLC(CSR *L, CSC *L_csc, double *correct_x, std::string name, int lp)
        : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
            L1_csr_ = L;
            L1_csc_ = L_csc;
            correct_x_ = correct_x;
            nthreads = lp;

            Scheduling_time1.start_timer();
            //Create the DAG
//            std::vector<int> unprune_DAG_ptr, unprune_DAG_set;
//            computingDAG_CSR(L1_csr_, unprune_DAG_ptr, unprune_DAG_set);
            GLC::partialSparsification(n_, nnz_, L1_csc_->p, L1_csc_->i, DAG_ptr, DAG_set);
            //Create the Levelset
            computingLevelSet_CSC(n_, L1_csc_->p, L1_csc_->i, level_ptr, level_set, levelNo);

            computingNode2Level(level_ptr, level_set, levelNo, n_, node_to_level);
            Scheduling_time1.measure_elapsed_time();
        };

        double getSchedulingTime() { return Scheduling_time1.elapsed_time + Scheduling_time2.elapsed_time; }

        int nodeExist(std::vector<int>& node_set, int node_id){
            int cnt = 0;
            for(auto &iter: node_set){
                if(iter == node_id){
                    return cnt;
                }
                cnt++;
            }
            return cnt;
        }

        double computeGroupCriticalPath(int group_id, std::vector<int>& group_ptr, std::vector<int>& group_set){
            int group_size = group_ptr[group_id + 1] - group_ptr[group_id];

            std::vector<int> node_set(group_set.data() + group_ptr[group_id], group_set.data() + group_ptr[group_id + 1]);
            std::sort(node_set.begin(), node_set.end());

            std::vector<int> in_degree(group_size, 0);
            std::vector<bool> visited(group_size, false);

            double lvl = 0;
            int outgoing_more_than_1 = 0;
            for(auto& iter : node_set){
                if( (DAG_ptr[iter + 1] - DAG_ptr[iter] - 1) > 1){
                    outgoing_more_than_1++;
                }
                for(int child_ptr = DAG_ptr[iter] + 1; child_ptr < DAG_ptr[iter + 1]; child_ptr++){
                    int child = DAG_set[child_ptr];
                    int index = nodeExist(node_set, child);
                    if(index != node_set.size()){
                        in_degree[index]++;
                    }
                }
            }
            assert(outgoing_more_than_1 < 2);

            bool all_node_visited = false;
            while(!all_node_visited){
                std::vector<int> visited_node_set;
                for(auto &iter: node_set){
                    int index = nodeExist(node_set, iter);
                    assert(index < group_size);
                    if(in_degree[index] == 0 && !visited[index]){
                        visited[index] = true;
                        visited_node_set.push_back(iter);
                    }
                }
                lvl++;
                if(lvl > group_size){
                    return -1;//Group has a cycle
                }

                for(auto& iter : visited_node_set){
                    for(int child_ptr = DAG_ptr[iter] + 1; child_ptr < DAG_ptr[iter + 1]; child_ptr++){
                        int child = DAG_set[child_ptr];
                        int index = nodeExist(node_set, child);
                        if(index != node_set.size()){
                            in_degree[index]--;
                        }
                    }
                }

                all_node_visited = true;
                for(int i = 0; i < visited.size(); i++){
                    if(!visited[i]){
                        all_node_visited = false;
                        break;
                    }
                }
            }
            assert(lvl > 0);
            return lvl; //return number of level = critical path
        }

        void computeWM(){

            std::vector<int> group_ptr, group_set;
            int ngroups;
            GLC::treeBasedGrouping(n_, DAG_ptr, DAG_set,
                                   ngroups, group_ptr, group_set, false);
            std::vector<int> group_DAG_ptr, group_DAG_set;
            GLC::buildGroupDAG(n_, ngroups, group_ptr.data(), group_set.data(), DAG_ptr.data(), DAG_set.data(),
                               group_DAG_ptr, group_DAG_set);

            std::vector<int> group_level_ptr, group_levelset;
            int nlevels;
            computingLevelSet_CSC(ngroups, group_DAG_ptr.data(), group_DAG_set.data(), group_level_ptr, group_levelset, nlevels);

            int wm = 0;
            WM.push_back(wm);
            for(int i = 0; i < nlevels; i++){
                double group_avg = 0;
                int cnt = 0;
                for(int j = group_level_ptr[i]; j < group_level_ptr[i + 1]; j++){
                    int group = group_levelset[j];
                    group_avg +=  computeGroupCriticalPath(group, group_ptr, group_set);
                    cnt++;
                }
                group_avg = group_avg / cnt;
                if(group_avg - floor(group_avg) > 0.5){
                    wm += group_avg;
                    wm++;
                    if(wm < levelNo + 1){
                        WM.push_back(wm);
                    }
                } else {
                    wm += group_avg;
                    if(wm < levelNo + 1){
                        WM.push_back(wm);
                    }
                }
            }
            if(WM.back() != levelNo + 1){
                WM.push_back(levelNo + 1);
            }

        }


        ~SpTrSv_LL_DLC() = default;
    };
    //================================ SpICh0 GLC UL ================================
    class SpTrSv_LL_GLC : public Sptrsv_LL_Serial {
    protected:
        int nthreads;
        timing_measurement Scheduling_time;
        int ngroups;
        std::vector<int> final_level_ptr, final_part_ptr, final_node_ptr;
        int final_level_no;
        std::vector<int> level_ptr, level_set;
        int nlevels;
        double avg_par, max_diff, var;
        int part_no;
        bool bin_pack;
        std::vector<double> cost;
        void build_set() override {
            Scheduling_time.start_timer();
            //Computing node Cost
            cost.resize(n_, 0);
            auto CSC_Lp = L1_csc_->p;
            auto CSC_Li = L1_csc_->i;
            auto CSR_Lp = L1_csr_->p;
            auto CSR_Li = L1_csr_->i;
            GLC::costComputation(n_, CSC_Lp, CSC_Li, CSR_Lp, CSR_Li,
                                 GLC::SpTrSv_LL, NULL, NULL, false, cost);

            std::vector<int> DAG_ptr(L1_csc_->p, L1_csc_->p + n_ + 1);
            std::vector<int> DAG_set(L1_csc_->i, L1_csc_->i + nnz_);
            GLC::GLC(n_, nnz_, DAG_ptr, DAG_set, cost, nthreads,
                     final_level_no, final_level_ptr, final_part_ptr, final_node_ptr,
                     false, false, bin_pack);
            Scheduling_time.measure_elapsed_time();
        }

        timing_measurement fused_code() override {
            timing_measurement t1;
            t1.start_timer();
            sptrsv_csr_lbc(n_, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_,
                           final_level_no, final_level_ptr.data(),
                           final_part_ptr.data(), final_node_ptr.data());
            t1.measure_elapsed_time();
            sym_lib::copy_vector(0, n_, x_in_, x_);
            return t1;
        }

    public:
        /*
         * @brief Class constructor
         * @param L Sparse Matrix with CSR format
         * @param L_csc Sparse Matrix with CSC format
         * @param correct_x The right answer for x
         * @param name The name of the algorithm
         * @param nt number of threads
         */
        SpTrSv_LL_GLC(CSR *L, CSC *L_csc,
                           double *correct_x, std::string name,
                           int nt, bool bin_pack)
                           : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
            L1_csr_ = L;
            L1_csc_ = L_csc;
            correct_x_ = correct_x;
            nthreads = nt;
            this->bin_pack = bin_pack;
        };
        //I'm not in the mood to delete codes
        void computeCharacteristics(double& first, double& second, double& third, double& avg){

            std::vector<double> max_array(final_level_no, 0);
            std::vector<double> whole_array(final_level_no, 0);
            for(int lvl = 0; lvl < final_level_no; lvl++){
                double total_nodes = 0;
                double max_size = 0;
                for (int w_ptr = final_level_ptr[lvl]; w_ptr < final_level_ptr[lvl + 1]; ++w_ptr) {
                    int cc_size = final_part_ptr[w_ptr + 1] - final_part_ptr[w_ptr];
                    if(cc_size > max_size){
                        max_size = cc_size;
                    }
                    total_nodes += cc_size;
                }
                if(lvl == 0){
                    first = max_size / total_nodes;
                }
                if(lvl == 1){
                    second = max_size / total_nodes;
                }
                if(lvl == 2){
                    third = max_size / total_nodes;
                }
                max_array[lvl] = max_size;
                whole_array[lvl] = total_nodes;
            }

            double max_cc_sizes = 0;
            for(auto& iter: max_array){
                max_cc_sizes += iter;
            }
            avg = max_cc_sizes / n_;

            #ifndef NDEBUG
            double sum = 0;
            for(auto &iter: whole_array){
                sum += iter;
            }
            assert(sum == n_);
            #endif
        }

        void computeStatistic(){
            double part_no = 0;
            for(int i = 0; i < final_part_ptr.size() - 1; i++){
                if(final_part_ptr[i + 1] != final_part_ptr[i]){
                    part_no++;
                }
            }
            this->part_no = part_no;
            avg_par = part_no/ final_level_no;
            std::vector<double> max_part(final_level_no, 0);
            std::vector<double> min_part(final_level_no, 0);
            std::vector<double> level_weight(final_level_no, 0);
            for(int lvl = 0; lvl < final_level_no; lvl++){
                double total_nodes = 0;
                double max_size = 0;
                double min_size = n_;
                for (int w_ptr = final_level_ptr[lvl]; w_ptr < final_level_ptr[lvl + 1]; ++w_ptr) {
                    double cc_size = final_part_ptr[w_ptr + 1] - final_part_ptr[w_ptr];
                    double cc_cost = 0;
                    for(int i = final_part_ptr[w_ptr]; i < final_part_ptr[w_ptr + 1]; i++){
                        int node = final_node_ptr[i];
                        cc_cost += cost[node];
                    }
                    cc_size = cc_cost;
                    if(cc_size > max_size){
                        max_size = cc_size;
                    }
                    if(cc_size < min_size){
                        min_size = cc_size;
                    }
                    total_nodes += cc_size;
                }
                max_part[lvl] = max_size;
                min_part[lvl] = min_size;
                level_weight[lvl] = total_nodes;
            }

            max_diff = 0;
            double total_weight = 0;
            for(int l = 0; l < final_level_no; l++){
                max_diff += level_weight[l] * (max_part[l] - min_part[l]);
                total_weight += level_weight[l];
            }
            max_diff = max_diff / total_weight;
            var = 0;
        }
        void getWaveStatistic(int& nlevels, int& part_no){
            part_no = this->part_no;
            nlevels = this->final_level_no;
        }
        double getSchedulingTime() { return Scheduling_time.elapsed_time;}
        void getStatistic(double& avg_par, double& max_diff, double& var){
            avg_par = this->avg_par;
            max_diff = this->max_diff;
            var = this->var;
        }
        ~SpTrSv_LL_GLC() = default;
    };
    //================================ SpICh0 GLC UL ================================
    class SpTrSv_LL_Tree_GLC : public Sptrsv_LL_Serial {
    protected:
        int nthreads;
        timing_measurement Scheduling_time;
        int ngroups;
        std::vector<int> final_level_ptr, final_part_ptr, final_node_ptr;
        int final_level_no;
        std::vector<int> DAG_ptr;
        std::vector<int> DAG_set;
        std::vector<int> group_set, group_ptr;
        std::vector<int> level_ptr, level_set;
        int nlevels;
        bool isLfactor;
        double avg_par, max_diff, var;
        int part_no;
        bool bin_pack;
        std::vector<double> cost;
        void build_set() override {
            Scheduling_time.start_timer();
            //Create the DAG
            if(isLfactor){
                std::cout << "Building Tree" << std::endl;
                GLC::buildETree(L1_csr_, nnz_, DAG_ptr, DAG_set);
            } else {
                #ifndef NDEBUG
                    std::cout << "Compute DAG" << std::endl;
                #endif
//                computingDAG_CSR(L1_csr_, L1_csc_->p, L1_csc_->i);
                GLC::partialSparsification(n_, nnz_, L1_csc_->p, L1_csc_->i, DAG_ptr, DAG_set);
            }
            GLC::treeBasedGrouping(n_, DAG_ptr, DAG_set,
                                   ngroups, group_ptr, group_set, isLfactor);
            std::vector<int> group_DAG_ptr, group_DAG_set;
            GLC::buildGroupDAG(n_, ngroups, group_ptr.data(), group_set.data(), DAG_ptr.data(), DAG_set.data(), group_DAG_ptr, group_DAG_set);

            //Computing node Cost
            cost.resize(ngroups, 0);
            auto CSC_Lp = L1_csc_->p;
            auto CSC_Li = L1_csc_->i;
            auto CSR_Lp = L1_csr_->p;
            auto CSR_Li = L1_csr_->i;
            GLC::costComputation(ngroups, CSC_Lp, CSC_Li, CSR_Lp, CSR_Li,
                                 GLC::SpTrSv_LL, group_ptr.data(), group_set.data(), true, cost);

            GLC::GLC(ngroups, group_DAG_ptr[ngroups], group_DAG_ptr, group_DAG_set, cost, nthreads,
                        final_level_no, final_level_ptr, final_part_ptr, final_node_ptr,false, false, bin_pack);

            if(ANALYSE_STATISTIC){
                std::cout << "Static is COMPUTING" << std::endl;
                computeStatistic();
            }

            GLC::ungroupingScheduleAndApplyOrdering(n_, final_level_no, final_level_ptr, final_part_ptr, final_node_ptr,
                                                    group_ptr, group_set);
            Scheduling_time.measure_elapsed_time();
        }

        timing_measurement fused_code() override {
            timing_measurement t1;
            t1.start_timer();
            sptrsv_csr_lbc(n_, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_,
                           final_level_no, final_level_ptr.data(),
                           final_part_ptr.data(), final_node_ptr.data());
            t1.measure_elapsed_time();
            sym_lib::copy_vector(0, n_, x_in_, x_);
            return t1;
        }

    public:
        /*
         * @brief Class constructor
         * @param L Sparse Matrix with CSR format
         * @param L_csc Sparse Matrix with CSC format
         * @param correct_x The right answer for x
         * @param name The name of the algorithm
         * @param nt number of threads
         */
        SpTrSv_LL_Tree_GLC(CSR *L, CSC *L_csc,
                           double *correct_x, std::string name,
                           int nt, bool isLfactor, bool bin_pack)
                           : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
            L1_csr_ = L;
            L1_csc_ = L_csc;
            correct_x_ = correct_x;
            nthreads = nt;
            this->isLfactor = isLfactor;
            this->bin_pack = bin_pack;
        };
        //I'm not in the mood to delete codes
        void computeCharacteristics(double& first, double& second, double& third, double& avg){

            std::vector<double> max_array(final_level_no, 0);
            std::vector<double> whole_array(final_level_no, 0);
            for(int lvl = 0; lvl < final_level_no; lvl++){
                double total_nodes = 0;
                double max_size = 0;
                for (int w_ptr = final_level_ptr[lvl]; w_ptr < final_level_ptr[lvl + 1]; ++w_ptr) {
                    int cc_size = final_part_ptr[w_ptr + 1] - final_part_ptr[w_ptr];
                    if(cc_size > max_size){
                        max_size = cc_size;
                    }
                    total_nodes += cc_size;
                }
                if(lvl == 0){
                    first = max_size / total_nodes;
                }
                if(lvl == 1){
                    second = max_size / total_nodes;
                }
                if(lvl == 2){
                    third = max_size / total_nodes;
                }
                max_array[lvl] = max_size;
                whole_array[lvl] = total_nodes;
            }

            double max_cc_sizes = 0;
            for(auto& iter: max_array){
                max_cc_sizes += iter;
            }
            avg = max_cc_sizes / n_;

            #ifndef NDEBUG
            double sum = 0;
            for(auto &iter: whole_array){
                sum += iter;
            }
            assert(sum == n_);
            #endif
        }

        void computeStatistic(){
            double part_no = 0;
            for(int i = 0; i < final_part_ptr.size() - 1; i++){
                if(final_part_ptr[i + 1] != final_part_ptr[i]){
                    part_no++;
                }
            }
            this->part_no = part_no;
            avg_par = part_no/ final_level_no;
            std::vector<double> max_part(final_level_no, 0);
            std::vector<double> min_part(final_level_no, 0);
            std::vector<double> level_weight(final_level_no, 0);
            for(int lvl = 0; lvl < final_level_no; lvl++){
                double total_nodes = 0;
                double max_size = 0;
                double min_size = n_;
                for (int w_ptr = final_level_ptr[lvl]; w_ptr < final_level_ptr[lvl + 1]; ++w_ptr) {
                    double cc_size = final_part_ptr[w_ptr + 1] - final_part_ptr[w_ptr];
                    double cc_cost = 0;
                    for(int i = final_part_ptr[w_ptr]; i < final_part_ptr[w_ptr + 1]; i++){
                        int node = final_node_ptr[i];
                        cc_cost += cost[node];
                    }
                    cc_size = cc_cost;
                    if(cc_size > max_size){
                        max_size = cc_size;
                    }
                    if(cc_size < min_size){
                        min_size = cc_size;
                    }
                    total_nodes += cc_size;
                }
                max_part[lvl] = max_size;
                min_part[lvl] = min_size;
                level_weight[lvl] = total_nodes;
            }

            max_diff = 0;
            double total_weight = 0;
            for(int l = 0; l < final_level_no; l++){
                max_diff += level_weight[l] * (max_part[l] - min_part[l]);
                total_weight += level_weight[l];
            }
            max_diff = max_diff / total_weight;
            var = 0;
        }
        void getWaveStatistic(int& nlevels, int& part_no){
            part_no = this->part_no;
            nlevels = this->final_level_no;
        }
        double getSchedulingTime() { return Scheduling_time.elapsed_time;}
        void getStatistic(double& avg_par, double& max_diff, double& var){
            avg_par = this->avg_par;
            max_diff = this->max_diff;
            var = this->var;
        }
        ~SpTrSv_LL_Tree_GLC() = default;
    };
    //================================ SpICh0 GLC UL ================================
    class SpTrSv_LL_Tree_GLC_BFS : public Sptrsv_LL_Serial {
    protected:
        int nthreads;
        timing_measurement Scheduling_time;
        int ngroups;
        std::vector<int> final_level_ptr, final_part_ptr, final_node_ptr;
        int final_level_no;
        std::vector<int> DAG_ptr;
        std::vector<int> DAG_set;
        std::vector<int> group_set, group_ptr;
        std::vector<int> level_ptr, level_set;
        std::vector<int> group_DAG_ptr, group_DAG_set;
        double avgGroupSize;
        int nlevels;
        double avg_par, max_diff, var;
        int part_no = 0;
        std::vector<double> cost;
        bool isLfactor;
        bool bin_pack;
        void build_set() override {
            Scheduling_time.start_timer();
            //Create the DAG
            if(isLfactor){
                std::cout << "Building Tree" << std::endl;
                GLC::buildETree(L1_csr_, nnz_, DAG_ptr, DAG_set);
            } else {
//                computingDAG_CSR(L1_csr_, DAG_ptr_not_prune, DAG_set_not_prune);
                GLC::partialSparsification(n_, nnz_, L1_csc_->p, L1_csc_->i, DAG_ptr, DAG_set);
            }

            GLC::treeBasedGroupingBFS(n_, DAG_ptr, DAG_set,
                                      ngroups, group_ptr, group_set, isLfactor);
            avgGroupSize = 0;
            for(int i = 0; i < ngroups; i++){
                avgGroupSize += (group_ptr[i + 1] - group_ptr[i]);
            }
            avgGroupSize = avgGroupSize / ngroups;
            GLC::buildGroupDAG(n_, ngroups, group_ptr.data(), group_set.data(),
                               DAG_ptr.data(), DAG_set.data(), group_DAG_ptr, group_DAG_set);


            //Computing node Cost
            cost.resize(ngroups, 0);
            auto CSC_Lp = L1_csc_->p;
            auto CSC_Li = L1_csc_->i;
            auto CSR_Lp = L1_csr_->p;
            auto CSR_Li = L1_csr_->i;
            GLC::costComputation(ngroups, CSC_Lp, CSC_Li, CSR_Lp, CSR_Li,
                                 GLC::SpTrSv_LL, group_ptr.data(), group_set.data(), true, cost);

            GLC::GLC(ngroups, group_DAG_ptr[ngroups], group_DAG_ptr, group_DAG_set, cost, nthreads,
                     final_level_no, final_level_ptr, final_part_ptr, final_node_ptr,
                     false, false, bin_pack);

            if(ANALYSE_STATISTIC){
                std::cout << "Static is COMPUTING" << std::endl;
                computeStatistic();
            }

            GLC::ungroupingScheduleAndApplyOrdering(n_, final_level_no, final_level_ptr, final_part_ptr, final_node_ptr,
                                                    group_ptr, group_set,
                                                    DAG_ptr.data(), DAG_set.data(),
                                                    true, false);

            Scheduling_time.measure_elapsed_time();
        }

        timing_measurement fused_code() override {
            timing_measurement t1;
            t1.start_timer();
            sptrsv_csr_lbc(n_, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_,
                           final_level_no, final_level_ptr.data(),
                           final_part_ptr.data(), final_node_ptr.data());
            t1.measure_elapsed_time();
            sym_lib::copy_vector(0, n_, x_in_, x_);
            return t1;
        }

    public:
        /*
         * @brief Class constructor
         * @param L Sparse Matrix with CSR format
         * @param L_csc Sparse Matrix with CSC format
         * @param correct_x The right answer for x
         * @param name The name of the algorithm
         * @param nt number of threads
         */
        SpTrSv_LL_Tree_GLC_BFS(CSR *L, CSC *L_csc,
                           double *correct_x, std::string name,
                           int nt, bool isLfactor, bool bin_pack)
                           : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
            L1_csr_ = L;
            L1_csc_ = L_csc;
            correct_x_ = correct_x;
            nthreads = nt;
            this->isLfactor = isLfactor;
            this->bin_pack = bin_pack;
        };

        void computeCharacteristics(double& first, double& second, double& third, double& avg){

            std::vector<double> max_array(final_level_no, 0);
            std::vector<double> whole_array(final_level_no, 0);
            for(int lvl = 0; lvl < final_level_no; lvl++){
                double total_nodes = 0;
                double max_size = 0;
                for (int w_ptr = final_level_ptr[lvl]; w_ptr < final_level_ptr[lvl + 1]; ++w_ptr) {
                    int cc_size = final_part_ptr[w_ptr + 1] - final_part_ptr[w_ptr];
                    if(cc_size > max_size){
                        max_size = cc_size;
                    }
                    total_nodes += cc_size;
                }
                if(lvl == 0){
                    first = max_size / total_nodes;
                }
                if(lvl == 1){
                    second = max_size / total_nodes;
                }
                if(lvl == 2){
                    third = max_size / total_nodes;
                }
                max_array[lvl] = max_size;
                whole_array[lvl] = total_nodes;
            }

            double max_cc_sizes = 0;
            for(auto& iter: max_array){
                max_cc_sizes += iter;
            }
            avg = max_cc_sizes / n_;

            #ifndef NDEBUG
            double sum = 0;
            for(auto &iter: whole_array){
                sum += iter;
            }
            assert(sum == n_);
            #endif
        }

        void computeStatistic(){
            double part_no = 0;
            for(int i = 0; i < final_part_ptr.size() - 1; i++){
                if(final_part_ptr[i + 1] != final_part_ptr[i]){
                    part_no++;
                }
            }
            this->part_no = part_no;
            avg_par = part_no/ final_level_no;
            std::vector<double> max_part(final_level_no, 0);
            std::vector<double> min_part(final_level_no, 0);
            std::vector<double> level_weight(final_level_no, 0);
            for(int lvl = 0; lvl < final_level_no; lvl++){
                double total_nodes = 0;
                double max_size = 0;
                double min_size = n_;
                for (int w_ptr = final_level_ptr[lvl]; w_ptr < final_level_ptr[lvl + 1]; ++w_ptr) {
                    double cc_size = final_part_ptr[w_ptr + 1] - final_part_ptr[w_ptr];
                    double cc_cost = 0;
                    for(int i = final_part_ptr[w_ptr]; i < final_part_ptr[w_ptr + 1]; i++){
                        int node = final_node_ptr[i];
                        cc_cost += cost[node];
                    }
                    cc_size = cc_cost;
                    if(cc_size > max_size){
                        max_size = cc_size;
                    }
                    if(cc_size < min_size){
                        min_size = cc_size;
                    }
                    total_nodes += cc_size;
                }
                max_part[lvl] = max_size;
                min_part[lvl] = min_size;
                level_weight[lvl] = total_nodes;
            }

            max_diff = 0;
            double total_weight = 0;
            for(int l = 0; l < final_level_no; l++){
                max_diff += level_weight[l] * (max_part[l] - min_part[l]);
                total_weight += level_weight[l];
            }
            max_diff = max_diff / total_weight;
            var = 0;
        }
        double getSchedulingTime() { return Scheduling_time.elapsed_time;}
        void getStatistic(double& avg_par, double& max_diff, double& var){
            avg_par = this->avg_par;
            max_diff = this->max_diff;
            var = this->var;
        }
        void getWaveStatistic(int& nlevels, int& part_no){
            part_no = this->part_no;
            nlevels = this->final_level_no;
        }
        double getPotentialGain(){
            double critical_path = 0;
            std::vector<double> cost(n_, 0);
            std::vector<double> cost_per_level(final_level_no, 0);
            std::vector<double> part_per_level(final_level_no, 0);
            std::vector<double> node_per_level(final_level_no, 0);
            std::vector<double> max_cost(final_level_no, 0);
            auto CSC_Lp = L1_csc_->p;
            auto CSC_Li = L1_csc_->i;
            auto CSR_Lp = L1_csr_->p;
            auto CSR_Li = L1_csr_->i;
            GLC::costComputation(n_, CSC_Lp, CSC_Li, CSR_Lp, CSR_Li,
                                 GLC::SpTrSv_LL, NULL, NULL, false, cost);
            double Total_cost = 0;
            double total_part_number = 0;

            for(int lvl = 0; lvl < final_level_no; lvl++){
                for(int prt_ptr = final_level_ptr[lvl]; prt_ptr < final_level_ptr[lvl + 1]; prt_ptr++){
                    if(final_part_ptr[prt_ptr+1] != final_part_ptr[prt_ptr]){
                        total_part_number++;
                        part_per_level[lvl]++;
                    }
                    double part_cost = 0;
                    for(int node_ptr = final_part_ptr[prt_ptr]; node_ptr < final_part_ptr[prt_ptr+1]; node_ptr++){
                        int node = final_node_ptr[node_ptr];
                        assert(!isnan(cost[node]));
                        Total_cost += cost[node];
                        cost_per_level[lvl] += cost[node];
                        part_cost += cost[node];
                        node_per_level[lvl]++;
                    }
                    if(part_cost > max_cost[lvl]){
                        max_cost[lvl] = part_cost;
                    }
                }
                critical_path += max_cost[lvl];
            }

            return 1 - (Total_cost / nthreads) / critical_path;
        }

        void computeMetrics(double& avgParallelism, double& WeightedAvgParallelism, double& CostAvgParallelism,
                            double& AvgGroupSize,
                            double& AvgCostPerVertexSpTrSv, double& AvgCostPerNNZSpTrSv,
                            double& AvgCostPerVertexSpILU0, double& AvgCostPerNNZSpLU0,
                            double& AvgCostPerVertexSpIC0, double& AvgCostPerNNZSpIC0,
                            double& PG){
            double Total_cost = 0;
            double total_part_number = 0;
            std::vector<double> cost_per_level(final_level_no, 0);
            std::vector<double> part_per_level(final_level_no, 0);
            std::vector<double> node_per_level(final_level_no, 0);
            std::vector<double> max_cost(final_level_no, 0);
            double critical_path = 0;
            std::vector<double> cost;
            cost.resize(n_, 0);
            auto CSC_Lp = L1_csc_->p;
            auto CSC_Li = L1_csc_->i;
            auto CSR_Lp = L1_csr_->p;
            auto CSR_Li = L1_csr_->i;
            GLC::costComputation(n_, CSC_Lp, CSC_Li, CSR_Lp, CSR_Li,
                                 GLC::SpTrSv_LL, NULL, NULL, false, cost);
            for(int lvl = 0; lvl < final_level_no; lvl++){
                for(int prt_ptr = final_level_ptr[lvl]; prt_ptr < final_level_ptr[lvl + 1]; prt_ptr++){
                    if(final_part_ptr[prt_ptr+1] != final_part_ptr[prt_ptr]){
                        total_part_number++;
                        part_per_level[lvl]++;
                    }
                    double part_cost = 0;
                    for(int node_ptr = final_part_ptr[prt_ptr]; node_ptr < final_part_ptr[prt_ptr+1]; node_ptr++){
                        int node = final_node_ptr[node_ptr];
                        assert(!isnan(cost[node]));
                        Total_cost += cost[node];
                        cost_per_level[lvl] += cost[node];
                        part_cost += cost[node];
                        node_per_level[lvl]++;
                    }
                    if(part_cost > max_cost[lvl]){
                        max_cost[lvl] = part_cost;
                    }
                }
                critical_path += max_cost[lvl];
            }

            PG = 1 - (Total_cost / nthreads) / critical_path;

            avgParallelism = total_part_number / final_level_no;
            for(int lvl = 0; lvl < final_level_no; lvl++){
                CostAvgParallelism += (part_per_level[lvl] * cost_per_level[lvl]);
                WeightedAvgParallelism += (part_per_level[lvl] * node_per_level[lvl]);
            }

            CostAvgParallelism = CostAvgParallelism / Total_cost;
            WeightedAvgParallelism = WeightedAvgParallelism / n_;
            AvgGroupSize = this->avgGroupSize;


            AvgCostPerVertexSpTrSv = Total_cost / n_;
            AvgCostPerNNZSpTrSv = Total_cost / nnz_;

            cost.clear();
            cost.resize(n_, 0);
            GLC::costComputation(n_, CSC_Lp, CSC_Li, CSR_Lp, CSR_Li,
                                 GLC::SpICh0_UL, NULL, NULL, false, cost);

            Total_cost = 0;
            for(auto &iter: cost){
                Total_cost += iter;
            }
            AvgCostPerVertexSpIC0 = Total_cost / n_;
            AvgCostPerNNZSpIC0 = Total_cost / nnz_;

            auto tmp = make_full(L1_csc_);
            auto CSR_A = csc_to_csr(tmp);
            delete tmp;
            cost.clear();
            cost.resize(n_, 0);
            GLC::costComputation(n_, NULL, NULL, CSR_A->p, CSR_A->i,
                                 GLC::SpICh0_UL, NULL, NULL, false, cost);

            Total_cost = 0;
            for(auto &iter: cost){
                Total_cost += iter;
            }
            AvgCostPerVertexSpILU0 = Total_cost / n_;
            AvgCostPerNNZSpLU0 = Total_cost / nnz_;
            delete CSR_A;
        }
        ~SpTrSv_LL_Tree_GLC_BFS() = default;
    };
    //================================ SpICh0 Tree BFS GLC UL ================================
    class SpTrSv_LL_Tree_GLC_BFS_P2P : public Sptrsv_LL_Serial {
    protected:
        int nthreads;
        timing_measurement Scheduling_time;
        int ngroups;
        std::vector<int> final_level_ptr, final_part_ptr, final_node_ptr;
        int final_level_no;
        std::vector<int> DAG_ptr;
        std::vector<int> DAG_set;
        std::vector<int> group_set, group_ptr;
        std::vector<int> level_ptr, level_set;
        std::vector<int> group_DAG_ptr, group_DAG_set;
        int nlevels;
        double avg_par, max_diff, var;
        int part_no = 0;
        std::vector<double> cost;
        bool isLfactor;
        bool bin_pack;
        SpMP::LevelSchedule *p2pScheduleWithTransitiveReduction;
        const int* invPerm;
        bool isfirst_time;

        void build_set() override {
            Scheduling_time.start_timer();

            if(isLfactor){
                std::cout << "Building Tree" << std::endl;
                GLC::buildETree(L1_csr_, nnz_, DAG_ptr, DAG_set);
            } else {
                //                computingDAG_CSR(L1_csr_, DAG_ptr_not_prune, DAG_set_not_prune);
                GLC::partialSparsification(n_, nnz_, L1_csc_->p, L1_csc_->i, DAG_ptr, DAG_set);
            }

            GLC::treeBasedGroupingBFS(n_, DAG_ptr, DAG_set,
                                      ngroups, group_ptr, group_set, isLfactor);

            GLC::buildGroupDAG(n_, ngroups, group_ptr.data(), group_set.data(), DAG_ptr.data(), DAG_set.data(), group_DAG_ptr, group_DAG_set);


            //Computing node Cost
            cost.resize(ngroups, 0);
            auto CSC_Lp = L1_csc_->p;
            auto CSC_Li = L1_csc_->i;
            auto CSR_Lp = L1_csr_->p;
            auto CSR_Li = L1_csr_->i;
            GLC::costComputation(ngroups, CSC_Lp, CSC_Li, CSR_Lp, CSR_Li,
                                 GLC::SpTrSv_LL, group_ptr.data(), group_set.data(), true, cost);

            GLC::GLC(ngroups, group_DAG_ptr[ngroups], group_DAG_ptr, group_DAG_set, cost, nthreads,
                     final_level_no, final_level_ptr, final_part_ptr, final_node_ptr,
                     false, false, bin_pack);

            if(ANALYSE_STATISTIC){
                std::cout << "Static is COMPUTING" << std::endl;
                computeStatistic();
            }

            GLC::ungroupingScheduleAndApplyOrdering(n_, final_level_no, final_level_ptr, final_part_ptr, final_node_ptr,
                                                    group_ptr, group_set,
                                                    DAG_ptr.data(), DAG_set.data(),
                                                    true, false);

            //Create Group set
            GLC::getFinalScheduleDAG(n_, final_level_no, final_level_ptr.data(),
                                     final_part_ptr.data(), final_node_ptr.data(), L1_csc_->p, L1_csc_->i,
                                     ngroups, group_ptr, group_set,
                                     DAG_ptr, DAG_set);
            if(!isfirst_time){
                delete p2pScheduleWithTransitiveReduction;
            }

            //Create DAG
            CSC* CSC_DAG = new CSC(DAG_ptr.size() - 1, DAG_ptr.size() - 1, DAG_set.size());
            std::copy(DAG_ptr.begin(), DAG_ptr.end(), CSC_DAG->p);
            std::copy(DAG_set.begin(), DAG_set.end(), CSC_DAG->i);
            std::fill(CSC_DAG->x, CSC_DAG->x + DAG_set.size(), 1);
            CSC_DAG->stype = -1;
            auto tmp = make_full(CSC_DAG);
            auto CSR_A = csc_to_csr(tmp);
            delete tmp;
            delete CSC_DAG;
            //Create SpMP schedule
            //Convert lbc CSR version to SpMP CSR version
            auto* A = new SpMP::CSR();
            Convert_LBCCSR_to_SpMP(CSR_A, A);

            /////////////////////////////////////////////////////////////////////////////
            // Construct schedules
            /////////////////////////////////////////////////////////////////////////////
            bool hasZeroDiag = A->hasZeroDiag();
            if (hasZeroDiag) {
                fprintf(stderr, "Matrix has a zero diagonal element. Can't do Gauss Seidel\n");
            }

            int *symRowPtr = NULL, *symColIdx = NULL, *symDiagPtr = NULL, *symExtPtr = NULL;
            bool wasSymmetric = getSymmetricNnzPattern(A, &symRowPtr, &symDiagPtr, &symExtPtr, &symColIdx);

            p2pScheduleWithTransitiveReduction = new SpMP::LevelSchedule;
            if (wasSymmetric) {
                p2pScheduleWithTransitiveReduction->constructTaskGraph(*A);
            }
            else {
                p2pScheduleWithTransitiveReduction->constructTaskGraph(
                        A->m, symRowPtr, symDiagPtr, symExtPtr, symColIdx,
                        SpMP::PrefixSumCostFunction(symRowPtr));

                FREE(symRowPtr);
                FREE(symColIdx);
                FREE(symDiagPtr);
                FREE(symExtPtr);
            }

            /////////////////////////////////////////////////////////////////////////////
            // Reorder matrix
            /////////////////////////////////////////////////////////////////////////////
            const int *perm = p2pScheduleWithTransitiveReduction->origToThreadContPerm;
            const int *invPerm = p2pScheduleWithTransitiveReduction->threadContToOrigPerm;
            assert(SpMP::isPerm(perm, A->m));
            assert(SpMP::isPerm(invPerm, A->m));
            this->invPerm = invPerm;

            delete A;
            delete CSR_A;
            isfirst_time = false;

            Scheduling_time.measure_elapsed_time();
        }

        timing_measurement fused_code() override {
            timing_measurement t1;
            auto Lp = L1_csr_->p;
            auto Li = L1_csr_->i;
            auto Lx = L1_csr_->x;
            const SpMP::LevelSchedule &schedule = *p2pScheduleWithTransitiveReduction;
            const int* perm = invPerm;
            double* x = x_in_;
            t1.start_timer();
            #pragma omp parallel
            {
                int nthreads = omp_get_num_threads();
                int tid = omp_get_thread_num();

                const int ntasks = schedule.ntasks;
                const short *nparents = schedule.nparentsForward;
                const std::vector<int>& threadBoundaries = schedule.threadBoundaries;
                const std::vector<int>& taskBoundaries = schedule.taskBoundaries;

                int nPerThread = (ntasks + nthreads - 1)/nthreads;
                int nBegin = std::min(nPerThread*tid, ntasks);
                int nEnd = std::min(nBegin + nPerThread, ntasks);

                volatile int *taskFinished = schedule.taskFinished;
                int **parents = schedule.parentsForward;

                memset((char *)(taskFinished + nBegin), 0, (nEnd - nBegin)*sizeof(int));

                synk::Barrier::getInstance()->wait(tid);

                for (int task = threadBoundaries[tid]; task < threadBoundaries[tid + 1]; ++task) {
                    SPMP_LEVEL_SCHEDULE_WAIT;

                    for (int c = taskBoundaries[task]; c < taskBoundaries[task + 1]; ++c) {
                        int group = perm[c];
                        for(int g_ptr = group_ptr[group]; g_ptr < group_ptr[group + 1]; g_ptr++){
                            int i = group_set[g_ptr];
                            for (int j = Lp[i]; j < Lp[i + 1] - 1; j++)
                            {
                                x[i] -= Lx[j] * x[Li[j]];
                            }
                            x[i] /= Lx[Lp[i + 1] - 1];
                        }

                    }
                    SPMP_LEVEL_SCHEDULE_NOTIFY;
                } // for each task
            } // omp parallel
            t1.measure_elapsed_time();
            sym_lib::copy_vector(0, n_, x_in_, x_);
            return t1;
        }

    public:
        /*
         * @brief Class constructor
         * @param L Sparse Matrix with CSR format
         * @param L_csc Sparse Matrix with CSC format
         * @param correct_x The right answer for x
         * @param name The name of the algorithm
         * @param nt number of threads
         */
        SpTrSv_LL_Tree_GLC_BFS_P2P(CSR *L, CSC *L_csc,
                               double *correct_x, std::string name,
                               int nt, bool isLfactor, bool bin_pack)
                               : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
            L1_csr_ = L;
            L1_csc_ = L_csc;
            correct_x_ = correct_x;
            nthreads = nt;
            this->isLfactor = isLfactor;
            isfirst_time = true;
            this->bin_pack = bin_pack;
        };

        void computeCharacteristics(double& first, double& second, double& third, double& avg){

            std::vector<double> max_array(final_level_no, 0);
            std::vector<double> whole_array(final_level_no, 0);
            for(int lvl = 0; lvl < final_level_no; lvl++){
                double total_nodes = 0;
                double max_size = 0;
                for (int w_ptr = final_level_ptr[lvl]; w_ptr < final_level_ptr[lvl + 1]; ++w_ptr) {
                    int cc_size = final_part_ptr[w_ptr + 1] - final_part_ptr[w_ptr];
                    if(cc_size > max_size){
                        max_size = cc_size;
                    }
                    total_nodes += cc_size;
                }
                if(lvl == 0){
                    first = max_size / total_nodes;
                }
                if(lvl == 1){
                    second = max_size / total_nodes;
                }
                if(lvl == 2){
                    third = max_size / total_nodes;
                }
                max_array[lvl] = max_size;
                whole_array[lvl] = total_nodes;
            }

            double max_cc_sizes = 0;
            for(auto& iter: max_array){
                max_cc_sizes += iter;
            }
            avg = max_cc_sizes / n_;

            #ifndef NDEBUG
            double sum = 0;
            for(auto &iter: whole_array){
                sum += iter;
            }
            assert(sum == n_);
            #endif
        }

        void computeStatistic(){
            double part_no = 0;
            for(int i = 0; i < final_part_ptr.size() - 1; i++){
                if(final_part_ptr[i + 1] != final_part_ptr[i]){
                    part_no++;
                }
            }
            this->part_no = part_no;
            avg_par = part_no/ final_level_no;
            std::vector<double> max_part(final_level_no, 0);
            std::vector<double> min_part(final_level_no, 0);
            std::vector<double> level_weight(final_level_no, 0);
            for(int lvl = 0; lvl < final_level_no; lvl++){
                double total_nodes = 0;
                double max_size = 0;
                double min_size = n_;
                for (int w_ptr = final_level_ptr[lvl]; w_ptr < final_level_ptr[lvl + 1]; ++w_ptr) {
                    double cc_size = final_part_ptr[w_ptr + 1] - final_part_ptr[w_ptr];
                    double cc_cost = 0;
                    for(int i = final_part_ptr[w_ptr]; i < final_part_ptr[w_ptr + 1]; i++){
                        int node = final_node_ptr[i];
                        cc_cost += cost[node];
                    }
                    cc_size = cc_cost;
                    if(cc_size > max_size){
                        max_size = cc_size;
                    }
                    if(cc_size < min_size){
                        min_size = cc_size;
                    }
                    total_nodes += cc_size;
                }
                max_part[lvl] = max_size;
                min_part[lvl] = min_size;
                level_weight[lvl] = total_nodes;
            }

            max_diff = 0;
            double total_weight = 0;
            for(int l = 0; l < final_level_no; l++){
                max_diff += level_weight[l] * (max_part[l] - min_part[l]);
                total_weight += level_weight[l];
            }
            max_diff = max_diff / total_weight;
            var = 0;
        }
        double getSchedulingTime() { return Scheduling_time.elapsed_time;}
        void getStatistic(double& avg_par, double& max_diff, double& var){
            avg_par = this->avg_par;
            max_diff = this->max_diff;
            var = this->var;
        }
        void getWaveStatistic(int& nlevels, int& part_no){
            part_no = this->part_no;
            nlevels = this->final_level_no;
        }

        ~SpTrSv_LL_Tree_GLC_BFS_P2P(){
            delete p2pScheduleWithTransitiveReduction;
        };
    };
    //================================ SpICh0 GLC UL ================================
    class SpTrSv_LL_Tree_GLC_DFS : public Sptrsv_LL_Serial {
    protected:
        int nthreads;
        timing_measurement Scheduling_time;
        int ngroups;
        std::vector<int> final_level_ptr, final_part_ptr, final_node_ptr;
        int final_level_no;
        std::vector<int> DAG_ptr;
        std::vector<int> DAG_set;
        std::vector<int> group_set, group_ptr;
        std::vector<int> level_ptr, level_set;
        std::vector<int> group_DAG_ptr, group_DAG_set;
        int nlevels;
        bool isLfactor;
        void build_set() override {
            Scheduling_time.start_timer();
            //Create the DAG
            std::vector<int> DAG_ptr_not_prune;
            std::vector<int> DAG_set_not_prune;
            this->isLfactor=false;
            if(isLfactor){
                GLC::buildETree(L1_csr_, nnz_, DAG_ptr, DAG_set);
            } else {
//                computingDAG_CSR(L1_csr_, DAG_ptr_not_prune, DAG_set_not_prune);
                GLC::partialSparsification(n_, nnz_, L1_csc_->p, L1_csc_->i, DAG_ptr, DAG_set);
            }

            GLC::chainGrouping(n_, DAG_ptr, DAG_set,
                                      ngroups, group_ptr, group_set);

            GLC::buildGroupDAG(n_, ngroups, group_ptr.data(), group_set.data(), DAG_ptr.data(), DAG_set.data(), group_DAG_ptr, group_DAG_set);


            //Computing node Cost
            std::vector<double> cost(ngroups, 0);
            auto CSC_Lp = L1_csc_->p;
            auto CSC_Li = L1_csc_->i;
            auto CSR_Lp = L1_csr_->p;
            auto CSR_Li = L1_csr_->i;
            GLC::costComputation(ngroups, CSC_Lp, CSC_Li, CSR_Lp, CSR_Li,
                                 GLC::SpTrSv_LL, group_ptr.data(), group_set.data(), true, cost);

            GLC::GLC(ngroups, group_DAG_ptr[ngroups], group_DAG_ptr, group_DAG_set, cost, nthreads,
                     final_level_no, final_level_ptr, final_part_ptr, final_node_ptr,false, false);

            GLC::ungroupingScheduleAndApplyOrdering(n_, final_level_no, final_level_ptr, final_part_ptr, final_node_ptr,
                                                    group_ptr, group_set, DAG_ptr.data(), DAG_set.data(),
                                                    true, false);

            Scheduling_time.measure_elapsed_time();
        }

        timing_measurement fused_code() override {
            timing_measurement t1;
            t1.start_timer();
            sptrsv_csr_lbc(n_, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_,
                           final_level_no, final_level_ptr.data(),
                           final_part_ptr.data(), final_node_ptr.data());
            t1.measure_elapsed_time();
            sym_lib::copy_vector(0, n_, x_in_, x_);
            return t1;
        }

    public:
        /*
         * @brief Class constructor
         * @param L Sparse Matrix with CSR format
         * @param L_csc Sparse Matrix with CSC format
         * @param correct_x The right answer for x
         * @param name The name of the algorithm
         * @param nt number of threads
         */
        SpTrSv_LL_Tree_GLC_DFS(CSR *L, CSC *L_csc,
                               double *correct_x, std::string name,
                               int nt, bool isLfactor)
                               : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
            L1_csr_ = L;
            L1_csc_ = L_csc;
            correct_x_ = correct_x;
            nthreads = nt;
            this->isLfactor = isLfactor;
        };

        void computeCharacteristics(double& first, double& second, double& third, double& avg){

            std::vector<double> max_array(final_level_no, 0);
            std::vector<double> whole_array(final_level_no, 0);
            for(int lvl = 0; lvl < final_level_no; lvl++){
                double total_nodes = 0;
                double max_size = 0;
                for (int w_ptr = final_level_ptr[lvl]; w_ptr < final_level_ptr[lvl + 1]; ++w_ptr) {
                    int cc_size = final_part_ptr[w_ptr + 1] - final_part_ptr[w_ptr];
                    if(cc_size > max_size){
                        max_size = cc_size;
                    }
                    total_nodes += cc_size;
                }
                if(lvl == 0){
                    first = max_size / total_nodes;
                }
                if(lvl == 1){
                    second = max_size / total_nodes;
                }
                if(lvl == 2){
                    third = max_size / total_nodes;
                }
                max_array[lvl] = max_size;
                whole_array[lvl] = total_nodes;
            }

            double max_cc_sizes = 0;
            for(auto& iter: max_array){
                max_cc_sizes += iter;
            }
            avg = max_cc_sizes / n_;

            #ifndef NDEBUG
            double sum = 0;
            for(auto &iter: whole_array){
                sum += iter;
            }
            assert(sum == n_);
            #endif
        }

        void groupAnalysis(){
            //Average Group size, max and min
            double max_size = 0;
            double min_size = n_;
            double avg_size = 0;
            for(int i = 0; i < ngroups; i++){
                int group_size = group_ptr[i + 1] - group_ptr[i];
                if(group_size > max_size){
                    max_size = group_size;
                }
                if(group_size < min_size){
                    min_size = group_size;
                }
                avg_size += group_size;
            }

            computingLevelSet_CSC(ngroups, group_DAG_ptr.data(), group_DAG_set.data(), level_ptr, level_set, nlevels);
            std::cout << "The max size is: " << max_size << " The min size is: " <<
            min_size << " The avg size is: " << avg_size / ngroups << std::endl;


            //Reusing Ratio
            //Computing node_to_lvl
            std::vector<int> node_to_lvl(n_, 0);
            for(int i = 0; i < nlevels; i++){
                for(int j = level_ptr[i]; j < level_ptr[i + 1]; j++){
                    int group = level_set[j];
                    for(int g = group_ptr[group]; g < group_ptr[group + 1]; g++){
                        int node = group_set[g];
                        node_to_lvl[node] = i;
                    }
                }
            }

            double reuse_edge = 0;
            double edge = 0;
            //Computing the reusing (edges inside a group which are inside a level);
            for(int i = 0; i < nlevels; i++){
                for(int j = level_ptr[i]; j < level_ptr[i + 1]; j++){
                    int group = level_set[j];
                    for(int g = group_ptr[group]; g < group_ptr[group + 1]; g++){
                        int node = group_set[g];
                        for(int child_ptr = L1_csc_->p[node] + 1; child_ptr < L1_csc_->p[node + 1]; child_ptr++){
                            int child = L1_csc_->i[child_ptr];
                            if(node_to_lvl[child] == node_to_lvl[node]){
                                reuse_edge++;
                            }
                            edge++;
                        }
                    }
                }
            }
            assert(edge == nnz_ - n_);
            std::cout << "The reuse ratio is: " << reuse_edge / edge << std::endl;

        }
        void getReuseRatio(){
            //Compute the node 2 level based on l-paritions
            std::vector<int> node_to_l_partition(n_, 0);
            for(int lvl = 0; lvl < final_level_no; lvl++){
                for(int part_ptr = final_level_ptr[lvl]; part_ptr < final_level_ptr[lvl + 1]; part_ptr++){
                    for(int node_ptr = final_part_ptr[part_ptr]; node_ptr < final_part_ptr[part_ptr + 1]; node_ptr++){
                        int node = final_node_ptr[node_ptr];
                        node_to_l_partition[node] = lvl;
                    }
                }
            }

            //Compute the reuse ratio as the number of edges inside a l-partition / total number of edges;
            auto graph_ptr = L1_csc_->p;
            auto graph_set = L1_csc_->i;
            double total_edges = 0, reuse_edges = 0;
            for(int node = 0; node < n_; node++){
                for(int child_ptr = graph_ptr[node] + 1; child_ptr < graph_ptr[node + 1]; child_ptr++){
                    int child = graph_set[child_ptr];
                    if(node_to_l_partition[child] == node_to_l_partition[node]){
                        reuse_edges++;
                    }
                    total_edges++;
                }
            }
            std::cout << "The reuse ratio is: " <<  reuse_edges / total_edges << std::endl;
        }
        void memoryFootPrint(){
            std::vector<long long> lvl_incomming_edges(final_level_no, 0);
            for(int lvl = 0; lvl < final_level_no; lvl++){
                for(int part_ptr = final_level_ptr[lvl]; part_ptr < final_level_ptr[lvl + 1]; part_ptr++){
                    for(int node_ptr = final_part_ptr[part_ptr]; node_ptr < final_part_ptr[part_ptr + 1]; node_ptr++){
                        int node = final_node_ptr[node_ptr];
                        lvl_incomming_edges[lvl] += (2 * (L1_csr_->p[node + 1] - L1_csr_->p[node]));
                    }
                }
            }

        }
        double getSchedulingTime() { return Scheduling_time.elapsed_time;}

        ~SpTrSv_LL_Tree_GLC_DFS() = default;
    };
    //================================ SpTrSv Block Vectorized LBC LL ================================
    class SpTrSv_LL_Blocked_GLC : public Sptrsv_LL_Serial {
    protected:
        int nthreads;
        timing_measurement Scheduling_time;
        int num_supernodes, num_dep_edges;
        std::vector<int> supernode_ptr;
        std::vector<int> final_level_ptr, final_part_ptr, final_node_ptr;
        int final_level_no;
        std::vector<int> DAG_ptr;
        std::vector<int> DAG_set;
        std::vector<double> cost;
        supernodaltools::SuperNodal* supernode_obj;
        void build_set() override {
            Scheduling_time.start_timer();
            supernodaltools::SuperNodal super_node_obj(supernodaltools::CSR_TYPE,
                                                       L1_csc_, L1_csr_, L1_csc_, 0);

            super_node_obj.getDependencyDAGCSCformat(DAG_ptr, DAG_set,num_supernodes, num_dep_edges);

            cost.resize(num_supernodes, 0);
            super_node_obj.getSuperNodeGroups(supernode_ptr);
            #pragma omp parallel for
            for (int super = 0; super < num_supernodes; ++super) {
                for(int i = supernode_ptr[super]; i < supernode_ptr[super + 1]; i++){
                    cost[super] += L1_csr_->p[i + 1] - L1_csr_->p[i];
                }
            }

            GLC::GLC(num_supernodes, num_dep_edges, DAG_ptr, DAG_set, cost, nthreads,
                     final_level_no, final_level_ptr, final_part_ptr, final_node_ptr,
                     true, false, false);
            Scheduling_time.measure_elapsed_time();
        }

        timing_measurement fused_code() override {
            //sym_lib::rhs_init(L1_csc_->n, L1_csc_->p, L1_csc_->i, L1_csc_->x, x_); // x is b
            timing_measurement t1;
            auto Lp = L1_csr_->p;
            auto Li = L1_csr_->i;
            auto Lx = L1_csr_->x;
            auto x = x_in_;

            t1.start_timer();
            #pragma omp parallel
            {
                std::vector<double> tempvec(n_, 0);
                for (int lvl = 0; lvl < final_level_no; ++lvl) {
                    #pragma omp for
                    for (int w_ptr = final_level_ptr[lvl]; w_ptr < final_level_ptr[lvl + 1]; ++w_ptr) {
                        for(int super_ptr = final_part_ptr[w_ptr]; super_ptr < final_part_ptr[w_ptr + 1]; super_ptr++){
                            int super = final_node_ptr[super_ptr];
                            int start_row = supernode_ptr[super];
                            int end_row = supernode_ptr[super + 1];
                            //ncol is the number of columns in the off-diagonal block
                            int num_off_diag_col = Lp[start_row + 1] - Lp[start_row] - 1;
                            assert(num_off_diag_col >= 0);
                            int nrows = end_row - start_row;
                            assert(nrows > 0);
                            //Solving the independent part
                            //Copy x[Li[col]] into a continues buffer
                            for (int col_ptr = Lp[supernode_ptr[super]], k = 0;
                            col_ptr < Lp[supernode_ptr[super]] + num_off_diag_col; col_ptr++, k++) {
                                tempvec[k] = x[Li[col_ptr]];
                            }
                            custom_blas::SpTrSv_MatVecCSR_BLAS(nrows, num_off_diag_col, &Lx[Lp[start_row]], tempvec.data(), &x[start_row]);
                            custom_blas::SpTrSv_LSolveCSR_BLAS(nrows, num_off_diag_col, &Lx[Lp[start_row] + num_off_diag_col], &x[start_row]);
                        }
                    }
                }
            }
            t1.measure_elapsed_time();
            sym_lib::copy_vector(0, n_, x_in_, x_);
            return t1;
        }

    public:
        /*
         * @brief Class constructor
         * @param L Sparse Matrix with CSR format
         * @param L_csc Sparse Matrix with CSC format
         * @param correct_x The right answer for x
         * @param name The name of the algorithm
         * @param nt number of threads
         */
        SpTrSv_LL_Blocked_GLC(CSR *L, CSC *L_csc,
                                      double *correct_x, std::string name,
                                      int nt, bool sort=true)
                                      : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
            L1_csr_ = L;
            L1_csc_ = L_csc;
            correct_x_ = correct_x;
            nthreads = nt;
        };

        double getSchedulingTime() { return Scheduling_time.elapsed_time; }

        ~SpTrSv_LL_Blocked_GLC() = default;
    };
    //================================ SpTrSv Block Vectorized LBC LL ================================
    class SpTrSv_LL_Blocked_Grouped_GLC : public Sptrsv_LL_Serial {
    protected:
        int nthreads;
        timing_measurement Scheduling_time;
        int num_supernodes, num_dep_edges;
        std::vector<int> supernode_ptr;
        std::vector<int> final_level_ptr, final_part_ptr, final_node_ptr;
        std::vector<int> groupPtr, groupSet;
        int ngroups;
        int final_level_no;
        std::vector<int> DAG_ptr;
        std::vector<int> DAG_set;
        std::vector<double> cost;
        supernodaltools::SuperNodal* supernode_obj;
        void build_set() override {
            Scheduling_time.start_timer();
            supernodaltools::SuperNodal super_node_obj(supernodaltools::CSR_TYPE,
                                                       L1_csc_, L1_csr_, L1_csc_, 0);

            super_node_obj.getDependencyDAGCSCformat(DAG_ptr, DAG_set,num_supernodes, num_dep_edges);

            CSC DAG(num_supernodes,num_supernodes, num_dep_edges);
            std::copy(DAG_ptr.begin(), DAG_ptr.end(), DAG.p);
            std::copy(DAG_set.begin(), DAG_set.end(), DAG.i);
            auto DAG_inv = csc_to_csr(&DAG);

            std::vector<int> DAG_inv_ptr(num_supernodes + 1);
            std::vector<int> DAG_inv_set(num_dep_edges);
            std::copy(DAG_inv->p, DAG_inv->p + num_supernodes + 1, DAG_inv_ptr.begin());
            std::copy(DAG_inv->i, DAG_inv->i + num_dep_edges, DAG_inv_set.begin());

            std::vector<int> group_DAG_ptr;
            std::vector<int> group_DAG_set;
            //Grouped DAG
            computingGroupedDAG_CSR(num_supernodes, DAG_inv_ptr, DAG_inv_set,
                                    group_DAG_ptr, group_DAG_set, groupPtr, groupSet, ngroups);

            cost.resize(ngroups, 0);
            super_node_obj.getSuperNodeGroups(supernode_ptr);
            #pragma omp parallel for
            for(int g = 0; g < ngroups; g++){
                for(int node_ptr = groupPtr[g]; node_ptr < groupPtr[g + 1]; node_ptr++){
                    int node = groupSet[node_ptr];
                    for(int i = supernode_ptr[node]; i < supernode_ptr[node + 1]; i++){
                        cost[g] += L1_csr_->p[i + 1] - L1_csr_->p[i];
                    }
                }
            }

            assert(L1_csc_ != nullptr);
            assert(L1_csr_ != nullptr);
            GLC::GLC(ngroups, group_DAG_ptr[ngroups], group_DAG_ptr, group_DAG_set, cost, nthreads,
                            final_level_no, final_level_ptr, final_part_ptr, final_node_ptr,
                            true, false, false);
            Scheduling_time.measure_elapsed_time();
        }

        timing_measurement fused_code() override {
            //sym_lib::rhs_init(L1_csc_->n, L1_csc_->p, L1_csc_->i, L1_csc_->x, x_); // x is b
            timing_measurement t1;
            auto Lp = L1_csr_->p;
            auto Li = L1_csr_->i;
            auto Lx = L1_csr_->x;
            auto x = x_in_;

            t1.start_timer();
            #pragma omp parallel
            {
                std::vector<double> tempvec(n_, 0);
                for (int lvl = 0; lvl < final_level_no; ++lvl) {
                    #pragma omp for
                    for (int w_ptr = final_level_ptr[lvl]; w_ptr < final_level_ptr[lvl + 1]; ++w_ptr) {
                        for(int group_ptr = final_part_ptr[w_ptr]; group_ptr < final_part_ptr[w_ptr + 1]; group_ptr++){
                            int group = final_node_ptr[group_ptr];
                            for(int super_ptr = groupPtr[group]; super_ptr < groupPtr[group + 1]; super_ptr++){
                                int super = groupSet[super_ptr];
                                int start_row = supernode_ptr[super];
                                int end_row = supernode_ptr[super + 1];
                                //ncol is the number of columns in the off-diagonal block
                                int num_off_diag_col = Lp[start_row + 1] - Lp[start_row] - 1;
                                assert(num_off_diag_col >= 0);
                                int nrows = end_row - start_row;
                                assert(nrows > 0);
                                //Solving the independent part
                                //Copy x[Li[col]] into a continues buffer
                                for (int col_ptr = Lp[supernode_ptr[super]], k = 0;
                                col_ptr < Lp[supernode_ptr[super]] + num_off_diag_col; col_ptr++, k++) {
                                    tempvec[k] = x[Li[col_ptr]];
                                }
                                custom_blas::SpTrSv_MatVecCSR_BLAS(nrows, num_off_diag_col, &Lx[Lp[start_row]], tempvec.data(), &x[start_row]);
                                custom_blas::SpTrSv_LSolveCSR_BLAS(nrows, num_off_diag_col, &Lx[Lp[start_row] + num_off_diag_col], &x[start_row]);
                            }
                        }
                    }
                }
            }
            t1.measure_elapsed_time();
            sym_lib::copy_vector(0, n_, x_in_, x_);
            return t1;
        }

    public:
        /*
         * @brief Class constructor
         * @param L Sparse Matrix with CSR format
         * @param L_csc Sparse Matrix with CSC format
         * @param correct_x The right answer for x
         * @param name The name of the algorithm
         * @param nt number of threads
         */
        SpTrSv_LL_Blocked_Grouped_GLC(CSR *L, CSC *L_csc,
                                         double *correct_x, std::string name,
                                         int nt, bool sort=true)
                : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
            L1_csr_ = L;
            L1_csc_ = L_csc;
            correct_x_ = correct_x;
            nthreads = nt;
        };

        double getSchedulingTime() { return Scheduling_time.elapsed_time; }

        ~SpTrSv_LL_Blocked_Grouped_GLC() = default;
    };
//    //================================ SpTrSv DAGP LL ================================
//    class SpTrSv_LL_DAGP : public Sptrsv_LL_Serial {
//    protected:
//        int nthreads;
//        timing_measurement Scheduling_time;
//        int ngroups;
//        int *final_level_ptr, *final_part_ptr, *final_node_ptr;
//        int final_level_no;
//        double avg_par, max_diff, var;
//        int nlevels;
//        int num_parts;
//        void build_set() override {
//            Scheduling_time.start_timer();
//            //num_parts = num_threads_;
//            int n_graphs = 1;
//            int *parts;
//            auto *sptrsv_spmv_oneb = convert_to_one_based(L1_csc_);
//            CSC *p_graph = dagp_partition( sptrsv_spmv_oneb, n_, num_parts, parts);
//            std::vector<int> h_level_ptr, h_par_ptr, h_partition;
//            int n_levels;
//            get_h_level_set(n_, p_graph, num_parts, parts, final_level_no, final_level_ptr,
//                            final_part_ptr, final_node_ptr);
//            std::vector<int> level_ptr(final_level_ptr, final_level_ptr + final_level_no + 1);
//            std::vector<int> part_ptr(final_part_ptr, final_part_ptr + num_parts + 1);
//            std::vector<int> node_ptr(final_node_ptr, final_node_ptr + n_);
//
//            if(ANALYSE_STATISTIC){
//                std::cout << "UNPACKING IS DEACTIVATED" << std::endl;
//                computeStatistic();
//                //            std::vector<int> final_level_ptr_copy(final_level_ptr, final_level_ptr + final_level_no + 1);
//                //            std::vector<int> final_part_ptr_copy(final_part_ptr, final_part_ptr + final_level_ptr[final_level_no] + 1);
//                //            std::vector<int> final_node_ptr_copy(final_node_ptr, final_node_ptr + n_);
//                assert(test_unique(n_, final_node_ptr));
//            }
//
//
//            //print_hlevel_set("hlevel: \n", n_h_levels, h_level_ptr, h_part_ptr,
//            //  h_partition );
//
//            free(parts);
//            delete p_graph;
//            delete sptrsv_spmv_oneb;
//            Scheduling_time.measure_elapsed_time();
//        }
//
//        timing_measurement fused_code() override {
//            timing_measurement t1;
//            t1.start_timer();
//            sptrsv_csr_lbc(n_, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_,
//                           final_level_no, final_level_ptr,
//                           final_part_ptr, final_node_ptr);
//            t1.measure_elapsed_time();
//            sym_lib::copy_vector(0, n_, x_in_, x_);
//            return t1;
//        }
//
//    public:
//        /*
//         * @brief Class constructor
//         * @param L Sparse Matrix with CSR format
//         * @param L_csc Sparse Matrix with CSC format
//         * @param correct_x The right answer for x
//         * @param name The name of the algorithm
//         * @param nt number of threads
//         */
//        SpTrSv_LL_DAGP(CSR *L, CSC *L_csc,
//                       double *correct_x, std::string name, int num_parts,
//                           int nt)
//                           : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
//            L1_csr_ = L;
//            L1_csc_ = L_csc;
//            correct_x_ = correct_x;
//            nthreads = nt;
//            this->num_parts = num_parts;
//        };
//
//
//        void computeStatistic(){
//            double part_no = final_level_ptr[final_level_no];
//            ngroups = part_no;
//            avg_par = part_no/ final_level_no;
//            std::vector<double> max_part(final_level_no, 0);
//            std::vector<double> min_part(final_level_no, 0);
//            std::vector<double> level_weight(final_level_no, 0);
//            for(int lvl = 0; lvl < final_level_no; lvl++){
//                double total_nodes = 0;
//                double max_size = 0;
//                double min_size = n_;
//                for (int w_ptr = final_level_ptr[lvl]; w_ptr < final_level_ptr[lvl + 1]; ++w_ptr) {
//                    int cc_size = final_part_ptr[w_ptr + 1] - final_part_ptr[w_ptr];
//                    if(cc_size > max_size){
//                        max_size = cc_size;
//                    }
//                    if(cc_size < min_size){
//                        min_size = cc_size;
//                    }
//                    total_nodes += cc_size;
//                }
//                max_part[lvl] = max_size;
//                min_part[lvl] = min_size;
//                level_weight[lvl] = total_nodes;
//            }
//
//            max_diff = 0;
//            double total_weight = 0;
//            for(int l = 0; l < final_level_no; l++){
//                max_diff += level_weight[l] * (max_part[l] - min_part[l]);
//                total_weight += level_weight[l];
//            }
//            max_diff = max_diff / total_weight;
//            var = 0;
//        }
//
//        double getSchedulingTime() { return Scheduling_time.elapsed_time;}
//        void getStatistic(double& avg_par, double& max_diff, double& var){
//            avg_par = this->avg_par;
//            max_diff = this->max_diff;
//            var = this->var;
//        }
//        void getWaveStatistic(int& nlevels, int& part_no){
//            part_no = ngroups;
//            nlevels = this->final_level_no;
//        }
//        ~SpTrSv_LL_DAGP() {
//            delete[] final_level_ptr;
//            delete[] final_part_ptr;
//            delete[] final_node_ptr;
//        }
//    };
    //================================ SpTrSv SpMP LL ================================
    class SpTrsv_LL_SpMP : public Sptrsv_LL_Serial {
    protected:
        int nlevels, nthreads;
        std::vector<int> DAG_ptr;
        std::vector<int> DAG_set;
        timing_measurement scheduling_time;
        CSR* LBC_A;
        SpMP::LevelSchedule *p2pScheduleWithTransitiveReduction;
        const int* invPerm;
        void build_set() override {
            //Convert lbc CSR version to SpMP CSR version
            SpMP::CSR* A = new SpMP::CSR();
            Convert_LBCCSR_to_SpMP(LBC_A, A);
            scheduling_time.start_timer();
            /////////////////////////////////////////////////////////////////////////////
            // Construct schedules
            /////////////////////////////////////////////////////////////////////////////
            bool hasZeroDiag = A->hasZeroDiag();
            if (hasZeroDiag) {
                fprintf(stderr, "Matrix has a zero diagonal element. Can't do Gauss Seidel\n");
            }

            int *symRowPtr = NULL, *symColIdx = NULL, *symDiagPtr = NULL, *symExtPtr = NULL;
            bool wasSymmetric = getSymmetricNnzPattern(A, &symRowPtr, &symDiagPtr, &symExtPtr, &symColIdx);

            p2pScheduleWithTransitiveReduction = new SpMP::LevelSchedule;
            if (wasSymmetric) {
                p2pScheduleWithTransitiveReduction->constructTaskGraph(*A);
            }
            else {
                p2pScheduleWithTransitiveReduction->constructTaskGraph(
                        A->m, symRowPtr, symDiagPtr, symExtPtr, symColIdx,
                        SpMP::PrefixSumCostFunction(symRowPtr));

                FREE(symRowPtr);
                FREE(symColIdx);
                FREE(symDiagPtr);
                FREE(symExtPtr);
            }

            /////////////////////////////////////////////////////////////////////////////
            // Reorder matrix
            /////////////////////////////////////////////////////////////////////////////

            const int *perm = p2pScheduleWithTransitiveReduction->origToThreadContPerm;
            const int *invPerm = p2pScheduleWithTransitiveReduction->threadContToOrigPerm;
            assert(SpMP::isPerm(perm, A->m));
            assert(SpMP::isPerm(invPerm, A->m));
            this->invPerm = invPerm;

            delete A;
            scheduling_time.measure_elapsed_time();
        }

//        timing_measurement fused_code() override {
//            timing_measurement t1;
//            auto A = *L;
//            int base = A.getBase();
//            const int *rowptr = A.rowptr - base;
//            const int *colidx = A.colidx - base;
//            const double *values = A.values - base;
//            const double *idiag = A.idiag - base;
//            const SpMP::LevelSchedule &schedule = *p2pScheduleWithTransitiveReduction;
//            const int* perm = invPerm;
//            double* b = x_in_;
//            t1.start_timer();
//            #pragma omp parallel
//            {
//                int nthreads = omp_get_num_threads();
//                int tid = omp_get_thread_num();
//
//                const int ntasks = schedule.ntasks;
//                const short *nparents = schedule.nparentsForward;
//                const std::vector<int>& threadBoundaries = schedule.threadBoundaries;
//                const std::vector<int>& taskBoundaries = schedule.taskBoundaries;
//                int nPerThread = (ntasks + nthreads - 1)/nthreads;
//                int nBegin = std::min(nPerThread*tid, ntasks);
//                int nEnd = std::min(nBegin + nPerThread, ntasks);
//
//                volatile int *taskFinished = schedule.taskFinished;
//                int **parents = schedule.parentsForward;
//
//                memset((char *)(taskFinished + nBegin), 0, (nEnd - nBegin)*sizeof(int));
//
//                synk::Barrier::getInstance()->wait(tid);
//
//                for (int task = threadBoundaries[tid]; task < threadBoundaries[tid + 1]; ++task) {
//                    SPMP_LEVEL_SCHEDULE_WAIT;
//                    for (int i = taskBoundaries[task]; i < taskBoundaries[task + 1]; ++i) {
//                        int row = perm[i] + base;
//                        double sum = b[row];
//                        for (int j = rowptr[row]; j < rowptr[row + 1]; ++j) {
//                            sum -= values[j]*b[colidx[j]];
//                        }
//                        b[row] = sum*idiag[row];
//                    }
//                    SPMP_LEVEL_SCHEDULE_NOTIFY;
//                } // for each task
//            } // omp parallel
//
//
//            t1.measure_elapsed_time();
//            sym_lib::copy_vector(0, n_, x_in_, x_);
//            return t1;
//        }
        timing_measurement fused_code() override {
            timing_measurement t1;
            auto Lp = L1_csr_->p;
            auto Li = L1_csr_->i;
            auto Lx = L1_csr_->x;
            const SpMP::LevelSchedule &schedule = *p2pScheduleWithTransitiveReduction;
            const int* perm = invPerm;
            double* x = x_in_;
            t1.start_timer();
            #pragma omp parallel
            {
                int nthreads = omp_get_num_threads();
                int tid = omp_get_thread_num();

                const int ntasks = schedule.ntasks;
                const short *nparents = schedule.nparentsForward;
                const std::vector<int>& threadBoundaries = schedule.threadBoundaries;
                const std::vector<int>& taskBoundaries = schedule.taskBoundaries;
                int nPerThread = (ntasks + nthreads - 1)/nthreads;
                int nBegin = std::min(nPerThread*tid, ntasks);
                int nEnd = std::min(nBegin + nPerThread, ntasks);

                volatile int *taskFinished = schedule.taskFinished;
                int **parents = schedule.parentsForward;

                memset((char *)(taskFinished + nBegin), 0, (nEnd - nBegin)*sizeof(int));

                synk::Barrier::getInstance()->wait(tid);

                for (int task = threadBoundaries[tid]; task < threadBoundaries[tid + 1]; ++task) {
                    SPMP_LEVEL_SCHEDULE_WAIT;
                    for (int c = taskBoundaries[task]; c < taskBoundaries[task + 1]; ++c) {
                        int row = perm[c];
                        for (int j = Lp[row]; j < Lp[row + 1] - 1; j++)
                        {
                            x[row] -= Lx[j] * x[Li[j]];
                        }
                        x[row] /= Lx[Lp[row + 1] - 1];
                    }
                    SPMP_LEVEL_SCHEDULE_NOTIFY;
                } // for each task
            } // omp parallel
            t1.measure_elapsed_time();
            sym_lib::copy_vector(0, n_, x_in_, x_);
            return t1;
        }

    public:
        /*
         * @brief Class constructor
         * @param L Sparse Matrix with CSR format
         * @param L_csc Sparse Matrix with CSC format
         * @param correct_x The right answer for x
         * @param name The name of the algorithm
         * @param nt number of threads
         */
        SpTrsv_LL_SpMP(CSR* A, CSR *L, CSC *L_csc,
                       double *correct_x, std::string name,
                       int nt)
                : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
            this -> LBC_A = A;
            L1_csr_ = L;
            L1_csc_ = L_csc;
            correct_x_ = correct_x;
            this->nthreads = nt;
            omp_set_num_threads(nt);
        };

        double getNumberOfP2P(){
            auto nparents = p2pScheduleWithTransitiveReduction->nparentsForward;
            const std::vector<int>& threadBoundaries = p2pScheduleWithTransitiveReduction->threadBoundaries;
            double incomming_edges = 0;
            for(int tid = 0; tid < nthreads; tid++){
                for (int task = threadBoundaries[tid]; task < threadBoundaries[tid + 1]; ++task) {
                    incomming_edges += nparents[task];
                }
            }
            return incomming_edges;
        }

        double getSchedulingTime(){return scheduling_time.elapsed_time;};
        ~SpTrsv_LL_SpMP() {
            delete p2pScheduleWithTransitiveReduction;
        };
    };
    //================================ SpTrSv LBC LL ================================
    class SpTrSv_LL_LBC : public Sptrsv_LL_Serial {
    protected:
        int nthreads;
        timing_measurement Scheduling_time;
        std::vector<int> DAG_ptr;
        std::vector<int> DAG_set;
        int levelNo;
        int *final_level_ptr, *final_part_ptr, *final_node_ptr;
        int final_level_no;
        int part_no;
        int lp_, cp_, ic_;
        bool first_time;
        double avg_par, max_diff, var;
        std::vector<double> cost;
        void build_set() override {
            Scheduling_time.start_timer();
            cost.resize(n_, 0);
            auto CSC_Lp = L1_csc_->p;
            auto CSC_Li = L1_csc_->i;
            auto CSR_Lp = L1_csr_->p;
            auto CSR_Li = L1_csr_->i;
            GLC::costComputation(n_, CSC_Lp, CSC_Li, CSR_Lp, CSR_Li,
                                 GLC::SpTrSv_LL,
                                 NULL, NULL, false, cost);

            get_coarse_levelSet_DAG_CSC_tree(n_, L1_csc_->p, L1_csc_->i, -1,
                                             final_level_no,
                                             final_level_ptr,part_no,
                                             final_part_ptr,final_node_ptr,
                                             lp_, ic_, cp_, cost.data());
            std::cout << "The nlevels for LBC is: " << final_level_no << std::endl;
            if(ANALYSE_STATISTIC){
                std::cout << "Static is COMPUTING" << std::endl;
                computeStatistic();
            }

            Scheduling_time.measure_elapsed_time();
        }

        timing_measurement fused_code() override {
            timing_measurement t1;
            t1.start_timer();
            sptrsv_csr_lbc(n_, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_,
                           final_level_no, final_level_ptr,
                           final_part_ptr, final_node_ptr);
            t1.measure_elapsed_time();
            sym_lib::copy_vector(0, n_, x_in_, x_);
            return t1;
        }

    public:
        /*
         * @brief Class constructor
         * @param L Sparse Matrix with CSR format
         * @param L_csc Sparse Matrix with CSC format
         * @param correct_x The right answer for x
         * @param name The name of the algorithm
         * @param nt number of threads
         */
        SpTrSv_LL_LBC(CSR *L, CSC *L_csc, double *correct_x, std::string name, int lp)
                : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
            L1_csr_ = L;
            L1_csc_ = L_csc;
            correct_x_ = correct_x;
            lp_=lp;
            nthreads = lp;
        };

        void setP2_P3(int p2, int p3) {
//            this->cp_ = p3;
//            this->ic_ = p2;

            this->cp_ = p3;
            this->ic_ = p2;

        }


        void computeStatistic(){
            this->part_no = part_no;
            avg_par = part_no/ final_level_no;
            std::vector<double> max_part(final_level_no, 0);
            std::vector<double> min_part(final_level_no, 0);
            std::vector<double> level_weight(final_level_no, 0);
            for(int lvl = 0; lvl < final_level_no; lvl++){
                double total_nodes = 0;
                double max_size = 0;
                double min_size = n_;
                for (int w_ptr = final_level_ptr[lvl]; w_ptr < final_level_ptr[lvl + 1]; ++w_ptr) {
                    double cc_size = final_part_ptr[w_ptr + 1] - final_part_ptr[w_ptr];
                    double cc_cost = 0;
                    for(int i = final_part_ptr[w_ptr]; i < final_part_ptr[w_ptr + 1]; i++){
                        int node = final_node_ptr[i];
                        cc_cost += cost[node];
                    }
                    cc_size = cc_cost;
                    if(cc_size > max_size){
                        max_size = cc_size;
                    }
                    if(cc_size < min_size){
                        min_size = cc_size;
                    }
                    total_nodes += cc_size;
                }
                max_part[lvl] = max_size;
                min_part[lvl] = min_size;
                level_weight[lvl] = total_nodes;
            }

            max_diff = 0;
            double total_weight = 0;
            for(int l = 0; l < final_level_no; l++){
                max_diff += level_weight[l] * (max_part[l] - min_part[l]);
                total_weight += level_weight[l];
            }
            max_diff = max_diff / total_weight;
            var = 0;
        }
        void getWaveStatistic(int& nlevels, int& part_no){
            part_no = this->part_no;
            nlevels = this->final_level_no;
        }
        double getSchedulingTime() { return Scheduling_time.elapsed_time;}
        void getStatistic(double& avg_par, double& max_diff, double& var){
            avg_par = this->avg_par;
            max_diff = this->max_diff;
            var = this->var;
        }

        ~SpTrSv_LL_LBC() {
            delete final_level_ptr;
            delete final_node_ptr;
            delete final_part_ptr;
        };
    };
    //================================ SpTrSv MKL LL ================================
    class SpTrSv_LL_MKL : public FusionDemo {
    protected:
        int *diagonal_locs_;
        sparse_matrix_t csrA;
        CSR *A_oneb_;
        int ierr;
        double alpha = 1.0;
        const double beta = 0;
        struct matrix_descr descrB, descrL;
        timing_measurement Scheduling_time;

        void build_set() override {
            Scheduling_time.start_timer();
            MKL_INT job[7] = {1, 0, 0, 0, 0, 1, 1};
            MKL_INT info = 0;
            const MKL_INT size_csr = n_;
            char transa;
            char matdescra[6];
            transa = 'n';
            matdescra[0] = 't';
            matdescra[1] = 'l';
            matdescra[2] = 'n';
            matdescra[3] = 'c';

            // Create handle with matrix stored in CSR format
            auto  statl = mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ZERO,
                                                  size_csr,  // number of rows
                                                  size_csr,  // number of cols
                                                  L1_csr_->p,
                                                  L1_csr_->p+1,
                                                  L1_csr_->i,
                                                  L1_csr_->x);
            if (statl != SPARSE_STATUS_SUCCESS) {
                printf("analysis failed with %d;", statl);
            }

            // Create matrix descriptor
            descrL.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
            descrL.mode = SPARSE_FILL_MODE_LOWER;
            descrL.diag = SPARSE_DIAG_NON_UNIT;

            if (mkl_sparse_set_sv_hint(csrA, SPARSE_OPERATION_NON_TRANSPOSE,
                                       descrL, 10000) != SPARSE_STATUS_SUCCESS) {
                printf("optimization failed with %d;", statl);
            }

            if (mkl_sparse_optimize(csrA) != SPARSE_STATUS_SUCCESS) {
                printf("optimization failed with %d;", statl);
            }

            Scheduling_time.measure_elapsed_time();
        }

        void setting_up() override {
            std::fill_n(x_, n_, 1.0);
        }

        timing_measurement fused_code() override {
            timing_measurement t1;
            auto *y = new double[n_]();

            t1.start_timer();
            auto stat = mkl_sparse_d_trsv(SPARSE_OPERATION_NON_TRANSPOSE,
                                          alpha,
                                          csrA,
                                          descrL,
                                          x_,
                                          y);

            t1.measure_elapsed_time();
            if (stat != SPARSE_STATUS_SUCCESS) {
                printf("analysis failed with %d;", stat);
            }
            std::copy(y,y+n_,x_);
            delete []y;
            return t1;
        }

    public:
        SpTrSv_LL_MKL (CSR *A, CSC *A_csc,
                       double *correct_x, std::string name, int nt) :
                       FusionDemo(A->n, A->nnz, name) {
            L1_csr_ = A;
            L1_csc_ = A_csc;
            correct_x_ = correct_x;
            MKL_Set_Num_Threads(nt);
        };

        double getSchedulingTime() { return Scheduling_time.elapsed_time; }
        ~SpTrSv_LL_MKL () override {
        };
    };

}

#endif
