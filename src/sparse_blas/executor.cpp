//
// Created by labuser (Bangtian Liu) on 9/16/20.
//

#include <executor.h>
#include <vector>


namespace sym_lib
{
    int fs_csr_executor_sgroup(int n, int *Lp, int *Li, double *Lx, double *b, double *x, int *groupPtr, int *groupSet,  int ngroup,
                               int levels, int *levelPtr, int *levelSet) {
     if (!Lp || !Li || !x) return (0);
#pragma omp parallel
     {
      for (int l = 0; l < levels; ++l) {
#pragma omp for schedule(auto)
       for (int li = levelPtr[l]; li < levelPtr[l + 1]; ++li) {
        int lidx = levelSet[li];
        for (int k1 = groupPtr[lidx]; k1 < groupPtr[lidx + 1]; ++k1) {
         int j1 = groupSet[k1];
//                        double tmp = b[j1];
         for (int j = Lp[j1]; j < Lp[j1 + 1] - 1; ++j) {
          x[j1] -= Lx[j] * x[Li[j]];
         }
         x[j1] /= Lx[Lp[j1 + 1] - 1];
        }
       }
      }
     };
        return 0;
    }

    int CSRGeneralGroupingExecuter(int n, int *Lp, int *Li, double *Lx, double *x,
                                   int* group_range, int levels, int *levelPtr, int *levelSet) {
        if (!Lp || !Li || !x) return (0);
        #pragma omp parallel default(none) shared(n, Lp, Li, Lx, x, group_range, \
        levels, levelPtr, levelSet)
        {
            for(int lvl = 0; lvl < levels; lvl++) {
                #pragma omp for schedule(auto)
                for(int li = levelPtr[lvl]; li < levelPtr[lvl + 1]; li++) {
                    int node = levelSet[li];
                    for (int k1 = group_range[node]; k1 < group_range[node + 1]; k1++) {
                        for (int j = Lp[k1]; j < Lp[k1 + 1] - 1; j++) {
                            x[k1] -= Lx[j] * x[Li[j]];
                        }
                        x[k1] /= Lx[Lp[k1 + 1] - 1];
                    }
                }
            }
        }
        return 0;
    }

    int CSRGeneralGroupingExecutorWRegion(const int *Lp, const int *Li, const double *Lx,
                                          double *x,
                                          const int* group_range,
                                          int levels,
                                          const int *levelPtr, const int *levelSet) {
        if (!Lp || !Li || !x) return (0);
        #pragma omp parallel default(none) shared(Lp, Li, Lx, x, group_range, \
        levels, levelPtr, levelSet)
        {
            for(int lvl = 0; lvl < levels - 1; lvl++) {
                #pragma omp for schedule(auto)
                for(int li = levelPtr[lvl]; li < levelPtr[lvl + 1]; li++) {
                    int node = levelSet[li];
                    for (int k1 = group_range[node]; k1 < group_range[node + 1]; k1++) {
                        for (int j = Lp[k1]; j < Lp[k1 + 1] - 1; j++) {
                            x[k1] -= Lx[j] * x[Li[j]];
                        }
                        x[k1] /= Lx[Lp[k1 + 1] - 1];
                    }
                }
            }
        }
        //Serial region
        for(int li = levelPtr[levels - 1]; li < levelPtr[levels]; li++) {
            int node = levelSet[li];
            for (int k1 = group_range[node]; k1 < group_range[node + 1]; k1++) {
                for (int j = Lp[k1]; j < Lp[k1 + 1] - 1; j++) {
                    x[k1] -= Lx[j] * x[Li[j]];
                }
                x[k1] /= Lx[Lp[k1 + 1] - 1];
            }
        }
        return 0;
    }


    int CSRLevelExecuter(int n, int *Lp, int *Li, double *Lx, double *x,
                                   int levels, int *levelPtr, int *levelSet) {
        if (!Lp || !Li || !x) return (0);
        #pragma omp parallel default(none) shared(n, Lp, Li, Lx, x,\
        levels, levelPtr, levelSet)
        {
            for(int lvl = 0; lvl < levels; lvl++) {
                #pragma omp for schedule(auto)
                for(int li = levelPtr[lvl]; li < levelPtr[lvl + 1]; li++) {
                    int k1 = levelSet[li];
                    for (int j = Lp[k1]; j < Lp[k1 + 1] - 1; j++) {
                        x[k1] -= Lx[j] * x[Li[j]];
                    }
                    x[k1] /= Lx[Lp[k1 + 1] - 1];
                }
            }
        }
        return 0;
    }

    void fs_csr_stat(int n, int *Lp, int *Li,  int &flops, int &access_nnz, int &reuse_nnz)
    {
//        int flops=0;
     int t_flops=0;
     int t_mem=0;
     int t_reuse=0;

     for (int i = 0; i < n; i++) {
      for (int j = Lp[i]; j < Lp[i + 1] - 1; j++) {
//                x[i] -= Lx[j] * x[Li[j]];
       t_flops+=2;
       t_mem+=3;
      }
//            x[i] /= Lx[Lp[i + 1] - 1];
      t_mem+=2;
      t_flops+=1;
     }


     for (int i = 1; i < n; ++i) {
      for (int j = Lp[i]; j < Lp[i+1] -1; ++j) {
       if(Li[j]==i-1){
        ++t_reuse;
       }
      }
     }


     flops=t_flops;
     access_nnz=t_mem;
     reuse_nnz=t_reuse;
    }


    void ic0_csc_stat(int n,  int * colPtr, int *rowIdx, int &flops, int &access_nnz, int &reuse_nnz)
    {
     int i, k,l,m;
     double temp;

     int t_flops=0;
     int t_mem=0;
     int t_reuse=0;


     for (i = 0; i < n ; i++){
      t_flops += colPtr[i+1]-colPtr[i]+1;
      t_mem +=(colPtr[i+1]-colPtr[i]);
      for (m = colPtr[i] + 1; m < colPtr[i+1]; m++) {
       for (k = colPtr[rowIdx[m]] ; k < colPtr[rowIdx[m]+1]; k++){
        for ( l = m; l < colPtr[i+1] ; l++){
         ++t_mem;
         if (rowIdx[l] == rowIdx[k] ){
//                            if(rowIdx[l+1] <= rowIdx[k]){
//                            val[k] -= val[m]* val[l]; //S3
          t_flops+=2;
          t_mem+=3;
//                            }

         }
        }
       }
      }
     }

     flops=t_flops;
     access_nnz=t_mem;
     reuse_nnz=t_reuse;
    }





    void fs_csr_levelset_stat(int *Lp, int *Li, int *groupPtr, int *groupSet,
                              int levels, int *levelPtr, int *levelSet, int *lcost)
    {

     int *cost = (int *)malloc(sizeof(int)*omp_get_max_threads());
     memset(cost, 0, sizeof(int)*omp_get_max_threads());

     int len = omp_get_max_threads();

     if (!Lp || !Li ) return ;
//#pragma omp parallel
     {
      for (int l = 0; l < levels; ++l) {
#pragma omp parallel for schedule(auto)
       for (int li = levelPtr[l]; li < levelPtr[l + 1]; ++li) {
        int tidx = omp_get_thread_num();
        int t_cost=0;
        int lidx = levelSet[li];

        for (int k1 = groupPtr[lidx]; k1 < groupPtr[lidx + 1]; ++k1) {
         int j1 = groupSet[k1];
//                        double tmp = b[j1];
         cost[tidx]+= (Lp[j1+1]-Lp[j1]);
//                        for (int j = Lp[j1]; j < Lp[j1 + 1] - 1; ++j) {
//                            tmp -= Lx[j] * x[Li[j]];
//                        }
//                        x[j1] = tmp / Lx[Lp[j1 + 1] - 1];
        }
       }


       std::vector<int> temp_cost;
       for (int i = 0; i < len; ++i) {
        for (int j = i+1; j < len; ++j) {
//                        if(i!=j){
         temp_cost.push_back(std::abs(cost[i]-cost[j]));
//                        }
        }
       }
       std::sort(temp_cost.begin(), temp_cost.end(), std::greater<int>());
       if(temp_cost.size()==0)
        lcost[l]=0;
       else
        lcost[l]=temp_cost[0];


       temp_cost.clear();
       memset(cost, 0, len*sizeof(int));
      }
     };
    }

    void sptrsv_csr_lbc_stat(int n, int *Lp, int *Li,
                             int level_no, int *level_ptr,
                             int *par_ptr, int *partition, int *lcost) {
     int *cost = (int *)malloc(sizeof(int)*omp_get_max_threads());
     memset(cost, 0, sizeof(int)*omp_get_max_threads());

     int len = omp_get_max_threads();

     for (int i1 = 0; i1 < level_no; ++i1) {
      {
#pragma omp  parallel for schedule(auto)
       for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1) {
        int tidx = omp_get_thread_num();

        for (int k1 = par_ptr[j1]; k1 < par_ptr[j1 + 1]; ++k1) {
         int i = partition[k1];
         cost[tidx] += (Lp[i+1] - Lp[i]);
//                        for (int j = Lp[i]; j < Lp[i + 1] - 1; j++) {
//                            x[i] -= Lx[j] * x[Li[j]];
//                        }
//                        x[i] /= Lx[Lp[i + 1] - 1];
        }
       }
      }


      std::vector<int> temp_cost;
      for (int i = 0; i < len; ++i) {
       for (int j = i+1; j < len; ++j) {
        temp_cost.push_back(std::abs(cost[i]-cost[j]));
       }
      }

      std::sort(temp_cost.begin(), temp_cost.end(), std::greater<int>());
      if(temp_cost.size()==0)
       lcost[i1]=0;
      else
       lcost[i1]=temp_cost[0];
      temp_cost.clear();
      memset(cost, 0, len*sizeof(int));
     }

     free(cost);
    }



    void sptrsv_csr_group_lbc_stat(int n, int *Lp, int *Li,
                                   int level_no, int *level_ptr,
                                   int *par_ptr, int *partition, int *groupPtr, int *groupSet, int *lcost) {
     int *cost = (int *)malloc(sizeof(int)*omp_get_max_threads());
     memset(cost, 0, sizeof(int)*omp_get_max_threads());

     int len = omp_get_max_threads();

     for (int i1 = 0; i1 < level_no; ++i1) {
      {
#pragma omp  parallel for schedule(auto)
       for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1) {
        int tidx = omp_get_thread_num();

        for (int k1 = par_ptr[j1]; k1 < par_ptr[j1 + 1]; ++k1) {
         int p = partition[k1];
         for (int k = groupPtr[p]; k < groupPtr[p+1]; ++k) {
          int i = groupSet[k];
          cost[tidx] += (Lp[i+1] - Lp[i]);
         }
        }
       }
      }


      std::vector<int> temp_cost;
      for (int i = 0; i < len; ++i) {
       for (int j = i+1; j < len; ++j) {
        temp_cost.push_back(std::abs(cost[i]-cost[j]));
       }
      }

      std::sort(temp_cost.begin(), temp_cost.end(), std::greater<int>());
      if(temp_cost.size()==0)
       lcost[i1]=0;
      else
       lcost[i1]=temp_cost[0];
      temp_cost.clear();
      memset(cost, 0, len*sizeof(int));
     }

     free(cost);
    }






}
