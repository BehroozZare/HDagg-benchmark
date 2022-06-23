//
// Created by kazem on 2020-04-23.
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


    //=============================== Right Looking ==============================
    bool spic0_csc_RL_serial(int n, int *Lp, int *Li, double *Lx)
    {
        for (int i = 0; i < n ; i++){
            Lx[Lp[i]] = Lx[Lp[i]] / sqrt(Lx[Lp[i]]);//S1
            for (int m = Lp[i] + 1; m < Lp[i + 1]; m++){
                Lx[m] = Lx[m] / Lx[Lp[i]];//S2
            }
            for (int m = Lp[i] + 1; m < Lp[i + 1]; m++) {
                for (int k = Lp[Li[m]] ; k < Lp[Li[m]+1]; k++){
                    for (int l = m; l < Lp[i + 1] ; l++){
                        if (Li[l] == Li[k] ){
                            Lx[k] -= Lx[m] * Lx[l]; //S3
                        }
                    }
                }
            }
        }
        return true;
    }

    void spico_csc_RL_levelset(int n, const int *Lp, const int *Li, double *Lx,
                            int levels, const int *levelPtr,
                            const int *levelSet)
    {
        #pragma omp parallel
        {
            for (int s = 0; s < levels; s++) {
                #pragma omp for schedule(auto)
                for (int k1 = levelPtr[s]; k1 < levelPtr[s + 1]; ++k1) {
                    int i = levelSet[k1];
                    Lx[Lp[i]] = Lx[Lp[i]] / sqrt(Lx[Lp[i]]);
                    for (int m = Lp[i] + 1; m < Lp[i + 1]; m++) {
                        Lx[m] = Lx[m] / Lx[Lp[i]];
                    }
                    for (int m = Lp[i] + 1; m < Lp[i + 1]; m++) {
                        for (int k = Lp[Li[m]]; k < Lp[Li[m] + 1]; k++) {
                            for (int l = m; l < Lp[i + 1]; l++) {
                                if (Li[l] == Li[k]) {
                                    #pragma omp atomic
                                    Lx[k] -= Lx[m] * Lx[l];
                                }
                            }
                        }
                    }
                }
            }
        };
    }

    void spico_csc_RL_group_levelset(int n, const int *Lp, const int *Li, double *Lx,
                                  int levels, const int *levelPtr,
                                  const int *levelSet, int *groupPtr, int *groupSet)
    {
        #pragma omp parallel
        {
            for (int l = 0; l < levels; ++l)
            {
                #pragma omp for schedule(auto)
                for (int li = levelPtr[l]; li < levelPtr[l + 1]; ++li)
                {
                    int lidx = levelSet[li];
                    for (int k1 = groupPtr[lidx]; k1 < groupPtr[lidx + 1]; ++k1)
                    {
                        int i = groupSet[k1];
                        double temp = Lx[Lp[i]];
                        Lx[Lp[i]] = Lx[Lp[i]] / sqrt(temp); //S1

                        for (int m = Lp[i] + 1; m < Lp[i + 1]; m++)
                        {
                            Lx[m] = Lx[m] / Lx[Lp[i]]; //S2
                        }

                        for (int m = Lp[i] + 1; m < Lp[i + 1]; m++)
                        {
                            for (int k = Lp[Li[m]]; k < Lp[Li[m] + 1]; k++)
                            {
                                for (int l = m; l < Lp[i + 1]; l++)
                                {
                                    if (Li[l] == Li[k])
                                    {
//                            if(rowIdx[l+1] <= rowIdx[k]){
#pragma omp atomic
                                        Lx[k] -= Lx[m] * Lx[l]; //S3
                                        //                            }
                                    }
                                    else if (Li[l] > Li[k])
                                        break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }


    void spico_csc_RL_lbc(int n, int *Lp, int *Li, double *Lx,
                       int level_no, int *level_ptr,
                       int *par_ptr, int *partition)
    {
        #pragma omp parallel
        {
            for (int i1 = 0; i1 < level_no; ++i1){
                #pragma omp for schedule(auto)
                for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1){
                    for (int k1 = par_ptr[j1]; k1 < par_ptr[j1 + 1]; ++k1){
                        int i = partition[k1];
                        Lx[Lp[i]] = Lx[Lp[i]] / sqrt(Lx[Lp[i]]);//S1
                        for (int m = Lp[i] + 1; m < Lp[i + 1]; m++) {
                            Lx[m] = Lx[m] / Lx[Lp[i]];//S2
                        }
                        for (int m = Lp[i] + 1; m < Lp[i + 1]; m++) {
                            for (int k = Lp[Li[m]]; k < Lp[Li[m] + 1]; k++) {
                                for (int l = m; l < Lp[i + 1]; l++) {
                                    if (Li[l] == Li[k]) {
                                        #pragma omp atomic
                                        Lx[k] -= Lx[m] * Lx[l]; //S3
                                    }
                                }
                            }
                        }
                    }
                }
            }
        };
    }


    void spico_csc_RL_group_lbc(int n, int *Lp, int *Li, double *Lx,
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
                            double temp = Lx[Lp[i]];
                            Lx[Lp[i]] = Lx[Lp[i]] / sqrt(temp); //S1

                            for (int m = Lp[i] + 1; m < Lp[i + 1]; m++)
                            {
                                Lx[m] = Lx[m] / Lx[Lp[i]]; //S2
                            }

                            for (int m = Lp[i] + 1; m < Lp[i + 1]; m++)
                            {
                                for (int k = Lp[Li[m]]; k < Lp[Li[m] + 1]; k++)
                                {
                                    for (int l = m; l < Lp[i + 1]; l++)
                                    {
                                        if (Li[l] == Li[k])
                                        {
//                            if(rowIdx[l+1] <= rowIdx[k]){
#pragma omp atomic
                                            Lx[k] -= Lx[m] * Lx[l]; //S3
                                            //                            }
                                        }
                                        else if (Li[l] > Li[k])
                                            break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        };
    }


    void spico_csc_group_levelset_v2(int n, const int *Lp, const int *Li, double *Lx,
                                  int levels, const int *levelPtr, const int *levelSet,
                                  int *groupPtr, int *groupSet)
    {
        #pragma omp parallel
        {
            for (int l = 0; l < levels; ++l) {
                #pragma omp for schedule(auto)
                for (int li = levelPtr[l]; li < levelPtr[l + 1]; ++li) {
                    int lidx = levelSet[li];
                    for (int k1 = groupPtr[lidx]; k1 < groupPtr[lidx + 1]; ++k1) {
                        int i = groupSet[k1];
                        Lx[Lp[i]] = Lx[Lp[i]] / sqrt(Lx[Lp[i]]);//S1
                        for (int m = Lp[i] + 1; m < Lp[i + 1]; m++) {
                            Lx[m] = Lx[m] / Lx[Lp[i]];//S2
                        }
                        for (int m = Lp[i] + 1; m < Lp[i + 1]; m++) {
                            for (int k = Lp[Li[m]]; k < Lp[Li[m] + 1]; k++) {
                                for (int l = m; l < Lp[i + 1]; l++) {
                                    if (Li[l] == Li[k]) {
                                        #pragma omp atomic
                                        Lx[k] -= Lx[m] * Lx[l]; //S3
                                    }
//                                    if(Li[l] > Li[k]){
//                                        break;
//                                    } else if (Li[l] == Li[k] ){
//                                        #pragma omp atomic
//                                        Lx[k] -= Lx[m] * Lx[l]; //S3
//                                        break;
//                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    };




    void spico_csc_lbc_reordered(int n, double *Lx, int *Lp, int *Li,
                       int level_no, int *level_ptr,
                       int *par_ptr, int *partition)
    {
        #pragma omp parallel
        {
            for (int i1 = 0; i1 < level_no; ++i1){
                #pragma omp for schedule(auto)
                for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1){
                    for (int k1 = par_ptr[j1]; k1 < par_ptr[j1 + 1]; ++k1){
                        int i = k1;
                        double temp = Lx[Lp[i]];
                        Lx[Lp[i]] = Lx[Lp[i]] / sqrt(temp); //S1

                        for (int m = Lp[i] + 1; m < Lp[i + 1]; m++){
                            Lx[m] = Lx[m] / Lx[Lp[i]]; //S2
                        }

                        for (int m = Lp[i] + 1; m < Lp[i + 1]; m++){
                            for (int k = Lp[Li[m]]; k < Lp[Li[m] + 1]; k++){
                                for (int l = m; l < Lp[i + 1]; l++){
                                    if (Li[l] == Li[k]){
                                        #pragma omp atomic
                                        Lx[k] -= Lx[m] * Lx[l]; //S3
                                    }
                                }
                            }
                        }
                    }
                }
            }
        };
    }
    //=============================== Left Looking ==============================
    bool spic0_csc_LL_serial(int n, const int *Lp, const int *Li, double *Lx,
                    const int* PrunePtr, const int* PruneSet)
    {
        bool isSPD = true;
        bool isPruneSetCorrect = true;
        for (int i = 0; i < n ; i++){
            //Gather value from other columns that affect column i
            //PrunePtr is row i in Left looking ICh0 -> iterate over columns that contribute to the current column i
            for (int m = PrunePtr[i]; m < PrunePtr[i + 1] - 1; m++) {
                int column_i_start = Lp[i];
                int column_j_start = 0;
                //Iterate over the rows to find a starting point
                for(int k = Lp[PruneSet[m]]; k < Lp[PruneSet[m] + 1]; k++){
                    if(Li[k] == i){
                        column_j_start = k;
                        isPruneSetCorrect = true;
                        break;
                    } else {
                        isPruneSetCorrect = false;
                    }
                }
                //Iterate over the column with information for current column i
                for (int l = column_i_start; l < Lp[i + 1] ; l++){
                    for (int k = column_j_start; k < Lp[PruneSet[m] + 1]; k++){
                        //Iterate over the column i values
                        if (Li[l] == Li[k]){
                            Lx[l] -= Lx[k] * Lx[column_j_start]; //S3
                            break;
                        }
                    }
                }
            }

            //Divide the column with the value of the diagonal element
            if(Lx[Lp[i]] < 0){
                isSPD = false;
            } else {
                Lx[Lp[i]] = Lx[Lp[i]]/sqrt(Lx[Lp[i]]);//S2
            }
            for (int m = Lp[i] + 1; m < Lp[i+1]; m++){
                Lx[m] = Lx[m] / Lx[Lp[i]];//S3
            }
        }
//        if(!isSPD){
//            std::cerr << "The matrix is not SPD" << std::endl;
//        }
//        if(!isPruneSetCorrect){
//            std::cerr << "Prune Set is not correct" << std::endl;
//        }
        return isSPD;
    }

    bool spic0_csc_LL_serial_test(int n, const int *Lp, const int *Li, double *Lx,
                             const int* PrunePtr, const int* PruneSet)
                             {
        bool isSPD = true;
        bool isPruneSetCorrect = true;
        for (int i = 0; i < n ; i++){
            //Gather value from other columns that affect column i
            //PrunePtr is row i in Left looking ICh0 -> iterate over columns that contribute to the current column i
            for (int m = PrunePtr[i]; m < PrunePtr[i + 1] - 1; m++) {
                const int j = PruneSet[m];
                int column_j_start;
                //Iterate over the rows to find the starting point
                for(int k = Lp[j]; k < Lp[j + 1]; k++){
                    if(Li[k] == i){
                        column_j_start = k;
                    }
                }
                int l1 = Lp[i];
                int l2 = column_j_start;
                while (l1 < Lp[i + 1] && l2 < Lp[j + 1]) {
                    if (Li[l1] == Li[l2]){     // matching column?
                        Lx[l1++] -= Lx[column_j_start] * Lx[l2++];
                    } else if (Li[l1] < Li[l2]){       // else proceed until we find matching columns
                        l1++;
                    } else {
                        l2++;
                    }
                }
            }

            //Divide the column with the value of the diagonal element
            if(Lx[Lp[i]] < 0){
                isSPD = false;
            } else {
                Lx[Lp[i]] = Lx[Lp[i]]/sqrt(Lx[Lp[i]]);//S2
            }
            for (int m = Lp[i] + 1; m < Lp[i+1]; m++){
                Lx[m] = Lx[m] / Lx[Lp[i]];//S3
            }
        }
        //        if(!isSPD){
        //            std::cerr << "The matrix is not SPD" << std::endl;
        //        }
        //        if(!isPruneSetCorrect){
        //            std::cerr << "Prune Set is not correct" << std::endl;
        //        }
        return isSPD;
     }


    bool spic0_csc_LL_Levelset(int n, const int *Lp, const int *Li, double *Lx,
                                 const int* PrunePtr, const int* PruneSet,
                                 int nLevels, const int* LevelPtr, const int* LevelSet) {
        bool isSPD = true;
        bool isPruneSetCorrect = true;
        #pragma omp parallel shared(isSPD, isPruneSetCorrect)
        {
            for (int i1 = 0; i1 < nLevels; ++i1) {
                #pragma omp  for schedule(static)
                for (int j1 = LevelPtr[i1]; j1 < LevelPtr[i1 + 1]; ++j1) {
                    int i = LevelSet[j1];
                    //Gather value from other columns that affect column i
                    //PrunePtr is row i in Left looking ICh0 -> iterate over columns that contribute to the current column i
                    for (int m = PrunePtr[i]; m < PrunePtr[i + 1] - 1; m++) {
                        int column_i_start = Lp[i];
                        int column_j_start = 0;
                        //Iterate over the column with information for current column i
                        for (int k = Lp[PruneSet[m]]; k < Lp[PruneSet[m] + 1]; k++) {
                            if (Li[k] == i) {
                                column_j_start = k;
                                isPruneSetCorrect = true;
                                break;
                            } else {
                                isPruneSetCorrect = false;
                            }
                        }

                        for (int l = column_i_start; l < Lp[i + 1]; l++) {
                            for (int k = column_j_start; k < Lp[PruneSet[m] + 1]; k++) {
                                //Iterate over the column i values
                                if (Li[l] == Li[k]) {
                                    Lx[l] -= Lx[k] * Lx[column_j_start]; //S3
                                    break;
                                }
                            }
                        }
                    }

                    //Divide the column with the value of the diagonal element
                    if (Lx[Lp[i]] < 0) {
                        isSPD = false;
                    } else {
                        Lx[Lp[i]] = Lx[Lp[i]] / sqrt(Lx[Lp[i]]);//S2
                    }
                    for (int m = Lp[i] + 1; m < Lp[i + 1]; m++) {
                        Lx[m] = Lx[m] / Lx[Lp[i]];//S3
                    }
                }
            }
        }
//        if(!isSPD){
//            std::cerr << "The matrix is not SPD" << std::endl;
//        }
//        if(!isPruneSetCorrect){
//            std::cerr << "Prune Set is not correct" << std::endl;
//        }
        return isSPD;
    }

    bool spic0_csc_LL_group_Levelset(int n, const int *Lp, const int *Li, double *Lx,
                                     const int* PrunePtr, const int* PruneSet,
                                     int levels, const int *levelPtr, const int *levelSet,
                                     int *groupPtr, int *groupSet)
    {
        bool isSPD = true;
        bool isPruneSetCorrect = true;
        #pragma omp parallel shared(isSPD, isPruneSetCorrect)
        {
            for (int l = 0; l < levels; ++l)
            {
                #pragma omp for schedule(auto)
                for (int li = levelPtr[l]; li < levelPtr[l + 1]; ++li)
                {
                    int lidx = levelSet[li];
                    for (int k1 = groupPtr[lidx]; k1 < groupPtr[lidx + 1]; ++k1)
                    {
                        int i = groupSet[k1];
                        //Gather value from other columns that affect column i
                        //PrunePtr is row i in Left looking ICh0 -> iterate over columns that contribute to the current column i
                        for (int m = PrunePtr[i]; m < PrunePtr[i + 1] - 1; m++) {
                            int column_i_start = Lp[i];
                            int column_j_start = 0;
                            //Iterate over the column with information for current column i
                            for (int k = Lp[PruneSet[m]]; k < Lp[PruneSet[m] + 1]; k++) {
                                if (Li[k] == i) {
                                    column_j_start = k;
                                    isPruneSetCorrect = true;
                                    break;
                                } else {
                                    isPruneSetCorrect = false;
                                }
                            }

                            for (int l = column_i_start; l < Lp[i + 1]; l++) {
                                for (int k = column_j_start; k < Lp[PruneSet[m] + 1]; k++) {
                                    //Iterate over the column i values
                                    if (Li[l] == Li[k]) {
                                        Lx[l] -= Lx[k] * Lx[column_j_start]; //S3
                                        break;
                                    }
                                }
                            }
                        }

                        //Divide the column with the value of the diagonal element
                        if (Lx[Lp[i]] < 0) {
                            isSPD = false;
                        } else {
                            Lx[Lp[i]] = Lx[Lp[i]] / sqrt(Lx[Lp[i]]);//S2
                        }
                        for (int m = Lp[i] + 1; m < Lp[i + 1]; m++) {
                            Lx[m] = Lx[m] / Lx[Lp[i]];//S3
                        }
                    }
                }
            }
        }
        return isSPD;
    }


    bool spic0_csc_LL_LBC(int n, const int *Lp, const int *Li, double *Lx,
                               const int* PrunePtr, const int* PruneSet,
                               int nLevels, const int* LevelPtr, const int* ParPtr, const int* Partitions) {
        bool isSPD = true;
        bool isPruneSetCorrect = true;
        #pragma omp parallel shared(isSPD, isPruneSetCorrect)
        {
            for (int i1 = 0; i1 < nLevels; ++i1){
                #pragma omp for schedule(auto)
                for (int j1 = LevelPtr[i1]; j1 < LevelPtr[i1 + 1]; ++j1){
                    for (int k1 = ParPtr[j1]; k1 < ParPtr[j1 + 1]; ++k1){
                        int i = Partitions[k1];
                        //Gather value from other columns that affect column i
                        //PrunePtr is row i in Left looking ICh0 -> iterate over columns that contribute to the current column i
                        for (int m = PrunePtr[i]; m < PrunePtr[i + 1] - 1; m++){
                            int column_i_start = Lp[i];
                            int column_j_start = 0;
                            //Iterate over the column with information for current column i
                            for (int k = Lp[PruneSet[m]]; k < Lp[PruneSet[m] + 1]; k++) {
                                if (Li[k] == i) {
                                    column_j_start = k;
                                    isPruneSetCorrect = true;
                                    break;
                                } else {
                                    isPruneSetCorrect = false;
                                }
                            }

                            for (int l = column_i_start; l < Lp[i + 1]; l++) {
                                for (int k = column_j_start; k < Lp[PruneSet[m] + 1]; k++) {
                                    //Iterate over the column i values
                                    if (Li[l] == Li[k]) {
                                        Lx[l] -= Lx[k] * Lx[column_j_start]; //S3
                                        break;
                                    }
                                }
                            }
                        }

                        //Divide the column with the value of the diagonal element
                        if (Lx[Lp[i]] < 0) {
                            isSPD = false;
                        } else {
                            Lx[Lp[i]] = Lx[Lp[i]] / sqrt(Lx[Lp[i]]);//S2
                        }
                        for (int m = Lp[i] + 1; m < Lp[i + 1]; m++) {
                            Lx[m] = Lx[m] / Lx[Lp[i]];//S3
                        }
                    }
                }
            }
        }
//        if(!isSPD){
//            std::cerr << "The matrix is not SPD" << std::endl;
//        }
//        if(!isPruneSetCorrect){
//            std::cerr << "Prune Set is not correct" << std::endl;
//        }
        return isSPD;
    }

    bool spic0_csc_LL_group_LBC(int n, const int *Lp, const int *Li, double *Lx,
                                     const int* PrunePtr, const int* PruneSet,
                                    int nLevels, const int* LevelPtr, const int* ParPtr, const int* Partitions,
                                     int *groupPtr, int *groupSet)
    {
        bool isSPD = true;
        bool isPruneSetCorrect = true;
        #pragma omp parallel shared(isSPD, isPruneSetCorrect)
        {
            for (int i1 = 0; i1 < nLevels; ++i1){
                #pragma omp for schedule(auto)
                for (int j1 = LevelPtr[i1]; j1 < LevelPtr[i1 + 1]; ++j1){
                    for (int k1 = ParPtr[j1]; k1 < ParPtr[j1 + 1]; ++k1){
                        int lidx = Partitions[k1];
                        for (int k11 = groupPtr[lidx]; k11 < groupPtr[lidx + 1]; ++k11) {
                            int i = groupSet[k11];
                            //Gather value from other columns that affect column i
                            //PrunePtr is row i in Left looking ICh0 -> iterate over columns that contribute to the current column i
                            for (int m = PrunePtr[i]; m < PrunePtr[i + 1] - 1; m++) {
                                int column_i_start = Lp[i];
                                int column_j_start = 0;
                                //Iterate over the column with information for current column i
                                for (int k = Lp[PruneSet[m]]; k < Lp[PruneSet[m] + 1]; k++) {
                                    if (Li[k] == i) {
                                        column_j_start = k;
                                        isPruneSetCorrect = true;
                                        break;
                                    } else {
                                        isPruneSetCorrect = false;
                                    }
                                }

                                for (int l = column_i_start; l < Lp[i + 1]; l++) {
                                    for (int k = column_j_start; k < Lp[PruneSet[m] + 1]; k++) {
                                        //Iterate over the column i values
                                        if (Li[l] == Li[k]) {
                                            Lx[l] -= Lx[k] * Lx[column_j_start]; //S3
                                            break;
                                        }
                                    }
                                }
                            }

                            //Divide the column with the value of the diagonal element
                            if (Lx[Lp[i]] < 0) {
                                isSPD = false;
                            } else {
                                Lx[Lp[i]] = Lx[Lp[i]] / sqrt(Lx[Lp[i]]);//S2
                            }
                            for (int m = Lp[i] + 1; m < Lp[i + 1]; m++) {
                                Lx[m] = Lx[m] / Lx[Lp[i]];//S3
                            }
                        }
                    }
                }
            }
        }
        return isSPD;
    }

    //=============================== Up looking Looking ==============================
    bool spic0_csr_uL_serial(int n, const int *Lp, const int *Li, double *Lx, double *new_data)
    {
        bool finish_flag = true;
        for (int i = 0; i < n; ++i) {
            for (int k = Lp[i]; k < Lp[i+1]; ++k) {
                const int j = Li[k];     // column
                double dp = sparse_dot_product(
                        Lp[i], Lp[i+1] - 1,       // i-th row minus diagonal
                        Lp[j], Lp[j+1] - 1,       // j-th row minus diagonal
                        Li, &new_data[0]);

                const double A_ij = Lx[k];

                if (j < i) {        // below diagonal?
                    const double L_jj = new_data[Lp[j+1] - 1];    // diagonal is last entry of j-th row
                    if(L_jj != 0){
                        new_data[k] = (A_ij - dp) / L_jj;
                    }
                } else if (j == i){    // on the diagonal?
                    if(A_ij - dp >= 0){
                        new_data[k] = std::sqrt(A_ij - dp);
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
        return finish_flag;
    }

    bool spic0_csr_UL_Levelset(int n, const int *Lp, const int *Li, double *Lx, double* new_data,
                               int nLevels, const int* LevelPtr, const int* LevelSet) {
        bool finish_flag = true;
        #pragma omp parallel
        {
            for (int i1 = 0; i1 < nLevels; ++i1) {
                #pragma omp  for schedule(static)
                for (int j1 = LevelPtr[i1]; j1 < LevelPtr[i1 + 1]; ++j1) {
                    int i = LevelSet[j1];
                    for (int k = Lp[i]; k < Lp[i+1]; ++k) {
                        const int j = Li[k];     // column
                        double dp = sparse_dot_product(
                                Lp[i], Lp[i+1] - 1,       // i-th row minus diagonal
                                Lp[j], Lp[j+1] - 1,       // j-th row minus diagonal
                                Li, &new_data[0]);

                        const double A_ij = Lx[k];

                        if (j < i) {        // below diagonal?
                            const double L_jj = new_data[Lp[j+1] - 1];    // diagonal is last entry of j-th row
                            if(L_jj != 0){
                                new_data[k] = (A_ij - dp) / L_jj;
                            }
                        } else if (j == i){    // on the diagonal?
                            if(A_ij - dp >= 0){
                                new_data[k] = std::sqrt(A_ij - dp);
                            }
                        } else                // above diagonal -- input should be triangular!
                            throw std::logic_error("Matrix passed in the wrong format - should be triangular");
                    }
                }
            }
        }

//        #ifndef NDEBUG
//                for (int l = 0; l < n; ++l) {
//            assert(Lx[Lp[l]-1] != 0);
//            if(Lx[Lp[l]-1] == 0){
//                finish_flag = false;
//                break;
//            }
//        }
//        #endif
        return finish_flag;
    }

    bool spic0_csr_UL_group_Levelset(int n, const int *Lp, const int *Li, double *Lx, double* new_data,
                               int nLevels, const int* LevelPtr, const int* LevelSet,
                               const int* group_ptr, const int* group_set) {
        bool finish_flag = true;
        #pragma omp parallel
        {
            for (int i1 = 0; i1 < nLevels; ++i1) {
                #pragma omp  for schedule(static)
                for (int j1 = LevelPtr[i1]; j1 < LevelPtr[i1 + 1]; ++j1) {
                    int g_ptr = LevelSet[j1];
                    for(int g = group_ptr[g_ptr]; g < group_ptr[g_ptr + 1]; g++){
                        int i = group_set[g];
                        for (int k = Lp[i]; k < Lp[i+1]; ++k) {
                            const int j = Li[k];     // column
                            double dp = sparse_dot_product(
                                    Lp[i], Lp[i+1] - 1,       // i-th row minus diagonal
                                    Lp[j], Lp[j+1] - 1,       // j-th row minus diagonal
                                    Li, &new_data[0]);

                            const double A_ij = Lx[k];

                            if (j < i) {        // below diagonal?
                                const double L_jj = new_data[Lp[j+1] - 1];    // diagonal is last entry of j-th row
                                if(L_jj != 0){
                                    new_data[k] = (A_ij - dp) / L_jj;
                                }
                            } else if (j == i){    // on the diagonal?
                                if(A_ij - dp >= 0){
                                    new_data[k] = std::sqrt(A_ij - dp);
                                }
                            } else                // above diagonal -- input should be triangular!
                            throw std::logic_error("Matrix passed in the wrong format - should be triangular");
                        }
                    }
                }
            }
        }

        //        #ifndef NDEBUG
        //                for (int l = 0; l < n; ++l) {
        //            assert(Lx[Lp[l]-1] != 0);
        //            if(Lx[Lp[l]-1] == 0){
        //                finish_flag = false;
        //                break;
        //            }
        //        }
        //        #endif
        return finish_flag;
    }

    bool spic0_csr_UL_LBC(int n, const int *Lp, const int *Li, double *Lx, double* new_data,
                          int nLevels, const int* LevelPtr, const int* ParPtr, const int* Partitions) {
        bool finish_flag = true;
        #pragma omp parallel
        {
            for (int i1 = 0; i1 < nLevels; ++i1){
                #pragma omp for schedule(auto)
                for (int j1 = LevelPtr[i1]; j1 < LevelPtr[i1 + 1]; ++j1){
                    for (int k1 = ParPtr[j1]; k1 < ParPtr[j1 + 1]; ++k1){
                        int i = Partitions[k1];
                        for (int k = Lp[i]; k < Lp[i+1]; ++k) {
                            const int j = Li[k];     // column
                            double dp = sparse_dot_product(
                                    Lp[i], Lp[i+1] - 1,       // i-th row minus diagonal
                                    Lp[j], Lp[j+1] - 1,       // j-th row minus diagonal
                                    Li, &new_data[0]);

                            const double A_ij = Lx[k];

                            if (j < i) {        // below diagonal?
                                const double L_jj = new_data[Lp[j+1] - 1];    // diagonal is last entry of j-th row
                                if(L_jj != 0){
                                    new_data[k] = (A_ij - dp) / L_jj;
                                }
                            } else if (j == i){    // on the diagonal?
                                if(A_ij - dp >= 0){
                                    new_data[k] = std::sqrt(A_ij - dp);
                                }
                            } else                // above diagonal -- input should be triangular!
                                throw std::logic_error("Matrix passed in the wrong format - should be triangular");
                        }
                    }
                }
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
        return finish_flag;
    }

} // namespace sym_lib