//
// Created by Kazem on 10/14/19.
//

#ifndef PROJECT_LBC_H
#define PROJECT_LBC_H

#include <cstddef>

namespace sym_lib{

 int get_coarse_levelSet_DAG_CSC(size_t n,
                                 int *lC,
                                 int *lR,
                                 int &finaLevelNo,
                                 int *&finaLevelPtr,
                                 int &partNo,
                                 int *&finalPartPtr,
                                 int *&finalNodePtr,
                                 int innerParts,
                                 int minLevelDist,
                                 int divRate,
                                 double *nodeCost);

 int get_coarse_levelSet_DAG_CSC(size_t n,
                                 int *lC,
                                 int *lR,
                                 int &finaLevelNo,
                                 int *&finaLevelPtr,
                                 int &partNo,
                                 int *&finalPartPtr,
                                 int *&finalNodePtr,
                                 int innerParts,
                                 int minLevelDist,
                                 int divRate,
                                 double *nodeCost,
                                 bool bin_pack);
 /// Comnputes coarsened level-set for CSC DAG using etree
 /// \param n
 /// \param lC
 /// \param lR
 /// \param finaLevelNo
 /// \param finaLevelPtr
 /// \param partNo
 /// \param finalPartPtr
 /// \param finalNodePtr
 /// \param innerParts
 /// \param minLevelDist
 /// \param divRate
 /// \param nodeCost
 /// \return
 int get_coarse_levelSet_DAG_CSC_tree(size_t n,
                                      int *lC,
                                      int *lR,
                                      int stype,
                                      int &finaLevelNo,
                                      int *&finaLevelPtr,
                                      int &partNo,
                                      int *&finalPartPtr,
                                      int *&finalNodePtr,
                                      int innerParts,
                                      int minLevelDist,
                                      int divRate,
                                      double *nodeCost);

 /// Computes coarsened level set by working on the DAG directly.
 int get_coarse_Level_set_DAG_CSC03(size_t n,
                                    int *lC,
                                    int *lR,
                                    int &finaLevelNo,
                                    int *&finaLevelPtr,
                                    int &partNo,
                                    int *&finalPartPtr,
                                    int *&finalNodePtr,
                                    int innerParts,
                                    int minLevelDist,
                                    int divRate,
                                    double *nodeCost);

 /*
  * @brief This function returns the merged levels schedule
  * @param n Number of nodes in the DAG
  * @param lC is the columns pointer in CSC data structure
  * @param lR is the row index in CSC format
  * @return finalLevelNo is the number of levels in the merged DAG
  * @return finalLevelPtr is the pointer array that show the start and end point of w-partitions in partNo array
  * @return partNo is the number of w-partitions
  * @return finalPartPtr is the w-partitions sorted by their levels
  * @return finalNodePtr is the nodes sorted by their w-partitions
  * @param innerParts
  * @param minLevelDist
  * @param divRate
  * @param NumThreads is the number of processors
  * @param nodeCost needed for bin-packing. It is a cost that assigned to each node based on a cost function
  * @param binPacking is a flag that determines whether we want to use bin-packing or not
  */
 int get_coarse_Level_set_DAG_CSC03_parallel(
  size_t n, int *lC, int *lR, int &finaLevelNo, int *&finaLevelPtr, int &partNo,
  int *&finalPartPtr, int *&finalNodePtr, int innerParts, int minLevelDist,
  int divRate, int numThreads, double *nodeCost, bool binPacking = false);


    /*
    * @brief This function returns the merged levels schedule. Note that it needs the levelset information
    * thus it can reuse one levelset computation in multiple merging schedule
    * @param n Number of nodes in the DAG
    * @param lC is the columns pointer in CSC data structure
    * @param lR is the row index in CSC format
    * @return finalLevelNo is the number of levels in the merged DAG
    * @return finalLevelPtr is the pointer array that show the start and end point of w-partitions in partNo array
    * @return partNo is the number of w-partitions
    * @return finalPartPtr is the w-partitions sorted by their levels
    * @return finalNodePtr is the nodes sorted by their w-partitions
    * @param innerParts
    * @param minLevelDist
    * @param divRate
    * @param NumThreads is the number of processors
    * @param nodeCost needed for bin-packing. It is a cost that assigned to each node based on a cost function
    * @param levelSet Nodes are sorted based on their Level in this array
    * @param levelPtr is the pointer array in CSC format for Levels
    * @param it is an array that for each index i it will return the level of that node/iteration i
    * @param binPacking is a flag that determines whether we want to use bin-packing or not
    */
int get_coarse_Level_set_DAG_CSC03_parallel_NOLEVELSET(
        size_t n, int *lC, int *lR,
        int &finaLevelNo, int*& finaLevelPtr,
        int &partNo, int*& finalPartPtr,
        int*& finalNodePtr,
        int innerParts, int minLevelDist,
        int divRate, int numThreads, double *nodeCost,
        int levelNo, int* levelSet, int* levelPtr,
        int* node2Level, bool binPacking = false);




    /*
    * @brief This function returns the merged levels schedule. Note that it needs the levelset information
    * thus it can reuse one levelset computation in multiple merging schedule
    * @param n Number of nodes in the DAG
    * @param lC is the columns pointer in CSC data structure
    * @param lR is the row index in CSC format
    * @return finalLevelNo is the number of levels in the merged DAG
    * @return finalLevelPtr is the pointer array that show the start and end point of w-partitions in partNo array
    * @return partNo is the number of w-partitions
    * @return finalPartPtr is the w-partitions sorted by their levels
    * @return finalNodePtr is the nodes sorted by their w-partitions
    * @param innerParts Number of inner part (when using Bin Packing) it is normally equal to number of threads
    * @param minLevelDist Is the P2 (Where we should start cutting)
    * @param divRate is the window of merge or P3
    * @param NumThreads is the number of processors
    * @param nodeCost needed for bin-packing. It is a cost that assigned to each node based on a cost function
    * @param levelSet Nodes are sorted based on their Level in this array
    * @param levelPtr is the pointer array in CSC format for Levels
    * @param it is an array that for each index i it will return the level of that node/iteration i
    * @param binPacking is a flag that determines whether we want to use bin-packing or not
    */
    int get_coarse_Level_set_DAG_CSC03_parallel_NOLEVELSET_V2(
            size_t n, int *lC, int *lR,
            int &finaLevelNo, int* finaLevelPtr,
            int &partNo, int* finalPartPtr, int* finalNodePtr,
            int innerParts, int minLevelDist,
            int divRate, int numThreads, double *nodeCost,
            int levelNo, int* levelSet, int* levelPtr,
            int* node2Level, bool binPacking = false);

    int skeletonCoarsening(
            size_t n_nodes, int* DAG_ptr, int* DAG_set,
            int &finaLevelNo, std::vector<int>& finaLevelPtr, int &partNo,
            std::vector<int>& finalPartPtr, std::vector<int>& finalNodePtr,
            std::vector<int>& WM, int numThreads,
            int& levelNo, int* levelSet, int* levelPtr,
            int* node2Level);

 int getCoarseLevelSet_DAG_BCSC02(size_t n,
                                  size_t *lC,
                                  size_t *Li_ptr,
                                  int* lR,
                                  const int *blk2col,
                                  const int *col2blk,
                                  int &finaLevelNo,
                                  int* &finaLevelPtr,
                                  int* &parLevelSet,
                                  int &partNo,
                                  int* &finalPartPtr,
                                  int* &finalNodePtr,
                                  int innerParts,
                                  int minLevelDist,
                                  int divRate,
                                  double *nodeCost);

 int get_coarse_levelSet_tree(size_t n,
                              const int *eTree,
                              int &finaLevelNo,
                              int* &finaLevelPtr,
                              int &partNo,
                              int* &finalPartPtr,
                              int* &finalNodePtr,
                              int innerParts,
                              int minLevelDist,
                              int divRate,
                              double *nodeCost);
}

#endif //PROJECT_LBC_H
