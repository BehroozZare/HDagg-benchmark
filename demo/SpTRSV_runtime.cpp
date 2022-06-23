//
// Created by behrooz on 1/17/21.
//


#include <iostream>
#include <sparse_io.h>
#include <test_utils.h>
#include <omp.h>
#include <metis_interface.h>
#include <csv_utils.h>

#include "SpTRSV_runtime.h"


//#include <spdlog/spdlog.h>
//#include <CLI/CLI.hpp>


using namespace sym_lib;

//struct CLIArgs {
//    bool gen_data = false; ///<@brief if true, it creates a sensitized data
//    std::string mat_address; ///< @brief matrix's address - should be exact address
//    std::string output_address; ///< @brief output's address - should be exact address

//    spdlog::level::level_enum logLevel = spdlog::level::trace;
//    const std::vector<std::pair<std::string, spdlog::level::level_enum>>
//    SPDLOG_LEVEL_NAMES_TO_LEVELS = {
//            { "trace", spdlog::level::trace },
//            { "debug", spdlog::level::debug },
//            { "info", spdlog::level::info },
//            { "warning", spdlog::level::warn },
//            { "error", spdlog::level::err },
//            { "critical", spdlog::level::critical },
//            { "off", spdlog::level::off }
//    };

//    CLIArgs(int argc, char* argv[])
//    {
//        CLI::App app{ "HDAGG" };
//
//        app.add_flag("--GenData", gen_data, "if true, it creates a sensitized data");
//        app.add_option("--Address", mat_address, "matrix's address");
//        app.add_option("-o,--output", output_address, "output folder name");
//
//        try {
//            app.parse(argc, argv);
//        }
//        catch (const CLI::ParseError& e) {
//            exit(app.exit(e));
//        }
//
//    }
//};





int main(int argc, char *argv[]) {
//    CLIArgs args(argc, argv);

    CSC *Lower_A_CSC, *A = NULLPNTR;
    CSR *Lower_A_CSR;
    size_t n, nnz;
    int *perm;
    std::string matrix_name;
    std::vector<timing_measurement> time_array;

    if (argc < 2) {
        PRINT_LOG("Not enough input args, switching to random mode.\n");
        n = 50;
        double density = 0.3;
        matrix_name = "Random_" + std::to_string(n);
        A = random_square_sparse(n, density);
        if (A == NULLPNTR)
            return -1;
        Lower_A_CSC = make_half(A->n, A->p, A->i, A->x);
        delete A;

        //        std::vector<std::vector<bool>> DAG(n, std::vector<bool>(n,false));
        //        addEdge(DAG, 0, 1);
        //        addEdge(DAG, 0, 3);
        //        addEdge(DAG, 0, 4);
        //        addEdge(DAG, 1, 2);
        //        addEdge(DAG, 3, 4);
        //        addEdge(DAG, 3, 8);
        //        addEdge(DAG, 3, 11);
        ////        addEdge(DAG, 4, 11);
        //        addEdge(DAG, 5, 7);
        //        addEdge(DAG, 5, 8);
        //        addEdge(DAG, 6, 7);
        //        addEdge(DAG, 6, 8);
        //        addEdge(DAG, 7, 8);
        //        addEdge(DAG, 7, 9);
        //        addEdge(DAG, 8, 9);
        //        addEdge(DAG, 8, 10);
        //        addEdge(DAG, 8, 11);
        //        addEdge(DAG, 8, 12);
        //        addEdge(DAG, 11, 12);
        //
        //        Lower_A_CSC = generateDAGfromEdgeList(n, DAG);
        //        std::vector<int> DAG_ptr(Lower_A_CSC->p, Lower_A_CSC->p + n + 1);
        //        std::vector<int> DAG_set(Lower_A_CSC->i, Lower_A_CSC->i + Lower_A_CSC->nnz);
        nnz = Lower_A_CSC->nnz;
    } else {
        std::cout << "Reading matrix " << argv[1] << std::endl;
        std::string f1 = argv[1];
        matrix_name = f1;
        Lower_A_CSC = read_mtx(f1);
        if (Lower_A_CSC == NULLPNTR) {
            std::cout << "Reading failed" << std::endl;
            return -1;
        }
        n = Lower_A_CSC->n;
    }

    if (argc >= 3) {
        std::cout << "Un acceptable argument .. it was num_threads previously" << std::endl;
    }

    /// Re-ordering matrix A
    //    std::cout << "METIS IS NOT ACTIVATED" << std::endl;
    if (true) {
        std::cout << "METIS IS ACTIVATED" << std::endl;
        //We only reorder A since dependency matters more in l-solve.
        A = make_full(Lower_A_CSC);
        delete Lower_A_CSC;
        metis_perm_general(A, perm);
        Lower_A_CSC = make_half(A->n, A->p, A->i, A->x);
        CSC *Lt = transpose_symmetric(Lower_A_CSC, perm);
        CSC *L1_ord = transpose_symmetric(Lt, NULLPNTR);
        delete Lower_A_CSC;
        Lower_A_CSC = L1_ord;
        Lower_A_CSR = csc_to_csr(Lower_A_CSC);
        delete Lt;
        delete A;
        delete[] perm;
    } else {
        CSC *tmp = make_half(Lower_A_CSC->n, Lower_A_CSC->p, Lower_A_CSC->i, Lower_A_CSC->x);
        delete Lower_A_CSC;
        Lower_A_CSC = tmp;
        Lower_A_CSR = csc_to_csr(Lower_A_CSC);
    }
    nnz = Lower_A_CSC->nnz;

    bool isLfactor = false;

    double *y_serial, *y_correct = new double[n], *y_blocked_correct = new double[n];
    double *y_serial_perfect, *y_correct_perfect = new double[n];
    std::cout << "Starting SpTrSv Runtime analysis" << std::endl;

    std::vector<std::string> Runtime_headers;
    Runtime_headers.emplace_back("Matrix_Name");
    Runtime_headers.emplace_back("Algorithm");
    Runtime_headers.emplace_back("Kernel");
    Runtime_headers.emplace_back("Core");
    Runtime_headers.emplace_back("Scheduling_Time");
    Runtime_headers.emplace_back("Executor_Runtime");
    Runtime_headers.emplace_back("FLOPS");
    Runtime_headers.emplace_back("MemTraffic");
    Runtime_headers.emplace_back("OI");
    Runtime_headers.emplace_back("wm");
    Runtime_headers.emplace_back("nlevel");
    Runtime_headers.emplace_back("PG");
    Runtime_headers.emplace_back("Profitable");

    if (ANALYSE_STATISTIC) {
        Runtime_headers.emplace_back("Parallelism");
        Runtime_headers.emplace_back("MaximalDeference");
        Runtime_headers.emplace_back("Deviation");
        Runtime_headers.emplace_back("Parts");
    }


    std::string delimiter = "/";
    size_t matrix_name_start_pos = matrix_name.find(delimiter) + 1;
    delimiter = ".";
    size_t matrix_name_end_pos = matrix_name.find(delimiter);
    if (argc < 2) {
        matrix_name_start_pos = 0;
        matrix_name_end_pos = matrix_name.size();
    }
    std::string Mat_name(matrix_name.begin() + matrix_name_start_pos, matrix_name.begin() + matrix_name_end_pos);
    std::string Data_name = "SpTrSv_Final_20";
    profiling_utils::CSVManager runtime_csv(Data_name,
                                            "some address", Runtime_headers,
                                            false);
#ifdef ANALYSE_STATISTIC
    std::vector<std::string> Parameter_headers;
    Parameter_headers.emplace_back("Matrix_Name");
    Parameter_headers.emplace_back("Algorithm");
    Parameter_headers.emplace_back("Kernel");
    Parameter_headers.emplace_back("Core");
    Parameter_headers.emplace_back("AvgParallelism");
    Parameter_headers.emplace_back("WeightedAvgParallelism");
    Parameter_headers.emplace_back("AvgGroupSize");

    Parameter_headers.emplace_back("AvgCostPerNNZSpTrSv");
    Parameter_headers.emplace_back("CostPerVertexSpTrSv");

    Parameter_headers.emplace_back("AvgCostPerNNZSpIC0");
    Parameter_headers.emplace_back("CostPerVertexSpIC0");

    Parameter_headers.emplace_back("CostPerVertexSpILU0");
    Parameter_headers.emplace_back("AvgCostPerNNZSpILU0");

    Parameter_headers.emplace_back("PG");


    profiling_utils::CSVManager parameter_csv("SpTrSv_Analysis",
                                              "some address", Parameter_headers,
                                              false);
#endif


    auto tmp = make_full(Lower_A_CSC);
    auto CSR_A = csc_to_csr(tmp);
    std::vector<int> Cores{64};

    std::vector<bool> choice{true};
    std::vector<bool> bin_pack{true};
    std::vector<int> WM{1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28};
    //====================================== Left Looking Section ==========================================
    //"********************* LL Serial  *********************"
    Sptrsv_LL_FLOPS_MEMORY LL_FLOPS(Lower_A_CSR, Lower_A_CSC, NULLPNTR, "LL FLOPS serial"); //seq
    auto LL_FLOPS_runtime = LL_FLOPS.evaluate();
    //"********************* LL Serial  *********************"
    Sptrsv_LL_Serial LL_serial(Lower_A_CSR, Lower_A_CSC, NULLPNTR, "LL serial"); //seq
    auto LL_serial_runtime = LL_serial.evaluate();
    y_serial = LL_serial.solution();
    copy_vector(0, n, y_serial, y_correct);
    std::cout << "Running LL Serial Code - The runtime:" << LL_serial_runtime.elapsed_time << std::endl;
    runtime_csv.addElementToRecord(Mat_name, "Matrix_Name");
    runtime_csv.addElementToRecord("Serial", "Algorithm");
    runtime_csv.addElementToRecord("LL", "Kernel");
    runtime_csv.addElementToRecord(1, "Core");
    runtime_csv.addElementToRecord(0, "Scheduling_Time");
    runtime_csv.addElementToRecord(LL_serial_runtime.elapsed_time, "Executor_Runtime");
    runtime_csv.addElementToRecord(LL_FLOPS.getFLOPS() / LL_serial_runtime.elapsed_time, "FLOPS");
    runtime_csv.addElementToRecord(LL_FLOPS.getBytes(), "MemTraffic");
    runtime_csv.addElementToRecord(LL_FLOPS.getFLOPS() / LL_FLOPS.getBytes(), "OI");
    runtime_csv.addElementToRecord(1, "wm");
    runtime_csv.addElementToRecord(1, "nlevel");
    runtime_csv.addElementToRecord(0, "Profitable");
    runtime_csv.addElementToRecord(1, "PG");

    if (ANALYSE_STATISTIC) {
        runtime_csv.addElementToRecord(1, "Parallelism");
        runtime_csv.addElementToRecord(1, "MaximalDeference");
        runtime_csv.addElementToRecord(1, "Deviation");
        runtime_csv.addElementToRecord(Lower_A_CSC->n, "Parts");
    }

    runtime_csv.addRecord();

    //"********************* LL Levelset *********************"
    for (auto &core: Cores) {
        SpTrsv_LL_LevelSet LL_lvl_obj(Lower_A_CSR, Lower_A_CSC, y_correct,
                                      "LL Levelset ", core);
        timing_measurement LL_lvl_runtime;
        omp_set_num_threads(core);
        LL_lvl_runtime = LL_lvl_obj.evaluate();
        std::cout << "Running LL Levelset Code with #core: " << core <<
        " - The runtime:" << LL_lvl_runtime.elapsed_time << std::endl;
        int part_no, nlevels;
        LL_lvl_obj.getWaveStatistic(nlevels, part_no);
        runtime_csv.addElementToRecord(Mat_name, "Matrix_Name");
        runtime_csv.addElementToRecord("Levelset", "Algorithm");
        runtime_csv.addElementToRecord("LL", "Kernel");
        runtime_csv.addElementToRecord(core, "Core");
        runtime_csv.addElementToRecord(LL_lvl_obj.getSchedulingTime(), "Scheduling_Time");
        runtime_csv.addElementToRecord(LL_lvl_runtime.elapsed_time, "Executor_Runtime");
        runtime_csv.addElementToRecord(LL_FLOPS.getFLOPS() / LL_lvl_runtime.elapsed_time, "FLOPS");
        runtime_csv.addElementToRecord(LL_FLOPS.getBytes(), "MemTraffic");
        runtime_csv.addElementToRecord(LL_FLOPS.getFLOPS() / LL_FLOPS.getBytes(), "OI");
        runtime_csv.addElementToRecord(1, "wm");
        runtime_csv.addElementToRecord(nlevels, "nlevel");
        double profitable = (LL_lvl_obj.getSchedulingTime()) /
                (LL_serial_runtime.elapsed_time - LL_lvl_runtime.elapsed_time);
        runtime_csv.addElementToRecord(profitable, "Profitable");
        runtime_csv.addElementToRecord(1, "PG");

        if (ANALYSE_STATISTIC) {
            double avg_par = part_no * 1.0 / nlevels;
            runtime_csv.addElementToRecord(avg_par, "Parallelism");
            runtime_csv.addElementToRecord(0, "MaximalDeference");
            runtime_csv.addElementToRecord(0, "Deviation");
            runtime_csv.addElementToRecord(part_no, "Parts");

        }

        runtime_csv.addRecord();
    }
    //"********************* LL Levelset *********************"
    for (auto &core: Cores) {
        SpTrsv_LL_Parallel_Levelset LL_lvl_obj(CSR_A, Lower_A_CSR, Lower_A_CSC, y_correct,
                                               "LL Levelset ", core);
        timing_measurement LL_lvl_runtime;
        omp_set_num_threads(core);
        LL_lvl_runtime = LL_lvl_obj.evaluate();
        std::cout << "Running LL parallel Levelset Code with #core: " << core <<
        " - The runtime:" << LL_lvl_runtime.elapsed_time << std::endl;
        int part_no, nlevels;
        LL_lvl_obj.getWaveStatistic(nlevels, part_no);
        runtime_csv.addElementToRecord(Mat_name, "Matrix_Name");
        runtime_csv.addElementToRecord("Parallel_Levelset", "Algorithm");
        runtime_csv.addElementToRecord("LL", "Kernel");
        runtime_csv.addElementToRecord(core, "Core");
        runtime_csv.addElementToRecord(LL_lvl_obj.getSchedulingTime(), "Scheduling_Time");
        runtime_csv.addElementToRecord(LL_lvl_runtime.elapsed_time, "Executor_Runtime");
        runtime_csv.addElementToRecord(LL_FLOPS.getFLOPS() / LL_lvl_runtime.elapsed_time, "FLOPS");
        runtime_csv.addElementToRecord(LL_FLOPS.getBytes(), "MemTraffic");
        runtime_csv.addElementToRecord(LL_FLOPS.getFLOPS() / LL_FLOPS.getBytes(), "OI");
        runtime_csv.addElementToRecord(1, "wm");
        runtime_csv.addElementToRecord(nlevels, "nlevel");
        double profitable = (LL_lvl_obj.getSchedulingTime()) /
                (LL_serial_runtime.elapsed_time - LL_lvl_runtime.elapsed_time);
        runtime_csv.addElementToRecord(profitable, "Profitable");
        runtime_csv.addElementToRecord(1, "PG");

        if (ANALYSE_STATISTIC) {
            double avg_par = part_no * 1.0 / nlevels;
            runtime_csv.addElementToRecord(avg_par, "Parallelism");
            runtime_csv.addElementToRecord(0, "MaximalDeference");
            runtime_csv.addElementToRecord(0, "Deviation");
            runtime_csv.addElementToRecord(part_no, "Parts");

        }

        runtime_csv.addRecord();
    }
    //"********************* LL Tree + Levelset *********************"
    for (auto &core: Cores) {
        SpTrSv_LL_Tree_Levelset_No_unPACK LL_lvl_obj(Lower_A_CSR, Lower_A_CSC, y_correct,
                                                     "LL Tree Levelset ", core);
        timing_measurement LL_lvl_runtime;
        omp_set_num_threads(core);
        LL_lvl_runtime = LL_lvl_obj.evaluate();
        std::cout << "Running LL Levelset + Tree Code with #core: " << core <<
        " - The runtime:" << LL_lvl_runtime.elapsed_time << std::endl;
        int part_no, nlevels;
        LL_lvl_obj.getWaveStatistic(nlevels, part_no);
        runtime_csv.addElementToRecord(Mat_name, "Matrix_Name");
        runtime_csv.addElementToRecord("Tree_Levelset", "Algorithm");
        runtime_csv.addElementToRecord("LL", "Kernel");
        runtime_csv.addElementToRecord(core, "Core");
        runtime_csv.addElementToRecord(LL_lvl_obj.getSchedulingTime(), "Scheduling_Time");
        runtime_csv.addElementToRecord(LL_lvl_runtime.elapsed_time, "Executor_Runtime");
        runtime_csv.addElementToRecord(LL_FLOPS.getFLOPS() / LL_lvl_runtime.elapsed_time, "FLOPS");
        runtime_csv.addElementToRecord(LL_FLOPS.getBytes(), "MemTraffic");
        runtime_csv.addElementToRecord(LL_FLOPS.getFLOPS() / LL_FLOPS.getBytes(), "OI");
        runtime_csv.addElementToRecord(1, "wm");
        runtime_csv.addElementToRecord(nlevels, "nlevel");
        double profitable = (LL_lvl_obj.getSchedulingTime()) /
                (LL_serial_runtime.elapsed_time - LL_lvl_runtime.elapsed_time);
        runtime_csv.addElementToRecord(profitable, "Profitable");
        runtime_csv.addElementToRecord(1, "PG");

        if (ANALYSE_STATISTIC) {
            double avg_par = part_no * 1.0 / nlevels;
            runtime_csv.addElementToRecord(avg_par, "Parallelism");
            runtime_csv.addElementToRecord(0, "MaximalDeference");
            runtime_csv.addElementToRecord(0, "Deviation");
            runtime_csv.addElementToRecord(part_no, "Parts");

        }

        runtime_csv.addRecord();
    }
    //"********************* LL GLC *********************"
    for (auto &core: Cores) {
        for (auto &&bin: bin_pack) {
            SpTrSv_LL_GLC LL_lvl_obj(Lower_A_CSR, Lower_A_CSC, y_correct,
                                     "LL GLC ", core, bin);
            timing_measurement LL_lvl_runtime;
            omp_set_num_threads(core);
            LL_lvl_runtime = LL_lvl_obj.evaluate();
            if (bin) {
                std::cout << "Running LL GLC with BIN Code with #core: "
                << core << " - The runtime:" << LL_lvl_runtime.elapsed_time << std::endl;
            } else {
                std::cout << "Running LL GLC Code with #core: "
                << core << " - The runtime:" << LL_lvl_runtime.elapsed_time << std::endl;
            }
            int part_no, nlevels;
            LL_lvl_obj.getWaveStatistic(nlevels, part_no);
            double avg_par = part_no * 1.0 / nlevels;
            runtime_csv.addElementToRecord(Mat_name, "Matrix_Name");
            if (bin) {
                runtime_csv.addElementToRecord("GLC_BIN", "Algorithm");
            } else {
                runtime_csv.addElementToRecord("GLC_NO_BIN", "Algorithm");
            }
            runtime_csv.addElementToRecord("LL", "Kernel");
            runtime_csv.addElementToRecord(core, "Core");
            runtime_csv.addElementToRecord(LL_lvl_obj.getSchedulingTime(), "Scheduling_Time");
            runtime_csv.addElementToRecord(LL_lvl_runtime.elapsed_time, "Executor_Runtime");
            runtime_csv.addElementToRecord(LL_FLOPS.getFLOPS() / LL_lvl_runtime.elapsed_time, "FLOPS");
            runtime_csv.addElementToRecord(LL_FLOPS.getBytes(), "MemTraffic");
            runtime_csv.addElementToRecord(LL_FLOPS.getFLOPS() / LL_FLOPS.getBytes(), "OI");
            runtime_csv.addElementToRecord(1, "wm");
            runtime_csv.addElementToRecord(nlevels, "nlevel");
            double profitable = (LL_lvl_obj.getSchedulingTime()) /
                    (LL_serial_runtime.elapsed_time - LL_lvl_runtime.elapsed_time);
            runtime_csv.addElementToRecord(profitable, "Profitable");
            runtime_csv.addElementToRecord(1, "PG");


            if (ANALYSE_STATISTIC) {
                runtime_csv.addElementToRecord(avg_par, "Parallelism");
                runtime_csv.addElementToRecord(1, "MaximalDeference");
                runtime_csv.addElementToRecord(1, "Deviation");
                runtime_csv.addElementToRecord(part_no, "Parts");
            }

            runtime_csv.addRecord();
        }
    }
    //"********************* LL Tree GLC *********************"
    for (auto &core: Cores) {
        for (auto &&bin: bin_pack) {
            SpTrSv_LL_Tree_GLC LL_lvl_obj(Lower_A_CSR, Lower_A_CSC, y_correct,
                                          "LL Tree GLC ", core, isLfactor, bin);
            timing_measurement LL_lvl_runtime;
            omp_set_num_threads(core);
            LL_lvl_runtime = LL_lvl_obj.evaluate();
            if (bin) {
                std::cout << "Running LL Tree GLC with BIN Code with #core: "
                << core << " - The runtime:" << LL_lvl_runtime.elapsed_time << std::endl;
            } else {
                std::cout << "Running LL Tree GLC Code with #core: "
                << core << " - The runtime:" << LL_lvl_runtime.elapsed_time << std::endl;
            }
            int part_no, nlevels;
            LL_lvl_obj.getWaveStatistic(nlevels, part_no);
            runtime_csv.addElementToRecord(Mat_name, "Matrix_Name");
            if (bin) {
                runtime_csv.addElementToRecord("Tree_GLC_BIN", "Algorithm");
            } else {
                runtime_csv.addElementToRecord("Tree_GLC_NO_BIN", "Algorithm");
            }
            runtime_csv.addElementToRecord("LL", "Kernel");
            runtime_csv.addElementToRecord(core, "Core");
            runtime_csv.addElementToRecord(LL_lvl_obj.getSchedulingTime(), "Scheduling_Time");
            runtime_csv.addElementToRecord(LL_lvl_runtime.elapsed_time, "Executor_Runtime");
            runtime_csv.addElementToRecord(LL_FLOPS.getFLOPS() / LL_lvl_runtime.elapsed_time, "FLOPS");
            runtime_csv.addElementToRecord(LL_FLOPS.getBytes(), "MemTraffic");
            runtime_csv.addElementToRecord(LL_FLOPS.getFLOPS() / LL_FLOPS.getBytes(), "OI");
            runtime_csv.addElementToRecord(1, "wm");
            runtime_csv.addElementToRecord(nlevels, "nlevel");
            double profitable = (LL_lvl_obj.getSchedulingTime()) /
                    (LL_serial_runtime.elapsed_time - LL_lvl_runtime.elapsed_time);
            runtime_csv.addElementToRecord(profitable, "Profitable");
            runtime_csv.addElementToRecord(1, "PG");


            if (ANALYSE_STATISTIC) {
                double var, avg_par, max_diff;
                LL_lvl_obj.getStatistic(avg_par, max_diff, var);
                runtime_csv.addElementToRecord(avg_par, "Parallelism");
                runtime_csv.addElementToRecord(max_diff, "MaximalDeference");
                runtime_csv.addElementToRecord(var, "Deviation");
                runtime_csv.addElementToRecord(part_no, "Parts");
            }

            runtime_csv.addRecord();
        }
    }
    //"********************* LL Tree BFS GLC *********************"
    for (auto &core: Cores) {
        for (auto &&bin: {true}) {
            SpTrSv_LL_Tree_GLC_BFS LL_lvl_obj(Lower_A_CSR, Lower_A_CSC, y_correct,
                                              "LL Tree BFS GLC ", core, isLfactor, bin);
            timing_measurement LL_lvl_runtime;
            omp_set_num_threads(core);
            LL_lvl_runtime = LL_lvl_obj.evaluate();
            if (bin) {
                std::cout << "Running LL Tree BFS GLC with BIN Code with #core: "
                << core << " - The runtime:" << LL_lvl_runtime.elapsed_time << std::endl;
            } else {
                std::cout << "Running LL Tree BFS GLC Code with #core: "
                << core << " - The runtime:" << LL_lvl_runtime.elapsed_time << std::endl;
            }
            int part_no, nlevels;
            LL_lvl_obj.getWaveStatistic(nlevels, part_no);

            runtime_csv.addElementToRecord(Mat_name, "Matrix_Name");
            if (bin) {
                runtime_csv.addElementToRecord("Tree_GLC_BFS_BIN", "Algorithm");
            } else {
                runtime_csv.addElementToRecord("Tree_GLC_BFS_NO_BIN", "Algorithm");
            }
            runtime_csv.addElementToRecord("LL", "Kernel");
            runtime_csv.addElementToRecord(core, "Core");
            runtime_csv.addElementToRecord(LL_lvl_obj.getSchedulingTime(), "Scheduling_Time");
            runtime_csv.addElementToRecord(LL_lvl_runtime.elapsed_time, "Executor_Runtime");
            runtime_csv.addElementToRecord(LL_FLOPS.getFLOPS() / LL_lvl_runtime.elapsed_time, "FLOPS");
            runtime_csv.addElementToRecord(LL_FLOPS.getBytes(), "MemTraffic");
            runtime_csv.addElementToRecord(LL_FLOPS.getFLOPS() / LL_FLOPS.getBytes(), "OI");
            runtime_csv.addElementToRecord(1, "wm");
            runtime_csv.addElementToRecord(nlevels, "nlevel");
            double profitable = (LL_lvl_obj.getSchedulingTime()) /
                    (LL_serial_runtime.elapsed_time - LL_lvl_runtime.elapsed_time);
            runtime_csv.addElementToRecord(profitable, "Profitable");
            runtime_csv.addElementToRecord(LL_lvl_obj.getPotentialGain(), "PG");
            runtime_csv.addRecord();

            if (ANALYSE_STATISTIC) {
                double var, avg_par, max_diff;
                LL_lvl_obj.getStatistic(avg_par, max_diff, var);

                runtime_csv.addElementToRecord(avg_par, "Parallelism");
                runtime_csv.addElementToRecord(max_diff, "MaximalDeference");
                runtime_csv.addElementToRecord(var, "Deviation");
                runtime_csv.addElementToRecord(part_no, "Parts");


                double avgParallelism, WeightedAvgParallelism, CostAvgParallelism,
                AvgGroupSize, PG;
                double AvgCostPerVertexSpTrSv, AvgCostPerNNZSpTrSv;
                double AvgCostPerVertexSpILU0, AvgCostPerNNZSpLU0;
                double AvgCostPerVertexSpIC0, AvgCostPerNNZSpIC0;

                LL_lvl_obj.computeMetrics(avgParallelism, WeightedAvgParallelism,
                                          CostAvgParallelism, AvgGroupSize,
                                          AvgCostPerVertexSpTrSv, AvgCostPerNNZSpTrSv,
                                          AvgCostPerVertexSpILU0, AvgCostPerNNZSpLU0,
                                          AvgCostPerVertexSpIC0, AvgCostPerNNZSpIC0,
                                          PG);
                parameter_csv.addElementToRecord(Mat_name, "Matrix_Name");
                if (bin) {
                    parameter_csv.addElementToRecord("Tree_GLC_BIN", "Algorithm");
                } else {
                    parameter_csv.addElementToRecord("Tree_GLC_NO_BIN", "Algorithm");
                }
                parameter_csv.addElementToRecord(1, "Kernel");
                parameter_csv.addElementToRecord(core, "Core");
                parameter_csv.addElementToRecord(avgParallelism, "AvgParallelism");
                parameter_csv.addElementToRecord(WeightedAvgParallelism, "WeightedAvgParallelism");
                parameter_csv.addElementToRecord(AvgGroupSize, "AvgGroupSize");

                parameter_csv.addElementToRecord(AvgCostPerVertexSpTrSv, "CostPerVertexSpTrSv");
                parameter_csv.addElementToRecord(AvgCostPerNNZSpTrSv, "AvgCostPerNNZSpTrSv");

                parameter_csv.addElementToRecord(AvgCostPerVertexSpIC0, "CostPerVertexSpIC0");
                parameter_csv.addElementToRecord(AvgCostPerNNZSpIC0, "AvgCostPerNNZSpIC0");

                parameter_csv.addElementToRecord(AvgCostPerVertexSpILU0, "CostPerVertexSpILU0");
                parameter_csv.addElementToRecord(AvgCostPerNNZSpLU0, "AvgCostPerNNZSpILU0");

                parameter_csv.addElementToRecord(PG, "PG");
                parameter_csv.addRecord();
            }

        }
    }
    //"********************* LL Tree BFS GLC P2P *********************"
    for (auto &core: Cores) {
        for (auto &&bin: {false}) {
            omp_set_num_threads(core);
            SpTrSv_LL_Tree_GLC_BFS_P2P LL_lvl_obj(Lower_A_CSR, Lower_A_CSC, y_correct,
                                                  "LL Tree BFS GLC P2P ", core, isLfactor, bin);
            timing_measurement LL_lvl_runtime;
            LL_lvl_runtime = LL_lvl_obj.evaluate();
            if (bin) {
                std::cout << "Running LL Tree BFS GLC P2P with BIN Code with #core: "
                << core << " - The runtime:" << LL_lvl_runtime.elapsed_time << std::endl;
            } else {
                std::cout << "Running LL Tree BFS GLC P2P Code with #core: "
                << core << " - The runtime:" << LL_lvl_runtime.elapsed_time << std::endl;
            }
            int part_no, nlevels;
            LL_lvl_obj.getWaveStatistic(nlevels, part_no);

            runtime_csv.addElementToRecord(Mat_name, "Matrix_Name");
            if (bin) {
                runtime_csv.addElementToRecord("Tree_GLC_BFS_P2P_BIN", "Algorithm");
            } else {
                runtime_csv.addElementToRecord("Tree_GLC_BFS_P2P_NO_BIN", "Algorithm");
            }
            runtime_csv.addElementToRecord("LL", "Kernel");
            runtime_csv.addElementToRecord(core, "Core");
            runtime_csv.addElementToRecord(LL_lvl_obj.getSchedulingTime(), "Scheduling_Time");
            runtime_csv.addElementToRecord(LL_lvl_runtime.elapsed_time, "Executor_Runtime");
            runtime_csv.addElementToRecord(LL_FLOPS.getFLOPS() / LL_lvl_runtime.elapsed_time, "FLOPS");
            runtime_csv.addElementToRecord(LL_FLOPS.getBytes(), "MemTraffic");
            runtime_csv.addElementToRecord(LL_FLOPS.getFLOPS() / LL_FLOPS.getBytes(), "OI");
            runtime_csv.addElementToRecord(1, "wm");
            runtime_csv.addElementToRecord(nlevels, "nlevel");
            double profitable = (LL_lvl_obj.getSchedulingTime()) /
                    (LL_serial_runtime.elapsed_time - LL_lvl_runtime.elapsed_time);
            runtime_csv.addElementToRecord(profitable, "Profitable");
            runtime_csv.addElementToRecord(1, "PG");

            if (ANALYSE_STATISTIC) {
                double var, avg_par, max_diff;
                LL_lvl_obj.getStatistic(avg_par, max_diff, var);

                runtime_csv.addElementToRecord(avg_par, "Parallelism");
                runtime_csv.addElementToRecord(max_diff, "MaximalDeference");
                runtime_csv.addElementToRecord(var, "Deviation");
                runtime_csv.addElementToRecord(part_no, "Parts");

            }

            runtime_csv.addRecord();
        }
    }
    //"********************* LL LBC *********************"
    std::vector<int> P3 = {5000};
    for (auto core: Cores) {
        for (auto p3: P3) {
            SpTrSv_LL_LBC LL_lvl_obj(Lower_A_CSR, Lower_A_CSC, y_correct, "LL LBC_Tree", core);
            omp_set_num_threads(core);
            LL_lvl_obj.setP2_P3(-1, p3);
            auto LL_lvl_runtime = LL_lvl_obj.evaluate();
            std::cout << "Running LL LBC_Tree Code with #core: " << core << " and P3=" << p3 <<
            " - The runtime:" << LL_lvl_runtime.elapsed_time << std::endl;
            int part_no, nlevels;
            LL_lvl_obj.getWaveStatistic(nlevels, part_no);

            runtime_csv.addElementToRecord(Mat_name, "Matrix_Name");
            runtime_csv.addElementToRecord("LBC", "Algorithm");
            runtime_csv.addElementToRecord("LL", "Kernel");
            runtime_csv.addElementToRecord(core, "Core");
            runtime_csv.addElementToRecord(LL_lvl_obj.getSchedulingTime(), "Scheduling_Time");
            runtime_csv.addElementToRecord(LL_lvl_runtime.elapsed_time, "Executor_Runtime");
            runtime_csv.addElementToRecord(LL_FLOPS.getFLOPS() / LL_lvl_runtime.elapsed_time, "FLOPS");
            runtime_csv.addElementToRecord(LL_FLOPS.getBytes(), "MemTraffic");
            runtime_csv.addElementToRecord(LL_FLOPS.getFLOPS() / LL_FLOPS.getBytes(), "OI");
            runtime_csv.addElementToRecord(p3, "wm");
            runtime_csv.addElementToRecord(nlevels, "nlevel");
            double profitable = (LL_lvl_obj.getSchedulingTime()) /
                    (LL_serial_runtime.elapsed_time - LL_lvl_runtime.elapsed_time);
            runtime_csv.addElementToRecord(profitable, "Profitable");
            runtime_csv.addElementToRecord(1, "PG");

            if (ANALYSE_STATISTIC) {
                double var, avg_par, max_diff;
                LL_lvl_obj.getStatistic(avg_par, max_diff, var);

                runtime_csv.addElementToRecord(avg_par, "Parallelism");
                runtime_csv.addElementToRecord(max_diff, "MaximalDeference");
                runtime_csv.addElementToRecord(var, "Deviation");
                runtime_csv.addElementToRecord(part_no, "Parts");

            }
            runtime_csv.addRecord();
        }
    }
    //"********************* MKL *********************"
    for (auto &core: Cores) {
        omp_set_num_threads(core);
        SpTrSv_LL_MKL LL_lvl_obj(Lower_A_CSR, Lower_A_CSC, y_correct, "LL MKL", core);
        timing_measurement LL_lvl_runtime;
        LL_lvl_runtime = LL_lvl_obj.evaluate();
        std::cout << "Running LL MKL Code with #core: " << core <<
        " - The runtime:" << LL_lvl_runtime.elapsed_time << std::endl;
        runtime_csv.addElementToRecord(Mat_name, "Matrix_Name");
        runtime_csv.addElementToRecord("MKL", "Algorithm");
        runtime_csv.addElementToRecord("LL", "Kernel");
        runtime_csv.addElementToRecord(core, "Core");
        runtime_csv.addElementToRecord(LL_lvl_obj.getSchedulingTime(), "Scheduling_Time");
        runtime_csv.addElementToRecord(LL_lvl_runtime.elapsed_time, "Executor_Runtime");
        runtime_csv.addElementToRecord(LL_FLOPS.getFLOPS() / LL_lvl_runtime.elapsed_time, "FLOPS");
        runtime_csv.addElementToRecord(LL_FLOPS.getBytes(), "MemTraffic");
        runtime_csv.addElementToRecord(LL_FLOPS.getFLOPS() / LL_FLOPS.getBytes(), "OI");
        runtime_csv.addElementToRecord(1, "wm");
        double profitable = (LL_lvl_obj.getSchedulingTime()) /
                (LL_serial_runtime.elapsed_time - LL_lvl_runtime.elapsed_time);
        runtime_csv.addElementToRecord(profitable, "Profitable");
        runtime_csv.addElementToRecord(1, "nlevel");
        runtime_csv.addElementToRecord(1, "PG");
        if (ANALYSE_STATISTIC) {
            runtime_csv.addElementToRecord(1, "Parallelism");
            runtime_csv.addElementToRecord(1, "MaximalDeference");
            runtime_csv.addElementToRecord(1, "Deviation");
            runtime_csv.addElementToRecord(1, "Parts");
        }

        runtime_csv.addRecord();
    }
    //    //"********************* DAGP *********************"
    //    for (auto &core: Cores) {
    //        std::vector<int> Parts{1000};
    //        for (auto &part: Parts) {
    //            SpTrSv_LL_DAGP LL_lvl_obj(Lower_A_CSR, Lower_A_CSC, y_correct, "LL DAGP", part, core);
    //            omp_set_num_threads(core);
    //            timing_measurement LL_lvl_runtime;
    //            LL_lvl_runtime = LL_lvl_obj.evaluate();
    //            std::cout << "Running LL DAGP Code with #core: " << core <<
    //                      " - The runtime:" << LL_lvl_runtime.elapsed_time << std::endl;
    //            int part_no, nlevels;
    //            LL_lvl_obj.getWaveStatistic(nlevels, part_no);
    //
    //            runtime_csv.addElementToRecord(Mat_name, "Matrix_Name");
    //            runtime_csv.addElementToRecord("DAGP", "Algorithm");
    //            runtime_csv.addElementToRecord("LL", "Kernel");
    //            runtime_csv.addElementToRecord(core, "Core");
    //            runtime_csv.addElementToRecord(LL_lvl_obj.getSchedulingTime(), "Scheduling_Time");
    //            runtime_csv.addElementToRecord(LL_lvl_runtime.elapsed_time, "Executor_Runtime");
    //            runtime_csv.addElementToRecord(LL_FLOPS.getFLOPS() / LL_lvl_runtime.elapsed_time, "FLOPS");
    //            runtime_csv.addElementToRecord(LL_FLOPS.getBytes(), "MemTraffic");
    //            runtime_csv.addElementToRecord(LL_FLOPS.getFLOPS() / LL_FLOPS.getBytes(), "OI");
    //            runtime_csv.addElementToRecord(part, "wm");
    //            runtime_csv.addElementToRecord(nlevels, "nlevel");
    //
    //            double profitable = (LL_lvl_obj.getSchedulingTime()) /
    //                                (LL_serial_runtime.elapsed_time - LL_lvl_runtime.elapsed_time);
    //            runtime_csv.addElementToRecord(profitable, "Profitable");
    //            runtime_csv.addElementToRecord(1, "PG");
    //
    //            if (ANALYSE_STATISTIC) {
    //                double var, avg_par, max_diff;
    //                LL_lvl_obj.getStatistic(avg_par, max_diff, var);
    //
    //                runtime_csv.addElementToRecord(avg_par, "Parallelism");
    //                runtime_csv.addElementToRecord(max_diff, "MaximalDeference");
    //                runtime_csv.addElementToRecord(var, "Deviation");
    //                runtime_csv.addElementToRecord(part_no, "Parts");
    //
    //            }
    //
    //            runtime_csv.addRecord();
    //        }
    //    }
    //"********************* SpMP *********************"
    for (auto &core: Cores) {
        omp_set_num_threads(core);
        SpTrsv_LL_SpMP LL_lvl_obj(CSR_A, Lower_A_CSR, Lower_A_CSC, y_correct, "LL SpMP", core); //seq
        auto LL_lvl_runtime = LL_lvl_obj.evaluate();
        LL_lvl_runtime = LL_lvl_obj.evaluate();
        std::cout << "Running LL SpMP Code with #core: " << core <<
        " - The runtime:" << LL_lvl_runtime.elapsed_time << std::endl;
        runtime_csv.addElementToRecord(Mat_name, "Matrix_Name");
        runtime_csv.addElementToRecord("SpMP", "Algorithm");
        runtime_csv.addElementToRecord("LL", "Kernel");
        runtime_csv.addElementToRecord(core, "Core");
        runtime_csv.addElementToRecord(LL_lvl_obj.getSchedulingTime(), "Scheduling_Time");
        runtime_csv.addElementToRecord(LL_lvl_runtime.elapsed_time, "Executor_Runtime");
        runtime_csv.addElementToRecord(LL_FLOPS.getFLOPS() / LL_lvl_runtime.elapsed_time, "FLOPS");
        runtime_csv.addElementToRecord(LL_FLOPS.getBytes(), "MemTraffic");
        runtime_csv.addElementToRecord(LL_FLOPS.getFLOPS() / LL_FLOPS.getBytes(), "OI");
        runtime_csv.addElementToRecord(1, "wm");
        runtime_csv.addElementToRecord(LL_lvl_obj.getNumberOfP2P(), "nlevel");
        double profitable = (LL_lvl_obj.getSchedulingTime()) /
                (LL_serial_runtime.elapsed_time - LL_lvl_runtime.elapsed_time);
        runtime_csv.addElementToRecord(profitable, "Profitable");
        runtime_csv.addElementToRecord(1, "PG");
        if (ANALYSE_STATISTIC) {
            runtime_csv.addElementToRecord(1, "Parallelism");
            runtime_csv.addElementToRecord(1, "MaximalDeference");
            runtime_csv.addElementToRecord(1, "Deviation");
            runtime_csv.addElementToRecord(1, "Parts");
        }
        runtime_csv.addRecord();
    }


    delete tmp;
    delete CSR_A;
    delete[] y_correct;
    delete Lower_A_CSR;
    delete Lower_A_CSC;

    delete[] y_correct_perfect;
    return 0;
}


