cmake_minimum_required( VERSION 3.0)
project(lbc)

include_directories(includes
        ${SPARSE_UTIL_INC}
        )

add_library (lbc STATIC
        lbc_csc.cpp
        lbc_utils.cpp
        lbc_tree.cpp
        lbc_csc_dag.cpp
        ${SPARSE_UTIL_SRC}/sparse_utilities.cpp
        ${SPARSE_UTIL_SRC}/etree.cpp
        ${SPARSE_UTIL_SRC}/utils.cpp ../hdagg/includes/hdagg.h)