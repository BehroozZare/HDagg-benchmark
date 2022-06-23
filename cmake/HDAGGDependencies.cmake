# Prepare dependencies
#
# For each third-party library, if the appropriate target doesn't exist yet,
# download it via external project, and add_subdirectory to build it alongside
# this project.


# Download and update 3rd_party libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)
include(HDAGGDownloadExternal)

################################################################################
# Required libraries
################################################################################

# CLI11
#if(NOT TARGET CLI11::CLI11)
#  download_cli11()
#  add_subdirectory(${HDAGG_EXTERNAL}/cli11)
#endif()

## spdlog
#if(NOT TARGET spdlog::spdlog)
#  download_spdlog()
#  add_library(spdlog INTERFACE)
#  add_library(spdlog::spdlog ALIAS spdlog)
#  target_include_directories(spdlog SYSTEM INTERFACE ${HDAGG_EXTERNAL}/spdlog/include)
#endif()

