include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
  set(HDAGG_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
  set(HDAGG_EXTRA_OPTIONS "")
endif()

function(custom_download_project name)
  download_project(
    PROJ         ${name}
    SOURCE_DIR   ${HDAGG_EXTERNAL}/${name}
    DOWNLOAD_DIR ${HDAGG_EXTERNAL}/.cache/${name}
    QUIET
    ${HDAGG_EXTRA_OPTIONS}
    ${ARGN}
  )
endfunction()

################################################################################

# CLI11
function(download_cli11)
    custom_download_project(cli11
        GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
        GIT_TAG        v2.2.0
    )
endfunction()


## Logger
#function(download_spdlog)
#    custom_download_project(spdlog
#            GIT_REPOSITORY https://github.com/gabime/spdlog.git
#            GIT_TAG        v1.9.2
#            )
#endfunction()