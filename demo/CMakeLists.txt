set(REQUIRED_CMAKE_VERSION "3.9")

################################################################################
# Testing SpTRSV kernel
add_executable(${PROJECT_NAME}_SpTRSV SpTRSV_runtime.cpp)
target_link_libraries(${PROJECT_NAME}_SpTRSV PUBLIC ${PROJECT_NAME}_lib)
# CLI11
#target_link_libraries(${PROJECT_NAME}_SpTRSV PUBLIC CLI11::CLI11)
################################################################################

