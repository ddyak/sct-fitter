cmake_minimum_required (VERSION 3.0)

################################
# Project setup
################################
set(PROJECT_NAME TreeFitter)
project(${PROJECT_NAME})

################################
# Compiler flags
################################
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3")

###############################
# Include directories
###############################

find_package(Eigen3 REQUIRED)

include_directories(
    ./treefitter/include/
    ./external/
    ./dataobjects/
    ./
)

###############################
# Source files
###############################

file(GLOB SOURCES 
    treefitter/src/*.cpp 
    dataobjects/*.cpp
    external/*.cpp
)

########################################
# Main -- separate executable
########################################
list(REMOVE_ITEM SOURCES main.cpp)

########################################
# Compile source files into a library
########################################
add_library(fit_lib ${SOURCES})
target_link_libraries(fit_lib Eigen3::Eigen)


########################################
# Main is separate (e.g. library client)
########################################
add_executable(main main.cpp)
target_link_libraries(main Eigen3::Eigen fit_lib)

########################################
# Testing
########################################

# google test is a git submodule for the project, and it is also cmake-based
add_subdirectory(./googletest)

enable_testing()
# 
# Include the gtest library. gtest_SOURCE_DIR is available due to
# 'project(gtest)' above.
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
 
########################################
# Test files
########################################
file(GLOB TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/tests/*.cpp)

########################################
# Unit Tests
########################################
add_executable(runUnitTests ${TEST_SRC_FILES})

########################################
# Standard linking to gtest stuff.
########################################
target_link_libraries(runUnitTests gtest gtest_main)

########################################
# Extra linking for the project.
########################################
target_link_libraries(runUnitTests Eigen3::Eigen fit_lib)

