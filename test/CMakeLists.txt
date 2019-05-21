cmake_minimum_required(VERSION 3.1)

### Preparing Google Test
find_package(GTest QUIET)
if(GTEST_FOUND)
  add_executable(gtest ALIAS GTEST::GTEST)
  add_executable(gtest_main ALIAS GTEST::MAIN)
else()
  message(STATUS "GTest is not found. Installing Gtest as the external library. ")
  include(${CMAKE_CURRENT_SOURCE_DIR}/install_gtest.cmake)
endif()

# https://stackoverflow.com/questions/23163075/how-to-compile-an-mpi-included-c-program-using-cmake
find_package(MPI REQUIRED)
add_definitions(-DOMPI_SKIP_MPICXX)

set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

if(MSVC)
  add_compile_options(/W4 /WX)
else()
  add_compile_options(-Wall -Wextra -pedantic -Werror)
endif()

file(GLOB SRCS *.cpp)

add_executable(test_mpi ${SRCS})
target_link_libraries(test_mpi PRIVATE gtest MPI::MPI_C)
target_compile_features(test_mpi PRIVATE cxx_std_11)

add_test(NAME test_mpi COMMAND ${MPIEXEC} -np 4 --hostfile ${CMAKE_CURRENT_SOURCE_DIR}/hosts test_mpi)
