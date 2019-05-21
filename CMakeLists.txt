cmake_minimum_required(VERSION 3.1)

project(mpi_wrapper)

#message("project_source_dir: ${PROJECT_SOURCE_DIR}")
#message("cmake_current_source_dir: ${CMAKE_CURRENT_SOURCE_DIR}")

add_library(mpi_wrapper INTERFACE)
target_sources(mpi_wrapper INTERFACE ${PROJECT_SOURCE_DIR}/include/mpi_wrapper.hpp)
target_include_directories(mpi_wrapper INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

option(ENABLE_TESTING "enable testing" ON)

if(ENABLE_TESTING)
  include_directories(./include)
  enable_testing()
  add_subdirectory(test)
endif()

