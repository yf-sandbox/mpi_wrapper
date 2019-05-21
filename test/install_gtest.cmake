find_package(Threads REQUIRED)

include(ExternalProject)
ExternalProject_Add(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG master
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gtest
  INSTALL_COMMAND ""
  )
ExternalProject_Get_Property(googletest source_dir binary_dir)
message("source_dir: ${source_dir}")
message("binary_dir: ${binary_dir}")

add_library(gtest IMPORTED STATIC GLOBAL)
add_dependencies(gtest googletest)
set_target_properties(gtest PROPERTIES
  #  "IMPORTED_LOCATION" "${binary_dir}/lib/gtest/libgtest.a"
  "IMPORTED_LOCATION" "${binary_dir}/lib/libgtest.a"
  "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
  )

add_library(gtest_main IMPORTED STATIC GLOBAL)
add_dependencies(gtest_main googletest)
set_target_properties(gtest_main PROPERTIES
  #  "IMPORTED_LOCATION" "${binary_dir}/lib/gtest/libgtest_main.a"
  "IMPORTED_LOCATION" "${binary_dir}/lib/libgtest_main.a"
  "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
  )

add_library(gmock IMPORTED STATIC GLOBAL)
add_dependencies(gmock googletest)
set_target_properties(gmock PROPERTIES
  #  "IMPORTED_LOCATION" "${binary_dir}/lib/gtest/libgmock.a"
  "IMPORTED_LOCATION" "${binary_dir}/lib/libgmock.a"
  "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
  )

add_library(gmock_main IMPORTED STATIC GLOBAL)
add_dependencies(gmock_main googletest)
set_target_properties(gmock_main PROPERTIES
  #  "IMPORTED_LOCATION" "${binary_dir}/lib/gtest/libgmock_main.a"
  "IMPORTED_LOCATION" "${binary_dir}/lib/libgmock_main.a"
  "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
  )

include_directories("${source_dir}/googletest/include"
  "${source_dir}/googlemock/include")
#target_include_directories(gtest gtest_main PUBLIC "${source_dir}/googletest/include")
#target_include_directories(gmock gmock_main PUBLIC "${source_dir}/googlemock/include")
