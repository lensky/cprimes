cmake_minimum_required(VERSION 3.3)
project(cprimes)

option(BUILD_SHARED_LIBS "Build as a shared library." OFF)
option(PROFILE_GPROF "Enable profiling with gprof." OFF)

if(CMAKE_COMPILER_IS_GNUCC)
  set(CMAKE_C_FLAGS
    "${CMAKE_C_FLAGS} -mpopcnt")
  set(CMAKE_C_FLAGS_RELEASE
    "${CMAKE_C_FLAGS_RELEASE} -O3 -mtune=native -march=native -funroll-loops -ffast-math")
  if(PROFILE_GPROF)
    set(CMAKE_C_FLAGS
      "${CMAKE_C_FLAGS} -pg")
    set(CMAKE_EXE_LINKER_FLAGS
      "${CMAKE_EXE_LINKER_FLAGS} -pg")
  endif()
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -mpopcnt")
  set(CMAKE_CXX_FLAGS_RELEASE
    "${CMAKE_CXX_FLAGS_RELEASE} -O3 -mtune=native -march=native -funroll-loops -ffast-math")
endif()

set(CMAKE_C_STANDARD 11)

find_library(
  LIBBITARRAY
  NAMES bitarray libbitarray
  )
find_path(
  BITARRAY_INCLUDE_DIR
  bitarray.h
  )

include_directories(include)
install(DIRECTORY include DESTINATION include)
add_subdirectory(src)

enable_testing()

add_subdirectory(tests)
add_subdirectory(bench)
