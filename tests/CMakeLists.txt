include(ExternalProject)
ExternalProject_Add(googletest
  PREFIX "${CMAKE_CURRENT_BINARY_DIR}/googletest"
  URL https://github.com/google/googletest/archive/release-1.7.0.zip
  URL_HASH SHA256=b58cb7547a28b2c718d1e38aee18a3659c9e3ff52440297e965f5edffe34b6d0
  CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
             -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY=ReleaseLibs
             -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG=DebugLibs
             -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE=ReleaseLibs
             -Dgtest_force_shared_crt=ON
  INSTALL_COMMAND ""
  )

ExternalProject_Get_Property(googletest binary_dir source_dir)

add_executable(test-c-primes test-primes.cpp)

target_link_libraries(test-c-primes optimized "${binary_dir}/ReleaseLibs/${CMAKE_STATIC_LIBRARY_PREFIX}gtest${CMAKE_STATIC_LIBRARY_SUFFIX}")
target_link_libraries(test-c-primes debug "${binary_dir}/DebugLibs/${CMAKE_STATIC_LIBRARY_PREFIX}gtest${CMAKE_STATIC_LIBRARY_SUFFIX}")
target_include_directories(test-c-primes PUBLIC "${source_dir}/include")

target_link_libraries(test-c-primes primes)

add_test(test-c-primes test-c-primes)
