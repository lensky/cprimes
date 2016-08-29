find_library(
  LIBCBENCH
  NAMES cbench libcbench
  )

find_path(
  CBENCH_INCLUDES
  cbench.h
  )

add_executable(bench-primes bench-primes.c)

target_link_libraries(bench-primes ${LIBCBENCH})
target_include_directories(bench-primes PUBLIC ${CBENCH_INCLUDES})

if (UNIX)
  target_link_libraries(bench-primes m)
endif (UNIX)

target_link_libraries(bench-primes primes)
