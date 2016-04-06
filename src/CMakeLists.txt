file(GLOB PRIMES_CFILES "*.c")
add_library(primes ${PRIMES_CFILES})
target_include_directories(primes PUBLIC ${BITARRAY_INCLUDE_DIR})
target_link_libraries(primes ${LIBBITARRAY})
