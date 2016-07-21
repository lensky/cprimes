#include "cprimes/types.h"
#include "cprimes/sieve.h"

#define DEF_BUFSIZE (256 * 1024 * 8)

size_t segmented_count_primes(siever* siever,
                              sieve* segment,
                              prime from,
                              prime to);

void segmented_list_primes(siever* siever,
                           sieve* segment,
                           prime from,
                           prime to,
                           prime** primes,
                           size_t* nprimes);
