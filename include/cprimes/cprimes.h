#include <stdbool.h>

#include "cprimes/types.h"
#include "cprimes/sieve.h"

#define DEF_BUFSIZE (256 * 1024 * 8)

typedef struct prime_generator {
    siever siever;

    bool listing_wheel;
    size_t index;

    size_t starting_ix;
} prime_generator;

void new_prime_generator(prime_generator* pg,
                         const primewheel* const pw,
                         prime init_limit);
void free_prime_generator(prime_generator* pg);
void copy_prime_generator(prime_generator* dest, prime_generator* source);

prime prime_generator_next(prime_generator* pg, sieve* segment);

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
