#ifndef PRIMEWHEEL_H
#define PRIMEWHEEL_H

#include <stddef.h>
#include <stdint.h>

#include "cprimes/types.h"

typedef struct primewheel {
    prime* primes;
    size_t nprimes;

    prime circumference; // C_w, product of primes

    prime* cop_deltas; // dimension L_w, numbers coprime to primes
    size_t spokes; // L_w, number of numbers coprime to primes below
                   // C_w
    size_t* ix_deltas; // dimension C_w, cop_deltas[ix_deltas[n]] = n

    size_t* dix_ixpart;
} primewheel;

void new_primewheel(primewheel* pw, size_t nprimes);
void free_primewheel(primewheel* pw);

prime ixtowcoprime(const primewheel *const pw, size_t ix);
size_t wcoprimetoix(const primewheel *const pw, prime cp);

size_t find_closest_ixl(const primewheel *const pw, prime n);
size_t find_closest_ixg(const primewheel *const pw, prime n);

#endif /* PRIMEWHEEL_H */
