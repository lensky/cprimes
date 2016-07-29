/**
    @file
    @brief Implementation of prime wheels to ignore multiples of small primes.

    When sieving primes, it is useful to "automatically" ignore
    multiples of small primes. This generalizes only checking odd
    numbers for primality to eliminating multiples of any number of
    consecutive primes. Clearly this method could work for arbitarily
    large numbers of consecutive primes, but in practice it appears
    performance trade-offs appear early, with wheels handling more
    than 5 primes causing significant slowdowns.
 */

#ifndef PRIMEWHEEL_H
#define PRIMEWHEEL_H

#include <stddef.h>
#include <stdint.h>

#include "cprimes/types.h"

/// Structure for use as a generic prime wheel.
typedef struct primewheel {
    /// L_w, number of numbers coprime to primes below C_w.
    /// Note \f$ L_w = \prod_p (p - 1) \f$.
    size_t spokes;

    /// L_w x (L_w - 1) array used to calculate offsets only dependent on `p % C_w`.
    /**
       Access this array by `incr_ix_mod[(L_w - 1) * ip + n]` where `ip` is
       the index of the prime modulo L_w and `n` is the "needle"
       value.
     */
    size_t* incr_ix_mod;

    /// L_w array for sum of `incr_ix_mod[(L_w - 1) * ip + n]` over n.
    size_t* incr_ix_mod_tot;

    /// L_w - 1 array used to calculate offsets only dependent on `floor(p / C_w)`.
    size_t* incr_ix_div;

    /// C_w, product of primes
    prime circumference;

    /// dimension L_w, numbers coprime to primes < C_w
    prime* wp_coprimes;

    /// dimension C_w, `wp_coprimes[ix_of_coprime[n]] == n`
    size_t* ix_of_coprime;

    /// Primes that are to be ignored in indexing
    prime* primes;
    size_t nprimes;
} primewheel;

void new_primewheel(primewheel* pw, size_t nprimes);
void free_primewheel(primewheel* pw);

/// Takes an array index to a number according to the wheel.
static inline prime w_coprime_of_ix(const primewheel *const pw, size_t ix) {
    return (ix / pw->spokes) * pw->circumference + pw->wp_coprimes[ix % pw->spokes];
}

/// Takes a number coprime to the wheel primes to an array index.
static inline size_t w_ix_of_coprime(const primewheel *const pw, prime cp) {
    return (cp / pw->circumference) * pw->spokes
        + pw->ix_of_coprime[cp % pw->circumference];
}

/// Finds the index corresponding to the greatest number less than or
/// equal to n.
size_t closest_w_ixle(const primewheel *const pw, prime n);

/// Finds the index corresponding to the least number greater than or
/// equal to n.
size_t closest_w_ixge(const primewheel *const pw, prime n);

#endif /* PRIMEWHEEL_H */
