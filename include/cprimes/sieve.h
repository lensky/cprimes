#ifndef SIEVE_H
#define SIEVE_H

#include <stdint.h>

#include "bitarray.h"

#include "cprimes/types.h"
#include "cprimes/primewheel.h"

typedef struct sieve {
    bitarray* sievebits;
    size_t ix_offset;
} sieve;

typedef struct sieving_prime {
    uint_fast64_t next_ix;
    uint_fast32_t needle;

    uint_fast32_t offset; ///< `ix_of_coprime[p % C_w]` (`== ix[p] %
                          ///  L_w` if `ix[1] == 0`)
    prime p_div_C; ///< `p / C_w == ix[p] / L_w`

    prime prime;
} sieving_prime;

static inline void copy_sieving_prime(sieving_prime* dest, sieving_prime* source) {
    *dest = *source;
}

typedef struct siever {
    const primewheel* pw;
    sieving_prime** sieving_primes;
    size_t n_sieving_primes;
    prime sieve_limit;
} siever;

void new_siever(siever* siever,
                const primewheel* const pw,
                prime sieve_limit);
void free_siever(siever* siever);
void copy_siever(siever* dest, siever* source);

void extend_siever(siever* siever,
                   sieve* segment,
                   prime new_sieve_limit);

void new_sieve(sieve* sv, size_t nbits, size_t ix_offset);
void free_sieve(sieve* sv);
void reset_sieve(sieve* sv);

void init_svp_for_sieving(const primewheel* const pw,
                          size_t ix_offset,
                          sieving_prime* svp);

void new_sieving_prime(const primewheel* const pw,
                       sieving_prime* svp,
                       prime p,
                       size_t ixp);

void sieving_primes_to_n(const primewheel *const pw,
                         prime n,
                         sieving_prime*** svps,
                         size_t* nsvps);

void do_mark_prime(size_t* ix_workspace,
                   const primewheel* const pw,
                   sieve* const sieve,
                   sieving_prime* svp);

/// Runs a segmented sieve from numbers `from` and `to` inclusive,
/// running `action` after every segment is processed.
void segmented_sieve(siever* siever,
                     sieve* segment,
                     prime from,
                     prime to,
                     void (*action)(const primewheel * const pw,
                                    sieve* sieve,
                                    void* acc),
                     void* acc);

#endif /* SIEVE_H */
