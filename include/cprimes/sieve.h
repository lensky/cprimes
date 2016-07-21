#ifndef SIEVE_H
#define SIEVE_H

#include "bitarray.h"

#include "cprimes/types.h"
#include "cprimes/primewheel.h"

typedef struct sieve {
    bitarray* sievebits;
    size_t start_ix;
} sieve;

typedef struct sieving_prime {
    size_t delta_p; // p * L_w
    size_t ixp2; // index of p^2
    size_t offset; // ix_deltas[p % C_w]
    prime n; // p / C_w
    prime prime;
} sieving_prime;

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

void extend_siever(siever* siever,
                   sieve* segment,
                   prime new_sieve_limit);

void new_sieve(sieve* sv, size_t nbits, size_t start_ix);
void free_sieve(sieve* sv);
void reset_sieve(sieve* sv);

void new_sieving_prime(const primewheel* const pw,
                       sieving_prime* svp,
                       prime p,
                       size_t ixp);

size_t sieve_starting_index(const sieve* const sieve,
                            const sieving_prime* const svp);

size_t mark_prime(const primewheel *const pw,
                  sieve* const sieve,
                  const sieving_prime* const svp,
                  size_t starting_ix,
                  size_t* delta_ixps // workspace array of dimension L_w
    );

void sieving_primes_to_n(const primewheel *const pw,
                         prime n,
                         sieving_prime*** svps,
                         size_t* nsvps);

void segmented_sieve(siever* siever,
                     sieve* segment,
                     prime from,
                     prime to,
                     void (*action)(const primewheel * const pw,
                                    sieve* sieve,
                                    void* acc),
                     void* acc);

#endif /* SIEVE_H */
