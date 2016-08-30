#include "bitarray.h"

#include "cprimes/types.h"
#include "cprimes/cvector.h"
#include "cprimes/sieve.h"

#include "cprimes/cprimes.h"

static inline void count_primes_in_sieve(const primewheel* const pw,
                                         sieve* sieve,
                                         void *acc) {
    size_t* count = (size_t*) acc;
    (*count) += count_true(sieve->sievebits);
}

size_t segmented_count_primes(siever* siever,
                              sieve* segment,
                              prime from,
                              prime to) {
    size_t count = 0;
    if (from > 2) {
        for (int i = siever->pw->nprimes - 1; i >= 0; i--) {
            if (from > siever->pw->primes[i]) {
                break;
            }
            count += 1;
        }
    } else {
        count = siever->pw->nprimes;
    }

    if (to <= siever->pw->primes[siever->pw->nprimes - 1]) {
        for (int i = siever->pw->nprimes - 1; i >= 0; i--) {
            if (to >= siever->pw->primes[i]) {
                break;
            }
            count -= 1;
        }
        return count;
    }

    if ((siever->n_sieving_primes > 0) &&
        (from <= siever->sieving_primes[siever->n_sieving_primes - 1]->prime)) {
        size_t svp_ix;
        prime p;

        for (svp_ix = 0;
             svp_ix < siever->n_sieving_primes;
             svp_ix++) {
            if (siever->sieving_primes[svp_ix]->prime >= from)
                break;
        }

        for (;
             svp_ix < siever->n_sieving_primes;
             svp_ix++) {
            p = siever->sieving_primes[svp_ix]->prime;
            if (p <= to) {
                count += 1;
            } else {
                return count;
            }
        }
        from = p + 2;
        if (from > to) {
            return count;
        }
    }

    segmented_sieve(siever,
                    segment,
                    from,
                    to,
                    &count_primes_in_sieve,
                    &count);
    return count;
}

static inline void push_prime_in_sieve(const primewheel * const pw,
                                       sieve* sieve,
                                       void *acc) {
    cvector *primesv = (cvector*)acc;
    for (size_t i = 0; i < sieve->sievebits->nbits; i++) {
        if (testbit(sieve->sievebits, i)) {
            prime* pp = (prime*) next_elt(primesv);
            *pp = w_coprime_of_ix(pw, i + sieve->ix_offset);
        }
    }
}

void segmented_list_primes(siever* siever,
                           sieve* segment,
                           prime from,
                           prime to,
                           prime** primes,
                           size_t* nprimes) {
    if (to <= 1) {
        primes = NULL;
        nprimes = 0;
        return;
    }

    cvector primesv;
    new_cvector(&primesv, sizeof(prime));

    size_t pwstart = 0;
    for (; pwstart < siever->pw->nprimes; pwstart++) {
        if (from <= siever->pw->primes[pwstart]) {
            break;
        }
    }
    for (; pwstart < siever->pw->nprimes; pwstart++) {
        prime p = siever->pw->primes[pwstart];
        if (to < p) {
            *primes = (prime*)primesv.data;
            *nprimes = primesv.length;
            return;
        }
        prime* pp = (prime*)next_elt(&primesv); *pp = p;
    }

    size_t svstart = 0;
    for (; svstart < siever->n_sieving_primes; svstart++) {
        if (from <= siever->sieving_primes[svstart]->prime) {
            break;
        }
    }
    prime p = 0;
    for (; svstart < siever->n_sieving_primes; svstart++) {
        p = siever->sieving_primes[svstart]->prime;
        if (to < p) {
            *primes = (prime*)primesv.data;
            *nprimes = primesv.length;
            return;
        }
        prime* pp = (prime*)next_elt(&primesv); *pp = p;
    }

    if (from < p + 2) {
        from = p + 2;
    }

    if (from > to) {
        *primes = (prime*) primesv.data;
        *nprimes = primesv.length;
        return;
    }

    segmented_sieve(siever,
                    segment,
                    from,
                    to,
                    &push_prime_in_sieve,
                    &primesv);

    *primes = (prime*) primesv.data;
    *nprimes = primesv.length;
    return;
}
