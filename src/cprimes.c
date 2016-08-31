#include <stdbool.h>

#include "bitarray.h"

#include "cprimes/primewheel.h"

#include "cprimes/types.h"
#include "cprimes/cvector.h"
#include "cprimes/sieve.h"

#include "cprimes/cprimes.h"

static inline prime prime_generator_starting_prime(const prime_generator* const pg) {
    prime starting;

    const siever* const svr = &(pg->siever);
    size_t old_n_svps = svr->n_sieving_primes;

    if (old_n_svps > 0) {
        starting = svr->sieving_primes[old_n_svps - 1]->prime;
    } else {
        const primewheel* const pw = svr->pw;
        starting = pw->primes[pw->nprimes - 1];
    }
    return starting + 2;
}

void new_prime_generator(prime_generator* pg,
                         const primewheel* const pw,
                         prime init_limit) {
    siever* svr = &(pg->siever);
    new_siever(svr, pw, init_limit);
    pg->listing_wheel = true;
    pg->index = 0;

    pg->starting_ix = closest_w_ixge(pw, prime_generator_starting_prime(pg));

    for (size_t i = 0; i < svr->n_sieving_primes; i++) {
        init_svp_for_sieving(pw, pg->starting_ix, svr->sieving_primes[i]);
    }
}

void free_prime_generator(prime_generator* pg) {
    free_siever(&(pg->siever));
}

void copy_prime_generator(prime_generator* dest, prime_generator* source) {
    *dest = *source;
    copy_siever(&(dest->siever), &(source->siever));
}

prime prime_generator_next(prime_generator* pg, sieve* segment) {
    if (pg->listing_wheel) {
        if (pg->index < pg->siever.pw->nprimes) {
            return pg->siever.pw->primes[pg->index++];
        } else {
            pg->listing_wheel = false;
            pg->index = 0;
            return prime_generator_next(pg, segment);
        }
    } else {
        siever* svr = &(pg->siever);
        size_t old_n_svps = svr->n_sieving_primes;

        if (pg->index < old_n_svps) {
            return svr->sieving_primes[pg->index++]->prime;
        }

        const primewheel* pw = svr->pw;
        size_t segment_bits = segment->sievebits->nbits;

        reset_sieve(segment);

        segment->ix_offset = pg->starting_ix;
        pg->starting_ix += segment_bits;
        size_t* ix_workspace = malloc(sizeof(*ix_workspace) * pw->spokes);

        for (size_t j = 0; j < old_n_svps; j++) {
            do_mark_prime(ix_workspace, pw, segment, svr->sieving_primes[j]);
        }

        prime ending = w_coprime_of_ix(pw, pg->starting_ix - 1);
        size_t k = 0;
        size_t tot_nprimes = 0;
        while (svr->sieve_limit < ending) {
            size_t ix_limit = closest_w_ixle(pw, svr->sieve_limit);
            size_t tmp_k = ix_limit - segment->ix_offset + 1;

            init_bitarray_view(segment->sievebits, tmp_k);
            size_t n_newprimes = count_true(segment->sievebits) - tot_nprimes;
            init_bitarray_view(segment->sievebits, segment_bits);

            if (n_newprimes > 0) {
                svr->sieving_primes = realloc(svr->sieving_primes,
                                              sizeof(*(svr->sieving_primes))
                                              * (svr->n_sieving_primes + n_newprimes));
                sieving_prime** new_svp_start = svr->sieving_primes + svr->n_sieving_primes;
                sieving_prime** svp_end = new_svp_start + n_newprimes;
                for (size_t j = k; new_svp_start < svp_end; j++) {
                    if (testbit(segment->sievebits, j)) {
                        sieving_prime* svp = malloc(sizeof(*svp));
                        size_t pix = j + segment->ix_offset;
                        prime p = w_coprime_of_ix(pw, pix);
                        new_sieving_prime(pw, svp, p, pix);

                        init_svp_for_sieving(pw, segment->ix_offset, svp);
                        do_mark_prime(ix_workspace, pw, segment, svp);

                        *new_svp_start = svp; new_svp_start++;
                    }
                }
                svr->n_sieving_primes += n_newprimes;
                tot_nprimes += n_newprimes;
            }
            k = tmp_k;
            svr->sieve_limit *= svr->sieve_limit;
        }
        free(ix_workspace);

        tot_nprimes = count_true(segment->sievebits) - tot_nprimes;
        svr->sieving_primes = realloc(svr->sieving_primes,
                                      sizeof(*(svr->sieving_primes))
                                      * (svr->n_sieving_primes + tot_nprimes));
        sieving_prime** new_svp_start = svr->sieving_primes + svr->n_sieving_primes;
        sieving_prime** svp_end = new_svp_start + tot_nprimes;
        svr->n_sieving_primes += tot_nprimes;
        for (; new_svp_start < svp_end; k++) {
            if (testbit(segment->sievebits, k)) {
                sieving_prime* svp = malloc(sizeof(*svp));
                size_t pix = k + segment->ix_offset;
                new_sieving_prime(pw, svp, w_coprime_of_ix(pw, pix), pix);
                init_svp_for_sieving(pw, pg->starting_ix, svp);
                *new_svp_start = svp; new_svp_start++;
            }
        }
        if (svr->sieve_limit < ending * ending) {
            svr->sieve_limit = ending * ending;
        }
        return prime_generator_next(pg, segment);
    }
}

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
