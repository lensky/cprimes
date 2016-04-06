#include <stdlib.h>
#include <math.h>
#include "bitarray.h"

#include "utils.h"
#include "cvector.h"
#include "primewheel.h"
#include "cprimes.h"

#include "sieve.h"

void new_sieve(sieve* sv, size_t nbits, size_t start_ix) {
    sv->sievebits = (bitarray*) malloc(sizeof(bitarray));
    new_bitarray_bits(sv->sievebits, nbits);
    reset_sieve(sv);
    sv->start_ix = start_ix;
}

void free_sieve(sieve* sv) {
    free_bitarray(sv->sievebits);
    free(sv->sievebits);
}

void reset_sieve(sieve* sv) {
    fill_bitarray(sv->sievebits, true);
}

void new_sieving_prime(const primewheel * const pw,
                       sieving_prime* svp,
                       prime p,
                       size_t ixp) {
    svp->prime = p;
    svp->n = ixp / pw->spokes;
    svp->offset = ixp % pw->spokes;
    svp->delta_p = p * pw->spokes;

    svp->ixp2 = svp->n * (p + pw->cop_deltas[svp->offset]) * pw->spokes
        + pw->dix_ixpart[svp->offset * 2 * pw->spokes];
}

void new_siever(siever* siever,
                const primewheel* const pw,
                prime sieve_limit) {
    sieving_primes_to_n(pw,
                        floor(sqrt(sieve_limit)),
                        &siever->sieving_primes,
                        &siever->n_sieving_primes);
    siever->pw = pw;
    siever->sieve_limit = sieve_limit;
}

void free_siever(siever* siever) {
    for (size_t i = 0; i < siever->n_sieving_primes; i++) {
        free(siever->sieving_primes[i]);
    }
    free_cvector_data((cvector_data_unit*) siever->sieving_primes);
}

void extend_siever(siever* siever,
                   sieve* segment,
                   prime new_sieve_limit) {
    prime sqlm = (prime) floor(sqrt((double)new_sieve_limit));

    prime* newps;
    size_t n_newps;

    segmented_list_primes(siever,
                          segment,
                          siever->sieving_primes[siever->n_sieving_primes - 1]->prime + 2,
                          sqlm,
                          &newps,
                          &n_newps);

    size_t old_n_svps = siever->n_sieving_primes;
    siever->n_sieving_primes += n_newps;

    siever->sieving_primes = (sieving_prime**)
        realloc(siever->sieving_primes,
                sizeof(sieving_prime*) * siever->n_sieving_primes);

    sieving_prime** n_svp_start = siever->sieving_primes + old_n_svps;

    const primewheel* pw = siever->pw;

    for (size_t i = 0; i < n_newps; i++) {
        sieving_prime* svp = (sieving_prime*) malloc(sizeof(sieving_prime));
        prime p = newps[i];
        new_sieving_prime(pw, svp, p, wcoprimetoix(pw, p));

        n_svp_start[i] = svp;
    }

    siever->sieve_limit = new_sieve_limit;
}

size_t sieve_starting_index(const sieve* const sieve,
                            const sieving_prime* const svp) {
    if (svp->ixp2 >= sieve->start_ix) {
        return svp->ixp2 - sieve->start_ix;
    }
    size_t mm = (sieve->start_ix - svp->ixp2) % svp->delta_p;
    if (mm == 0) {
        return 0;
    }
    return svp->delta_p - mm;
}

static size_t mark_prime_above(const primewheel* const pw,
                               sieve* const sieve,
                               const sieving_prime* const svp,
                               size_t starting_ix,
                               size_t* delta_ixps) {
    size_t L_w = pw->spokes;
    bitarray* bits = sieve->sievebits;
    size_t nbits = bits->nbits;
    size_t m = starting_ix + svp->delta_p;
    size_t mdixmax = m + delta_ixps[L_w - 1];
    while (mdixmax < nbits) {
        for (size_t i = 0; i < L_w; i++) {
            clearbit(bits, m + delta_ixps[i]);
        }
        m += svp->delta_p;
        mdixmax += svp->delta_p;
    }

    if (m < nbits) {
        clearbit(bits, m);
        for (size_t i = 0; i < L_w; i++) {
            size_t mdix = m + delta_ixps[i];
            if (mdix < nbits) {
                clearbit(bits, mdix);
            } else {
                break;
            }
        }
        m += svp->delta_p;
    }

    return m;
}

size_t mark_prime_ixp2(const primewheel* const pw,
                       sieve* const sieve,
                       const sieving_prime* const svp,
                       size_t starting_ix,
                       size_t* delta_ixps) {
    if (starting_ix >= sieve->sievebits->nbits) {
        return starting_ix;
    }

    delta_ixps[0] = 0;
    clearbit(sieve->sievebits, starting_ix);
    size_t L_w = pw->spokes;
    for (size_t i = 1; i < L_w; i++) {
        size_t j = 2 * svp->offset * L_w + 2 * i;
        size_t dp = pw->dix_ixpart[j + 1];
        delta_ixps[i] = svp->n * dp * L_w + pw->dix_ixpart[j];
        size_t mdix = starting_ix + delta_ixps[i];
        if (mdix < sieve->sievebits->nbits) {
            clearbit(sieve->sievebits, mdix);
        }
    }

    return mark_prime_above(pw, sieve, svp, starting_ix, delta_ixps);
}

size_t mark_prime(const primewheel* const pw,
                  sieve* const sieve,
                  const sieving_prime* const svp,
                  size_t starting_ix,
                  size_t* delta_ixps) {
    if (sieve->start_ix <= svp->ixp2) {
        return mark_prime_ixp2(pw, sieve, svp, starting_ix, delta_ixps);
    }

    uint_fast32_t m = starting_ix;

    uint_fast32_t L_w = pw->spokes;

    if (m < sieve->sievebits->nbits) {
        delta_ixps[0] = 0;
        clearbit(sieve->sievebits, m);
        for (uint_fast32_t i = 1; i < L_w; i++) {
            uint_fast32_t j = 2 * svp->offset * L_w + 2 * i;
            uint_fast32_t dp = pw->dix_ixpart[j + 1];
            delta_ixps[i] = svp->n * dp * L_w + pw->dix_ixpart[j];
            uint_fast32_t mdix = m + delta_ixps[i];
            if (mdix >= svp->delta_p) {
                size_t mmdix = mdix - svp->delta_p;
                if (mmdix < sieve->sievebits->nbits) {
                    clearbit(sieve->sievebits, mmdix);
                }
            }
            if (mdix < sieve->sievebits->nbits) {
                clearbit(sieve->sievebits, mdix);
            }
        }
    } else {
        for (uint_fast32_t i = 1; i < L_w; i++) {
            uint_fast32_t j = 2 * svp->offset * L_w + 2 * i;
            uint_fast32_t dp = pw->dix_ixpart[j + 1];
            uint_fast32_t dixps = svp->n * dp * L_w + pw->dix_ixpart[j];
            uint_fast32_t mdix = m + dixps;
            if (mdix >= svp->delta_p) {
                uint_fast32_t mmdix = mdix - svp->delta_p;
                if (mmdix < sieve->sievebits->nbits) {
                    clearbit(sieve->sievebits, mmdix);
                }
            }
        }
        return m;
    }

    return mark_prime_above(pw, sieve, svp, starting_ix, delta_ixps);
}

void sieving_primes_to_n(const primewheel * const pw,
                         prime n,
                         sieving_prime*** svps,
                         size_t* nsvps) {
    size_t nix = find_closest_ixl(pw, n);
    if (nix <= 1) {
        *svps = NULL;
        *nsvps = 0;
        return;
    }

    sieve sieve;
    new_sieve(&sieve, nix, 1);

    cvector pvsvps;
    new_cvector(&pvsvps, sizeof(sieving_prime*));

    prime max_sieve = (size_t) floor(sqrt((double) n));
    size_t max_sieve_ixp1 = find_closest_ixl(pw, max_sieve);
    sieving_prime** svp;
    size_t* delta_ixps = (size_t*) malloc(sizeof(size_t) * pw->spokes);
    size_t i;

    for (i = 0; i < max_sieve_ixp1; i++) {
        if (testbit(sieve.sievebits, i)) {
            prime p = ixtowcoprime(pw, i + 1);
            svp = (sieving_prime**) next_elt(&pvsvps);
            *svp = (sieving_prime*) malloc(sizeof(sieving_prime));
            new_sieving_prime(pw, *svp, p, i + 1);

            mark_prime_ixp2(pw, &sieve, *svp, (*svp)->ixp2 - 1, delta_ixps);
        }
    }
    free(delta_ixps);

    for (; i < nix; i++) {
        if (testbit(sieve.sievebits, i)) {
            prime p = ixtowcoprime(pw, i + 1);
            svp = (sieving_prime**) next_elt(&pvsvps);
            *svp = (sieving_prime*) malloc(sizeof(sieving_prime));
            new_sieving_prime(pw, *svp, p, i + 1);
        }
    }

    *svps = (sieving_prime**)pvsvps.data;
    *nsvps = pvsvps.length;

    free_sieve(&sieve);
}

void segmented_sieve(siever* siever,
                     sieve* segment,
                     prime from,
                     prime to,
                     void (*action)(const primewheel * const pw,
                                    sieve* sieve,
                                    void* acc),
                     void* acc) {
    if (from < siever->pw->cop_deltas[1]) {
        from = siever->pw->cop_deltas[1];
    }

    if (to > siever->sieve_limit) {
        extend_siever(siever, segment, to);
    }

    size_t from_ix = find_closest_ixg(siever->pw, from);
    segment->start_ix = from_ix;
    size_t total_bits = (find_closest_ixl(siever->pw, to) + 1) - from_ix;

    size_t segment_bits = segment->sievebits->nbits;

    size_t full_segments = total_bits / segment_bits;
    size_t partial_segment_bits = total_bits % segment_bits;

    size_t* delta_ixps = (size_t*) malloc(sizeof(size_t)
                                          * siever->pw->spokes);

    const primewheel * const pw = siever->pw;

    if (full_segments > 0) {
        size_t* starting_ixs = (size_t*) malloc(sizeof(size_t)
                                                * siever->n_sieving_primes);

        reset_sieve(segment);
        for (size_t j = 0; j < siever->n_sieving_primes; j++) {
            sieving_prime* svp = siever->sieving_primes[j];;
            starting_ixs[j] = mark_prime(siever->pw,
                                        segment,
                                        svp,
                                        sieve_starting_index(segment, svp),
                                        delta_ixps) - segment_bits;
        }
        action(pw, segment, acc);
        segment->start_ix += segment_bits;

        for (size_t i = 1; i < full_segments; i++) {
            reset_sieve(segment);
            for (size_t j = 0; j < siever->n_sieving_primes; j++) {
                starting_ixs[j] = mark_prime(siever->pw,
                                            segment,
                                            siever->sieving_primes[j],
                                            starting_ixs[j],
                                            delta_ixps) - segment_bits;
            }
            action(pw, segment, acc);
            segment->start_ix += segment_bits;
        }
        init_bitarray_view(segment->sievebits, partial_segment_bits);
        reset_sieve(segment);
        for (size_t j = 0; j < siever->n_sieving_primes; j++) {
            mark_prime(siever->pw,
                       segment,
                       siever->sieving_primes[j],
                       starting_ixs[j],
                       delta_ixps);
        }
        action(pw, segment, acc);
        init_bitarray_view(segment->sievebits, segment_bits);
        free(starting_ixs);
    } else {
        init_bitarray_view(segment->sievebits, partial_segment_bits);
        reset_sieve(segment);
        for (size_t j = 0; j < siever->n_sieving_primes; j++) {
            sieving_prime* svp = siever->sieving_primes[j];
            mark_prime(siever->pw,
                       segment,
                       svp,
                       sieve_starting_index(segment, svp),
                       delta_ixps);
        }
        action(pw, segment, acc);
    }
    init_bitarray_view(segment->sievebits, segment_bits);

    free(delta_ixps);
}
