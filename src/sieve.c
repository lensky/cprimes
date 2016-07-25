#include <stdlib.h>
#include <math.h>
#include "bitarray.h"

#include "cprimes/utils.h"
#include "cprimes/cvector.h"
#include "cprimes/primewheel.h"
#include "cprimes/cprimes.h"

#include "cprimes/sieve.h"

void new_sieve(sieve* sv, size_t nbits, size_t ix_offset) {
    sv->sievebits = (bitarray*) malloc(sizeof(bitarray));
    new_bitarray_bits(sv->sievebits, nbits);
    reset_sieve(sv);
    sv->ix_offset = ix_offset;
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
    svp->p_div_C = ixp / pw->spokes;
    svp->offset = ixp % pw->spokes;

    svp->needle = svp->offset;
    svp->next_ix = wcoprimetoix(pw, p*p);
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
                          siever->sieving_primes[siever->n_sieving_primes - 1]->prime
                          + 2,
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

/// Calculates the increment to the index for the sieving prime given a needle.
// TODO: Look at speeding this function up
inline size_t calc_svp_ix_increment(const primewheel* const pw,
                                    const sieving_prime* svp,
                                    uint_fast32_t needle) {
    if (needle < (pw->spokes - 1)) {
        return pw->incr_ix_mod[(pw->spokes - 1) * svp->offset
                               + needle]
            + svp->p_div_C * pw->incr_ix_div[needle];
    } else {
        return pw->spokes * (svp->prime - svp->p_div_C
                             * (pw->wp_coprimes[pw->spokes - 1]
                                 - 1))
            - pw->incr_ix_mod_tot[svp->offset];
    }
}

/// Marks the prime in the current sieve and prepares for sieving with
/// a consecutive buffer.
/**
   This is the most called function in the library. It marks
   occurences of multiples of the prime in `svp` in sieve, then
   prepares `svp` for marking in a consecutive sieve.

   @param[in] ix_workspace auxilliary space of length pw->spokes
   @param[in] pw primewheel
   @param[in] sieve current sieve
   @param[in,out] svp sieving prime, will have next_ix appropriately updated
 */
inline void do_mark_prime(size_t* ix_workspace,
                          const primewheel* const pw,
                          sieve* const sieve,
                          sieving_prime* svp) {
    bitarray* bits = sieve->sievebits;
    // Useful for large primes. We lose nothing here as this check
    // must be performed anyway.
    if (svp->next_ix >= bits->nbits) {
        svp->next_ix -= bits->nbits;
        return;
    }

    size_t dp = pw->spokes * svp->prime;
    size_t tmp_ix = svp->next_ix;
    size_t old_needle = svp->needle;

    CLEARBIT(bits->bits, tmp_ix);
    ix_workspace[0] = tmp_ix + dp;

    uint_fast32_t i;
    for (i = 1; i < pw->spokes; i++) {
        tmp_ix += calc_svp_ix_increment(pw, svp, svp->needle);
        svp->needle = svp->needle < pw->spokes - 1 ? svp->needle + 1 : 0;

        if (tmp_ix < bits->nbits) {
            CLEARBIT(bits->bits, tmp_ix);
            ix_workspace[i] = tmp_ix + dp;
        } else {
            svp->next_ix = tmp_ix - bits->nbits;
            return;
        }
    }
    svp->needle = old_needle;

    // Note: this separate loop without checks in the inner portion is
    // FAR more efficient than the alternative.
    while (ix_workspace[pw->spokes - 1] < bits->nbits) {
        for (i = 0; i < pw->spokes; i++) {
            CLEARBIT(bits->bits, ix_workspace[i]);
            ix_workspace[i] += dp;
        }
    }

    for (i = 0; i < pw->spokes; i++) {
        if (ix_workspace[i] < bits->nbits) {
            CLEARBIT(bits->bits, ix_workspace[i]);
        } else {
            svp->next_ix = ix_workspace[i] - bits->nbits;
            svp->needle += i;
            if (svp->needle >= pw->spokes) { svp->needle -= pw->spokes; }
            return;
        }
    }
}

/// Initializes sieving prime for sieving, starting with the given sieve.
// At the moment, this function calculates the initial index and needle position.
void init_svp_for_sieving(const primewheel* const pw,
                          const sieve* const sieve,
                          sieving_prime* svp) {
    svp->needle = svp->offset;

    size_t ixp2 = wcoprimetoix(pw, svp->prime * svp->prime);
    size_t delta_p = pw->spokes * svp->prime;

    if (ixp2 < sieve->ix_offset) {
        svp->next_ix = delta_p
            - ((sieve->ix_offset - ixp2) % delta_p);

        size_t tmp_needle = svp->needle > 0 ? svp->needle - 1 : pw->spokes - 1;
        size_t incr = calc_svp_ix_increment(pw, svp, tmp_needle);
        while (svp->next_ix >= incr) {
            svp->next_ix -= incr;
            svp->needle = tmp_needle;

            tmp_needle = tmp_needle > 0 ? tmp_needle - 1 : pw->spokes - 1;
            incr = calc_svp_ix_increment(pw, svp, tmp_needle);
        }
    } else {
        svp->next_ix = ixp2 - sieve->ix_offset;
    }
}

void sieving_primes_to_n(const primewheel * const pw,
                         prime n,
                         sieving_prime*** svps,
                         size_t* nsvps) {
    size_t nix = closest_w_ixle(pw, n);
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
    size_t max_sieve_ixp1 = closest_w_ixle(pw, max_sieve);
    sieving_prime** svp;
    size_t i;

    size_t* ix_workspace = (size_t*) malloc(sizeof(size_t) * pw->spokes);

    for (i = 0; i < max_sieve_ixp1; i++) {
        if (testbit(sieve.sievebits, i)) {
            prime p = ixtowcoprime(pw, i + 1);
            svp = (sieving_prime**) next_elt(&pvsvps);
            *svp = (sieving_prime*) malloc(sizeof(sieving_prime));
            new_sieving_prime(pw, *svp, p, i + 1);
            (*svp)->next_ix -= sieve.ix_offset;

            do_mark_prime(ix_workspace, pw, &sieve, *svp);
        }
    }

    free(ix_workspace);

    for (; i < nix; i++) {
        if (testbit(sieve.sievebits, i)) {
            prime p = ixtowcoprime(pw, i + 1);
            svp = (sieving_prime**) next_elt(&pvsvps);
            *svp = (sieving_prime*) malloc(sizeof(sieving_prime));
            new_sieving_prime(pw, *svp, p, i + 1);
            (*svp)->next_ix -= sieve.ix_offset;
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
    if (from < siever->pw->wp_coprimes[1]) {
        from = siever->pw->wp_coprimes[1];
    }

    if (to > siever->sieve_limit) {
        extend_siever(siever, segment, to);
    }

    size_t from_ix = closest_w_ixge(siever->pw, from);
    segment->ix_offset = from_ix;
    size_t total_bits = (closest_w_ixle(siever->pw, to) + 1) - from_ix;

    size_t segment_bits = segment->sievebits->nbits;

    size_t full_segments = total_bits / segment_bits;
    size_t partial_segment_bits = total_bits % segment_bits;

    const primewheel * const pw = siever->pw;

    size_t* ix_workspace = (size_t*) malloc(sizeof(size_t) * pw->spokes);

    if (full_segments > 0) {
        reset_sieve(segment);
        for (size_t j = 0; j < siever->n_sieving_primes; j++) {
            sieving_prime* svp = siever->sieving_primes[j];;
            init_svp_for_sieving(siever->pw, segment, svp);

            do_mark_prime(ix_workspace, siever->pw, segment, svp);
        }
        action(pw, segment, acc);
        segment->ix_offset += segment_bits;

        for (size_t i = 1; i < full_segments; i++) {
            reset_sieve(segment);
            for (size_t j = 0; j < siever->n_sieving_primes; j++) {
                do_mark_prime(ix_workspace,
                              siever->pw,
                              segment,
                              siever->sieving_primes[j]);
            }
            action(pw, segment, acc);
            segment->ix_offset += segment_bits;
        }
        init_bitarray_view(segment->sievebits, partial_segment_bits);
        reset_sieve(segment);
        for (size_t j = 0; j < siever->n_sieving_primes; j++) {
            do_mark_prime(ix_workspace,
                          siever->pw,
                          segment,
                          siever->sieving_primes[j]);
        }
        action(pw, segment, acc);
        init_bitarray_view(segment->sievebits, segment_bits);
    } else {
        init_bitarray_view(segment->sievebits, partial_segment_bits);
        reset_sieve(segment);
        for (size_t j = 0; j < siever->n_sieving_primes; j++) {
            sieving_prime* svp = siever->sieving_primes[j];
            init_svp_for_sieving(siever->pw, segment, svp);
            do_mark_prime(ix_workspace, siever->pw, segment, svp);
        }
        action(pw, segment, acc);
    }
    init_bitarray_view(segment->sievebits, segment_bits);
    free(ix_workspace);
}
