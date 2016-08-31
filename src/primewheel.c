#include <stdlib.h>
#include <string.h>

#include "cprimes/utils.h"
#include "cprimes/primewheel.h"

#define PW_ALLOC_PRIMES(pw) (pw->primes = malloc((sizeof(*(pw->primes)) * pw->nprimes)))

#define PW_DIM_IX_MOD(pw) (pw->spokes * (pw->spokes - 1))
#define PW_DIM_IX_MOD_TOT(pw) (pw->spokes)
#define PW_DIM_IX_DIV(pw) (pw->spokes - 1)

// Does all allocation except pw->primes
static inline void alloc_primewheel(primewheel* pw) {
    pw->wp_coprimes = (prime*) malloc(sizeof(prime) * pw->spokes);
    pw->ix_of_coprime = (size_t*) malloc(sizeof(prime) * pw->circumference);

    pw->incr_ix_mod = (size_t*) malloc(PW_DIM_IX_MOD(pw) * sizeof(size_t));
    pw->incr_ix_mod_tot = (size_t*) malloc(PW_DIM_IX_MOD_TOT(pw) * sizeof(size_t));
    pw->incr_ix_div = (size_t*) malloc(PW_DIM_IX_DIV(pw) * sizeof(size_t));

}

void new_primewheel(primewheel* pw, size_t nprimes) {
    pw->nprimes = nprimes;

    PW_ALLOC_PRIMES(pw);

    pw->primes[0] = 2;
    pw->circumference = 2;
    pw->spokes = 1;

    prime p = 3;

    size_t i;

    for (i = 1; i < nprimes; i++) {
        pw->primes[i] = p;
        pw->circumference *= p;
        pw->spokes *= (p - 1);
        p = find_next_coprime(pw->primes, i, p + 2);
    }

    alloc_primewheel(pw);

    pw->wp_coprimes[0] = 1;
    pw->wp_coprimes[1] = p;

    pw->ix_of_coprime[1] = 0;
    pw->ix_of_coprime[p] = 1;

    prime q;
    i = 2;
    for (q = p + 2; q < pw->circumference; q += 2) {
        if (is_coprime(pw->primes, nprimes, q)) {
            pw->wp_coprimes[i] = q;
            pw->ix_of_coprime[q] = i++; // Note we increment i here
        }
    }

    for (i = 0; i < pw->spokes; i++) {
        prime cp = pw->wp_coprimes[i];
        pw->incr_ix_mod_tot[i] = 0;
        for (size_t n = 0; n < pw->spokes - 1; n++) {
            prime cn = pw->wp_coprimes[n];
            prime cn1 = pw->wp_coprimes[n + 1];
            prime cpcn = cp * cn;
            prime cpcn1 = cp * cn1;
            pw->incr_ix_mod[(pw->spokes - 1) * i + n] =
                (pw->ix_of_coprime[cpcn1 % pw->circumference]
                 + pw->spokes * ((cpcn1 / pw->circumference)
                                 - (cpcn / pw->circumference)))
                - pw->ix_of_coprime[cpcn % pw->circumference];
            pw->incr_ix_mod_tot[i] += pw->incr_ix_mod[(pw->spokes - 1) * i + n];
        }

        if (i < pw->spokes - 1) {
            pw->incr_ix_div[i] = pw->spokes * (pw->wp_coprimes[i + 1] - cp);
        }
    }
}

void copy_primewheel(primewheel* pw_dest, primewheel* pw_source) {
    *pw_dest = *pw_source;

    PW_ALLOC_PRIMES(pw_dest);
    alloc_primewheel(pw_dest);

    memcpy(pw_dest->primes, pw_source->primes,
           sizeof(*(pw_dest->primes)) * pw_dest->nprimes);
    memcpy(pw_dest->wp_coprimes, pw_source->wp_coprimes,
           sizeof(*(pw_dest->wp_coprimes)) * pw_dest->spokes);
    memcpy(pw_dest->ix_of_coprime, pw_source->ix_of_coprime,
           sizeof(*(pw_dest->ix_of_coprime)) * pw_dest->circumference);

    memcpy(pw_dest->incr_ix_mod, pw_source->incr_ix_mod,
           sizeof(*(pw_dest->incr_ix_mod)) * PW_DIM_IX_MOD(pw_dest));
    memcpy(pw_dest->incr_ix_mod_tot, pw_source->incr_ix_mod_tot,
           sizeof(*(pw_dest->incr_ix_mod_tot)) * PW_DIM_IX_MOD_TOT(pw_dest));
    memcpy(pw_dest->incr_ix_div, pw_source->incr_ix_div,
           sizeof(*(pw_dest->incr_ix_div)) * PW_DIM_IX_DIV(pw_dest));
}

void free_primewheel(primewheel* pw) {
    free(pw->primes);
    free(pw->wp_coprimes);
    free(pw->ix_of_coprime);

    free(pw->incr_ix_mod);
    free(pw->incr_ix_mod_tot);
    free(pw->incr_ix_div);
}

size_t closest_w_ixle(const primewheel *const pw, prime n) {
    if (!(n & 1)) {
        n -= 1;
    }

    prime ndcw = n / pw->circumference;
    prime nrcw = n % pw->circumference;

    size_t i;
    for (i = 0; (i < pw->spokes) && (pw->wp_coprimes[i] <= nrcw); i++) {}
    return (ndcw * pw->spokes + i) - 1;
}

size_t closest_w_ixge(const primewheel *const pw, prime n) {
    n |= 1;
    prime ndcw = n / pw->circumference;
    prime nrcw = n % pw->circumference;

    size_t i;
    for (i = 0; (i < pw->spokes) && (pw->wp_coprimes[i] < nrcw); i++) {}
    return ndcw * pw->spokes + i;
}
