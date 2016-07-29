#include <stdlib.h>

#include "cprimes/utils.h"
#include "cprimes/primewheel.h"

void new_primewheel(primewheel* pw, size_t nprimes) {
    pw->nprimes = nprimes;

    pw->primes = (prime*) malloc(sizeof(prime) * nprimes);

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

    pw->wp_coprimes = (prime*) malloc(sizeof(prime) * pw->spokes);
    pw->ix_of_coprime = (size_t*) malloc(sizeof(prime) * pw->circumference);

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

    pw->incr_ix_mod = (size_t*) malloc(pw->spokes * (pw->spokes - 1) * sizeof(size_t));
    pw->incr_ix_mod_tot = (size_t*) malloc(pw->spokes * sizeof(size_t));
    pw->incr_ix_div = (size_t*) malloc((pw->spokes - 1) * sizeof(size_t));

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
