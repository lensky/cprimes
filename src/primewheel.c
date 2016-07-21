#include <stdlib.h>

#include "cprimes/utils.h"
#include "cprimes/primewheel.h"

void new_primewheel(primewheel* pw, size_t nprimes) {
    pw->nprimes = nprimes;

    pw->primes = (prime*) malloc(sizeof(prime) * nprimes);

    pw->primes[0] = 2;
    pw->circumference = 2;

    prime p = 3;

    size_t i;

    for (i = 1; i < nprimes; i++) {
        pw->primes[i] = p;
        pw->circumference *= p;
        p = find_next_coprime(pw->primes, i, p + 2);
    }

    /* TODO: Bug for nprimes <= 1 */
    pw->spokes = 2;
    pw->cop_deltas = (prime*) malloc(sizeof(prime) * (pw->circumference / 2));
    pw->ix_deltas = (size_t*) malloc(sizeof(prime) * pw->circumference);

    pw->cop_deltas[0] = 1;
    pw->cop_deltas[1] = p;

    pw->ix_deltas[1] = 0;
    pw->ix_deltas[p] = 1;

    prime q;
    for (q = p + 2; q < pw->circumference; q += 2) {
        if (is_coprime(pw->primes, nprimes, q)) {
            pw->cop_deltas[pw->spokes] = q;
            pw->ix_deltas[q] = pw->spokes;
            pw->spokes++;
        }
    }
    pw->cop_deltas = (prime*) realloc(pw->cop_deltas, sizeof(prime) * pw->spokes);

    prime dixi0;
    prime dp;
    pw->dix_ixpart = (size_t*) malloc(pw->spokes * pw->spokes * 2 * sizeof(size_t));
    for (i = 0; i < pw->spokes; i++) {
        q = pw->cop_deltas[i];
        dixi0 = wcoprimetoix(pw, q * q);
        pw->dix_ixpart[i * 2 * pw->spokes] = dixi0;
        pw->dix_ixpart[i * 2 * pw->spokes + 1] = 0;
        for (size_t j = 1; j < pw->spokes; j++) {
            dp = pw->cop_deltas[(i + j) % pw->spokes];
            if (dp < q) {
                dp += pw->circumference;
            }
            pw->dix_ixpart[i * 2 * pw->spokes + 2 * j] =
                wcoprimetoix(pw, q * dp) - dixi0;
            pw->dix_ixpart[i * 2 * pw->spokes + 2 * j + 1] =
                dp - q;
        }
    }
}

void free_primewheel(primewheel* pw) {
    free(pw->primes);
    free(pw->cop_deltas);
    free(pw->ix_deltas);
    free(pw->dix_ixpart);
}

prime ixtowcoprime(const primewheel *const pw, size_t ix) {
    return (ix / pw->spokes) * pw->circumference + pw->cop_deltas[ix % pw->spokes];
}

size_t wcoprimetoix(const primewheel *const pw, prime cp) {
    return (cp / pw->circumference) * pw->spokes
        + pw->ix_deltas[cp % pw->circumference];
}

size_t find_closest_ixl(const primewheel *const pw, prime n) {
    if (!(n & 1)) {
        n -= 1;
    }

    prime ndcw = n / pw->circumference;
    prime nrcw = n % pw->circumference;

    size_t i;
    for (i = 0; (i < pw->spokes) && (pw->cop_deltas[i] <= nrcw); i++) {}
    return (ndcw * pw->spokes + i) - 1;
}

size_t find_closest_ixg(const primewheel *const pw, prime n) {
    n |= 1;
    prime ndcw = n / pw->circumference;
    prime nrcw = n % pw->circumference;

    size_t i;
    for (i = 0; (i < pw->spokes) && (pw->cop_deltas[i] < nrcw); i++) {}
    return ndcw * pw->spokes + i;
}
