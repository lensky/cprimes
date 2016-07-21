#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "cprimes/types.h"
#include "cprimes/utils.h"

bool is_coprime(const prime *const ps, size_t nps, prime p) {
    size_t fsqrtp = (size_t) floor(sqrt((double) p));
    prime q;
    for (size_t i = 0; i < nps; i++) {
        q = ps[i];
        if (q > fsqrtp) {
            return true;
        }
        if (p % q == 0) {
            return false;
        }
    }
    return true;
}

size_t find_next_coprime(const prime *const ps, size_t nps, prime p) {
    p |= 1;
    while (1) {
        if (is_coprime(ps, nps, p)) {
            return p;
        }
        p += 2;
    }
}
