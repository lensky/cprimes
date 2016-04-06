#ifndef UTILS_H
#define UTILS_H

#include <stdint.h>
#include <stdbool.h>

#include "types.h"

// Assumes ps is ordered list of primes starting from 2.
size_t find_next_coprime(const prime *const ps, size_t nps, prime p);
bool is_coprime(const prime *const ps, size_t nps, prime p);


#endif /* UTILS_H */
