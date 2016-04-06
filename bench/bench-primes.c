#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "cbench.h"

#include "bitarray.h"
#include "primewheel.h"
#include "sieve.h"
#include "cprimes.h"

const size_t unsegmented_maxp = 1e7;

const size_t segmented_maxp = 1e9;

const size_t bufsize = DEF_BUFSIZE;
const size_t wheel_size = 4;

CBENCH_DEF(SegmentedCount) {
    primewheel pw;
    new_primewheel(&pw, wheel_size);

    siever svr;
    new_siever(&svr, &pw, segmented_maxp);

    sieve segment;
    new_sieve(&segment, bufsize, 0);

    CBENCH_INIT("Segmented prime sieve (count) to %lg", (double)segmented_maxp) {
        CBENCH_MEASURE_START
        segmented_count_primes(&svr, &segment, 0, segmented_maxp);
        CBENCH_MEASURE_STOP
    }

    free_sieve(&segment);
    free_siever(&svr);
    free_primewheel(&pw);
}

CBENCH_DEF(UnsegmentedCount) {
    size_t nsvps;
    sieving_prime** svps;

    primewheel pw;
    new_primewheel(&pw, wheel_size);

    int maxp = unsegmented_maxp;

    CBENCH_INIT("Unsegmented prime sieve (count) to %lg", (double)maxp) {
        CBENCH_MEASURE_START
        sieving_primes_to_n(&pw, maxp, &svps, &nsvps);
        CBENCH_MEASURE_STOP
        for (size_t i = 0; i < nsvps; i++) {
            free(svps[i]);
        }
    }

    free_primewheel(&pw);
}

int main(int argc, char *argv[])
{
    printf("sizeof(bitunit) = %llu\n", sizeof(bitunit));
    printf("buffer size: %llu\n", bufsize);
    printf("wheel size: %llu\n\n", wheel_size);

    do_bench_print_summary(CBENCH_UnsegmentedCount,
                           3,
                           3,
                           NULL);

    do_bench_print_summary(CBENCH_SegmentedCount,
                           4,
                           4,
                           NULL);
    return 0;
}
