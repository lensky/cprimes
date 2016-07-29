extern "C" {
    #include "cprimes/primewheel.h"
    #include "cprimes/sieve.h"
    #include "cprimes/cprimes.h"
    #include "bitarray.h"
}
#include <cstdlib>
#include <ctime>

#include <iostream>

#include "gtest/gtest.h"

prime KNOWN_LOW_PRIMES[] = {2,3,5,7,11};

TEST(WheelTest, CorrectPrimes) {
    primewheel tw;
    new_primewheel(&tw, 4);
    for (size_t i = 0; i < tw.nprimes; i++) {
        ASSERT_EQ(KNOWN_LOW_PRIMES[i], tw.primes[i]);
    }
    free_primewheel(&tw);
}

prime W3_COPS[] = {1,7,11,13,17,19,23,29};

TEST(WheelTest, CorrectDeltas) {
    primewheel tw;
    new_primewheel(&tw, 3);
    ASSERT_EQ(8, tw.spokes);
    for (size_t i = 0; i < tw.spokes; i++) {
        ASSERT_EQ(W3_COPS[i], tw.wp_coprimes[i]);
    }
}

const size_t ncixtests = 10;

TEST(WheelTest, FindClosestIx) {
    primewheel pw;
    new_primewheel(&pw, 4);
    for (size_t i = 0; i < ncixtests; i++) {
        prime r = rand();
        size_t ix = closest_w_ixle(&pw, r);
        ASSERT_LE(w_coprime_of_ix(&pw, ix), r);
        ASSERT_LT(r, w_coprime_of_ix(&pw, ix + 1));
    }
    free_primewheel(&pw);
}

const size_t nixmtests = 10;

TEST(WheelTest, IxMapping) {
    primewheel pw;
    new_primewheel(&pw, 4);

    for (size_t i = 0; i < nixmtests; i++) {
        size_t p = w_coprime_of_ix(&pw, i);
        size_t pix = w_ix_of_coprime(&pw, p);

        ASSERT_EQ(i, pix);
    }

    free_primewheel(&pw);
}

const prime P_LIMITS[] = {1e2, 1e3, 1e4, 1e5,   1e6,  1e7,     1e8,     1e9};
const prime P_COUNTS[] = {25, 168, 1229, 9592, 78498, 664579, 5761455, 50847534};

#define SINGLE_SIEVE_LIMIT 1e7
#define SINGLE_SIEVE_IX 5

TEST(SingleSieveTest, Simple) {
    primewheel pw;
    new_primewheel(&pw, 3);

    sieving_prime** svps;
    size_t nsvps;

    size_t ssix = 0;
    for (; ssix < sizeof(P_LIMITS) / sizeof(P_LIMITS[0]); ssix++) {
        if (P_LIMITS[ssix] >= SINGLE_SIEVE_LIMIT) {
            break;
        }
    }
    sieving_primes_to_n(&pw, P_LIMITS[ssix], &svps, &nsvps);

    ASSERT_EQ(P_COUNTS[ssix] - pw.nprimes, nsvps);

    for (size_t i = 0; i < ssix; i++) {
        size_t count = pw.nprimes;
        for (size_t j = 0; j < nsvps; j++) {
            if (svps[j]->prime <= P_LIMITS[i]) {
                count++;
            }
        }
        ASSERT_EQ(P_COUNTS[i], count);
    }

    for (size_t i = 0; i < nsvps; i++) {
        free(svps[i]);
    }
    free(svps);
    free_primewheel(&pw);
}

TEST(SegmentedSieveTest, BasicCount) {
    primewheel pw;
    new_primewheel(&pw, 4);

    siever svr;
    new_siever(&svr, &pw, 1e9);

    sieve segment;
    new_sieve(&segment, DEF_BUFSIZE, 0);

    for (size_t i = 0; i < sizeof(P_COUNTS) / sizeof(P_COUNTS[0]); i++) {
        ASSERT_EQ(P_COUNTS[i],
                  segmented_count_primes(&svr, &segment, 0, P_LIMITS[i]));
    }

    free_sieve(&segment);
    free_siever(&svr);
    free_primewheel(&pw);
}

TEST(SegmentedSieveTest, ExtendSiever) {
    primewheel pw;
    new_primewheel(&pw, 4);

    siever svr1, svr2;
    new_siever(&svr1, &pw, 1e7);
    new_siever(&svr2, &pw, 1e9);

    sieve segment; new_sieve(&segment, DEF_BUFSIZE, 0);

    extend_siever(&svr1, &segment, 1e9);

    ASSERT_EQ(svr2.n_sieving_primes, svr1.n_sieving_primes);
    for (size_t i = 0; i < svr2.n_sieving_primes; i++) {
        sieving_prime *svp1, *svp2;
        svp1 = svr1.sieving_primes[i];
        svp2 = svr2.sieving_primes[i];
        ASSERT_EQ(svp2->prime, svp1->prime);
        ASSERT_EQ(svp2->offset, svp1->offset);
    }

    free_siever(&svr1); free_siever(&svr2);
    free_primewheel(&pw);
}

TEST(SegmentedSieveTest, List) {
    primewheel pw;
    new_primewheel(&pw, 4);

    siever svr;
    new_siever(&svr, &pw, 1e9);

    sieve segment;
    new_sieve(&segment, DEF_BUFSIZE, 0);

    prime* primes;
    size_t nprimes = 1;

    segmented_list_primes(&svr, &segment, 0, P_LIMITS[7], &primes, &nprimes);
    ASSERT_EQ(P_COUNTS[7], nprimes);
    ASSERT_EQ(17329489, primes[1111110]);

    free(primes);
    free_sieve(&segment);
    free_siever(&svr);
    free_primewheel(&pw);
}

TEST(SegmentedSieveTest, SubtractCount) {
    primewheel pw; new_primewheel(&pw, 4);
    siever svr; new_siever(&svr, &pw, 1e9);
    sieve segment; new_sieve(&segment, DEF_BUFSIZE, 0);

    size_t ncounts = sizeof(P_COUNTS) / sizeof(P_COUNTS[0]);
    size_t large_count = segmented_count_primes(&svr,
                                                &segment,
                                                0,
                                                P_LIMITS[ncounts - 1]);

    for (size_t i = 0; i < ncounts - 1; i++) {
        ASSERT_EQ(large_count - P_COUNTS[i],
                  segmented_count_primes(&svr, &segment,
                                         P_LIMITS[i], P_LIMITS[ncounts - 1]))
            << "failed limit: P_LIMITS[" << i << "] = "
            << std::scientific << (double)P_LIMITS[i];
    }

    free_sieve(&segment);
    free_siever(&svr);
    free_primewheel(&pw);
}

int main(int argc, char *argv[])
{
    srand(time(NULL));
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
