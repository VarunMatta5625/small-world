#include "rng.h"

/*
 * LCG parameters (classic textbook choice):
 *   next = (a * x + c) mod m
 * with a = 1103515245, c = 12345, m = 2^31
 * This is NOT cryptographic but is simple, portable, and deterministic.
 */

static int rng_next_internal(RNG* rng) {
    long a = 1103515245L;
    long c = 12345L;
    long m = 2147483648L; /* 2^31 */
    long s = (long)rng->state;
    s = (a * s + c) % m;
    rng->state = (int)s;
    return rng->state;
}

void rng_init(RNG* rng, int seed) {
    if (!rng) return;
    if (seed == 0) seed = 1;
    rng->state = seed;
}

int rng_next(RNG* rng) {
    if (!rng) return 0;
    return rng_next_internal(rng);
}

float rng_float(RNG* rng) {
    int x = rng_next(rng);
    return (float)x / 2147483648.0f;
}

int rng_range(RNG* rng, int n) {
    if (!rng || n <= 0) return 0;
    float f = rng_float(rng);
    int v = (int)(f * (float)n);
    if (v >= n) v = n - 1;
    if (v < 0) v = 0;
    return v;
}
