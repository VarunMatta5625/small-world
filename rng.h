#ifndef RNG_H
#define RNG_H
/*
 * Simple Linear Congruential Generator (LCG) for reproducible randomness.
 * Plain C: int state, no bitwise operations required by API users.
 *
 * Note: Implementation uses modulo arithmetic with 2^31 (see rng.c).
 */

typedef struct {
    int state; /* must be non-zero */
} RNG;

/* Initialize RNG with a seed (seed==0 is promoted to 1). */
void rng_init(RNG* rng, int seed);

/* Next raw integer in [0, 2^31). */
int rng_next(RNG* rng);

/* Random float in [0,1). */
float rng_float(RNG* rng);

/* Random integer in [0, n). If n<=0, returns 0. */
int rng_range(RNG* rng, int n);

#endif /* RNG_H */
