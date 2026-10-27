/*
 * Copyright (C) 2026 Oriane Crouzet, Antoine Genitrini, Mehdi Naima
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version. See the LICENSE file.
 */
/*
 * alias.h -- Integer-Alias Sampler: exact sampling from a discrete distribution
 * given by non-negative integer weights of arbitrary size.
 *
 * Typical use:
 *
 *     alias_dist d;
 *     alias_dist_read("weights.dist", &d, err, sizeof err);
 *     alias_report r;
 *     alias_analyze(&d, &r);                    // chooses the arithmetic backend
 *     alias_table *t = alias_build(&d, ALIAS_SPEED, r.backend);
 *     alias_rng_start();
 *     uint32_t x = alias_sample(t);             // index of an object, 0 <= x < n
 *     alias_rng_stop();
 *     alias_table_free(t);
 *     alias_dist_free(&d);
 *
 * See TECHNICAL.md for the algorithms and the data structures.
 */

#ifndef ALIAS_H
#define ALIAS_H

#include <stddef.h>
#include <stdint.h>
#include <gmp.h>

/* What the sampling step optimizes. Both modes sample the exact distribution. */
typedef enum {
    ALIAS_SPEED,   /* running time: one 64-bit random word per attempt (v2) */
    ALIAS_ALEA     /* random bits: Lumbroso's bit-optimal uniform and Bernoulli (v1) */
} alias_mode;

/* Arithmetic used to build the table and to resolve the rare exact slow paths. */
typedef enum {
    ALIAS_BACKEND_AUTO,  /* chosen by alias_analyze */
    ALIAS_BACKEND_U64,   /* uint64_t: total weight < 2^64 */
    ALIAS_BACKEND_U128,  /* unsigned __int128: total weight < 2^128 */
    ALIAS_BACKEND_GMP    /* GMP integers: no bound */
} alias_backend;

/* A distribution read from a .dist file. The weights are divided by their gcd,
 * then multiplied by n if their sum is smaller than n (so that the cell size
 * cs = floor(w / n) is positive). Neither step changes the distribution. */
typedef struct {
    uint32_t n;            /* number of objects */
    mpz_t *weights;        /* n weights (after normalization) */
    mpz_t total;           /* sum of the weights (after normalization) */
    uint32_t scale;        /* 1, or n if the weights were multiplied by n */
    mpz_t declared_total;  /* first line of the file (not used for sampling) */
    mpz_t gcd;             /* gcd of the original weights */
    int has_entropy;       /* the file has a third field */
    double declared_entropy;
} alias_dist;

/* Result of the analysis run before building the table. */
typedef struct {
    uint32_t n;
    uint32_t n_zero;               /* objects of weight 0 (never sampled) */
    size_t total_bits;             /* bit length of the total weight w */
    size_t max_weight_bits;        /* bit length of the largest weight */
    size_t gcd_bits;               /* bit length of the gcd removed from the weights */
    uint32_t scale;                /* factor applied after the gcd (1 or n) */
    size_t cs_bits;                /* bit length of the cell size cs = floor(w / n) */
    int declared_total_matches;    /* first line of the file == sum of the original weights */
    double entropy;                /* Shannon entropy (bits), computed from the weights */
    uint32_t below_double;         /* objects with n * w_i / w < 2^-53: never drawn by a
                                      double-precision Alias sampler */
    alias_backend backend;         /* smallest exact backend */
} alias_report;

/* Reads a .dist file: the total weight, then n followed by the n weights, then
 * optionally the entropy. Returns 0 on success, -1 on error (message in err). */
int alias_dist_read(const char *path, alias_dist *d, char *err, size_t errlen);
void alias_dist_free(alias_dist *d);

/* Analyzes the distribution and fills r, including the chosen backend. */
void alias_analyze(const alias_dist *d, alias_report *r);

/* Smallest backend able to represent every intermediate value exactly. */
alias_backend alias_min_backend(const alias_dist *d);

const char *alias_backend_name(alias_backend b);
const char *alias_mode_name(alias_mode m);

typedef struct alias_table alias_table;

/* Builds the alias table (Algorithm 1 of the paper). backend must be
 * ALIAS_BACKEND_AUTO or at least alias_min_backend(d). Returns NULL on error. */
alias_table *alias_build(const alias_dist *d, alias_mode mode, alias_backend backend);
void alias_table_free(alias_table *t);

/* Draws one object, 0 <= index < n. Requires alias_rng_start(). Not thread-safe. */
uint32_t alias_sample(const alias_table *t);

/* Information on a built table. */
typedef struct {
    alias_mode mode;
    alias_backend backend;
    uint32_t n;            /* objects */
    uint32_t q;            /* cells (n or n + 1) */
    uint32_t aliased;      /* cells with an alias */
    int has_virtual;       /* a virtual object completes the last cell */
    double reject_prob;    /* probability that an attempt yields the virtual object */
    size_t bytes;          /* memory of the table */
} alias_table_info;

void alias_info(const alias_table *t, alias_table_info *info);

/* Random-bit stream shared by all tables (one consumer thread). */
void alias_rng_start(void);
void alias_rng_stop(void);
uint64_t alias_rng_bits(void);   /* random bits consumed since alias_rng_start */

#endif
