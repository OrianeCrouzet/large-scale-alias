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
 * alias_internal.h -- data structures shared by the construction and the
 * sampling step (not part of the public interface).
 */

#ifndef ALIAS_INTERNAL_H
#define ALIAS_INTERNAL_H

#include <stdint.h>
#include <gmp.h>

#include "alias.h"

typedef unsigned __int128 u128;

/* Number of binary digits of T_i / cs stored in a cell. */
#define ALEA_DIGITS 31     /* ALIAS_ALEA: digits 1..31, plus a "finite expansion" flag */
#define SPEED_DIGITS 32    /* ALIAS_SPEED: digits 1..32 */

/*
 * One cell of the alias table (16 bytes).
 *
 * A cell holds a primary object and, possibly, an alias. Given the cell, the
 * primary object is returned with probability T_i / cs, the alias otherwise.
 * A cell without alias has alias == primary.
 *
 * digits, ALIAS_ALEA:  bits 31..1 = digits d_1..d_31 of T_i / cs (d_k at bit 32-k),
 *                      bit 0 = 1 if T_i 2^31 mod cs == 0 (finite expansion).
 * digits, ALIAS_SPEED: D_i = floor(T_i 2^32 / cs); UINT32_MAX for a cell without alias.
 * thr: position of the exact threshold T_i in the backend storage
 *      (index for u64 / u128, byte offset for GMP); cells with an alias only.
 */
typedef struct {
    uint32_t digits;
    uint32_t primary;
    uint32_t alias;
    uint32_t thr;
} alias_cell;

struct alias_table {
    alias_mode mode;
    alias_backend backend;
    uint32_t n;            /* objects */
    uint32_t q;            /* cells */
    uint32_t q_nbits;      /* bit length of q - 1 (first block of the Fast Dice Roller) */
    uint32_t q_reject;     /* (2^32 - q) mod q (Lemire's rejection threshold) */
    uint32_t virtual_obj;  /* index of the virtual object (n), UINT32_MAX if absent */
    uint32_t m;            /* cells with an alias: cells 0..m-1 */
    alias_cell *cells;     /* q cells */
    double reject_prob;    /* (cs - r) / (q cs): probability of drawing the virtual object */

    /* Exact thresholds T_i of cells 0..m-1 and the cell size cs, used only by
     * the slow paths (probability <= 2^-31 per attempt). One representation is used. */
    uint64_t cs64;
    uint64_t *thr64;
    u128 cs128;
    u128 *thr128;
    mpz_t cs_mpz;          /* initialized for every backend (alias_info) */
    unsigned char *thr_bytes;   /* GMP: little-endian significant bytes, cells in order */
    uint32_t thr_size;          /* GMP: bytes used in thr_bytes */
};

/* Slow paths, exposed for the tests.
 * alias_slow_tail: ALIAS_ALEA, the 31 first coins were 0 and the expansion is
 *   infinite; continues Lumbroso's Bernoulli from T_i 2^31 mod cs.
 * alias_slow_tie: ALIAS_SPEED, the 32 low bits of the word equal D_i; returns
 *   Bernoulli((T_i 2^32 mod cs) / cs) with fresh bits.
 * Both use the random-bit stream. */
uint32_t alias_slow_tail(const alias_table *t, uint32_t i);
uint32_t alias_slow_tie(const alias_table *t, uint32_t i);

/* T_i 2^k mod cs as a GMP integer (tests). */
void alias_threshold_remainder(const alias_table *t, uint32_t i, unsigned k, mpz_t out);

#endif
