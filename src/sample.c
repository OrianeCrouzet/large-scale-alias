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
 * sample.c -- sampling step, common to the three backends.
 *
 * Each attempt draws a cell uniformly among the q cells, then a Bernoulli
 * variable choosing between the primary object and the alias of the cell.
 * Attempts are repeated while they yield the virtual object.
 *
 *   ALIAS_ALEA  (v1): Fast Dice Roller for the cell and Lumbroso's Bernoulli,
 *                     both optimal in random bits (Lumbroso 2013).
 *   ALIAS_SPEED (v2): one 64-bit word per attempt: Lemire's multiply-shift on
 *                     the high half for the cell, comparison of the low half
 *                     with 32 digits of T_i / cs for the Bernoulli.
 *
 * Neither mode uses the exact thresholds, except in slow paths of probability
 * at most 2^-31 per attempt (alias_slow_tail, alias_slow_tie).
 */

#include "alias_internal.h"
#include "bitstream.h"

/* Fast Dice Roller (Lumbroso) on [0, n); nbits = bit length of n - 1. */
static inline uint32_t fdr_uniform(uint32_t n, uint32_t nbits) {
    if (n <= 1) {
        return 0;
    }
    uint64_t bound = 1ull << nbits;
    uint64_t x = bs_bits(nbits);
    while (bound < n || x >= n) {
        if (bound >= n) {
            bound -= n;
            x -= n;
        }
        bound <<= 1;
        x = (x << 1) | bs_bit();
    }
    return (uint32_t)x;
}

/* Lumbroso's Bernoulli(T_i / cs) on the 31 stored digits: reads coins until the
 * first 1, at position k, and returns d_k. If the expansion ends at digit L, at
 * most L coins are read (the k = L case returns c_L = 1, and L zeros return 0).
 * Returns 2 when the 31 coins are 0 and the expansion goes on (slow path). */
static inline uint32_t bernoulli_digits(uint32_t bits) {
    uint32_t finite = bits & 1;
    uint32_t cap = finite ? 32 - (uint32_t)__builtin_ctz(bits & ~1u) : ALEA_DIGITS;
    uint32_t k = bs_coins_until_one(cap);
    if (k) {
        return (bits >> (32 - k)) & 1;
    }
    return finite ? 0 : 2;
}

static uint32_t sample_alea(const alias_table *t) {
    for (;;) {
        uint32_t i = fdr_uniform(t->q, t->q_nbits);
        const alias_cell *c = &t->cells[i];
        uint32_t result = c->primary;
        if (c->alias != c->primary) {
            uint32_t b = bernoulli_digits(c->digits);
            if (__builtin_expect(b == 2, 0)) {
                b = alias_slow_tail(t, i);
            }
            result = b ? c->primary : c->alias;
        }
        if (result != t->virtual_obj) {
            return result;
        }
    }
}

/* Lemire's rejection (probability < q / 2^32): fresh 32-bit words until accepted. */
static __attribute__((noinline)) uint32_t lemire_retry(const alias_table *t, uint64_t m) {
    while ((uint32_t)m < t->q_reject) {
        m = bs_bits(32) * (uint64_t)t->q;
    }
    return (uint32_t)(m >> 32);
}

static uint32_t sample_speed(const alias_table *t) {
    for (;;) {
        uint64_t w = bs_word();
        uint64_t m = (w >> 32) * (uint64_t)t->q;
        uint32_t i = (uint32_t)(m >> 32);
        if (__builtin_expect((uint32_t)m < t->q, 0)) {
            i = lemire_retry(t, m);
        }
        const alias_cell *c = &t->cells[i];
        uint32_t lo = (uint32_t)w;
        /* V = (lo + U) / 2^32, U uniform on [0, 1): lo < D_i implies V < T_i / cs,
         * lo > D_i implies V >= T_i / cs, lo = D_i is decided exactly. */
        uint32_t b = lo < c->digits;
        if (__builtin_expect(lo == c->digits, 0)) {
            b = alias_slow_tie(t, i);
        }
        uint32_t result = b ? c->primary : c->alias;
        if (result != t->virtual_obj) {
            return result;
        }
    }
}

uint32_t alias_sample(const alias_table *t) {
    return (t->mode == ALIAS_SPEED) ? sample_speed(t) : sample_alea(t);
}
