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
 * build_native.h -- construction of the alias table on a machine integer type.
 *
 * Included by build.c once per type, with the following macros defined:
 *   NT            the unsigned integer type (uint64_t or u128)
 *   NT_SUFFIX     suffix of the generated function names (u64 or u128)
 *   NT_FROM_MPZ   NT NT_FROM_MPZ(const mpz_t x), x < 2^(bits of NT)
 *   NT_TO_MPZ     void NT_TO_MPZ(mpz_t out, NT v)
 *   NT_CS, NT_THR the fields of alias_table receiving cs and the thresholds
 *
 * Every value handled here is at most the total weight w < 2^(bits of NT),
 * and the doublings 2t are applied to t < cs <= w / 2, so nothing overflows.
 */

#define CAT2(a, b) a##_##b
#define CAT(a, b) CAT2(a, b)
#define FN(name) CAT(name, NT_SUFFIX)

/* Digits d_1..d_k of t / cs (0 < t < cs), packed with d_j at bit 32 - j;
 * *rem = t 2^k mod cs. Binary long division. */
static uint32_t FN(digits)(NT t, NT cs, unsigned k, NT *rem) {
    uint32_t bits = 0;
    for (unsigned j = 1; j <= k; j++) {
        t <<= 1;
        if (t >= cs) {
            t -= cs;
            bits |= 1u << (32 - j);
        }
    }
    *rem = t;
    return bits;
}

/* Cell (x2, x): x2 of weight tw < cs is returned with probability tw / cs. */
static void FN(set_cell)(alias_table *t, NT *thr, uint32_t *back, NT tw, NT cs,
                         uint32_t x2, uint32_t x) {
    if (tw == 0) {
        /* Bernoulli(0) = 0 without bits: the cell always yields x. */
        alias_cell *c = &t->cells[--*back];
        c->primary = x;
        c->alias = x;
        c->digits = (t->mode == ALIAS_SPEED) ? UINT32_MAX : 0;
        c->thr = 0;
        return;
    }
    alias_cell *c = &t->cells[t->m];
    c->primary = x2;
    c->alias = x;
    c->thr = t->m;
    thr[t->m] = tw;
    NT rem;
    if (t->mode == ALIAS_SPEED) {
        c->digits = FN(digits)(tw, cs, SPEED_DIGITS, &rem);
    } else {
        c->digits = FN(digits)(tw, cs, ALEA_DIGITS, &rem) | (rem == 0);
    }
    t->m++;
}

static void FN(set_alone)(alias_table *t, uint32_t *back, uint32_t x) {
    alias_cell *c = &t->cells[--*back];
    c->primary = x;
    c->alias = x;
    c->digits = (t->mode == ALIAS_SPEED) ? UINT32_MAX : 0;
    c->thr = 0;
}

/* Algorithm 1: stacks L (weight > cs) and S (weight <= cs), LIFO. */
static int FN(build)(const alias_dist *d, alias_table *t) {
    uint32_t n = d->n;
    NT *D = malloc(((size_t)n + 1) * sizeof(NT));
    NT w = 0;
    for (uint32_t i = 0; i < n; i++) {
        D[i] = NT_FROM_MPZ(d->weights[i]);
        w += D[i];
    }
    NT cs = w / n;                       /* cell size, >= 1 (alias_dist_read) */
    NT qq = w / cs;
    NT r = w % cs;
    uint32_t nobj = n;
    t->virtual_obj = UINT32_MAX;
    if (r > 0) {
        D[n] = cs - r;                   /* virtual object completing the last cell */
        t->virtual_obj = n;
        nobj = n + 1;
        qq++;
    }
    if (qq >= UINT32_MAX) {
        free(D);
        return -1;
    }
    t->q = (uint32_t)qq;
    t->cells = malloc((size_t)t->q * sizeof(alias_cell));
    NT *thr = malloc((size_t)t->q * sizeof(NT));
    t->m = 0;
    uint32_t back = t->q;

    uint32_t *small = malloc((size_t)nobj * sizeof(uint32_t));
    uint32_t *large = malloc((size_t)nobj * sizeof(uint32_t));
    uint32_t n_small = 0, n_large = 0;
    for (uint32_t i = 0; i < nobj; i++) {
        if (D[i] > cs) {
            large[n_large++] = i;
        } else {
            small[n_small++] = i;
        }
    }
    while (n_large > 0) {
        uint32_t x = large[--n_large];
        if (n_small > 0) {
            uint32_t x2 = small[--n_small];
            if (D[x2] != cs) {
                FN(set_cell)(t, thr, &back, D[x2], cs, x2, x);
                D[x] -= cs - D[x2];      /* x fills the rest of the cell */
            } else {
                FN(set_alone)(t, &back, x2);
            }
        } else {
            FN(set_alone)(t, &back, x);
            D[x] -= cs;
        }
        if (D[x] > cs) {
            large[n_large++] = x;
        } else {
            small[n_small++] = x;
        }
    }
    while (n_small > 0) {
        FN(set_alone)(t, &back, small[--n_small]);
    }
    free(small);
    free(large);
    free(D);
    if (back != t->m) {
        free(thr);
        return -1;                       /* cannot happen (Theorem 1) */
    }

    t->NT_CS = cs;
    t->NT_THR = realloc(thr, ((size_t)t->m + 1) * sizeof(NT));
    NT_TO_MPZ(t->cs_mpz, cs);
    t->reject_prob = (r > 0) ? (double)(cs - r) / ((double)t->q * (double)cs) : 0.0;
    return 0;
}

/* T_i 2^k mod cs (0 < T_i < cs). */
static NT FN(remainder)(const alias_table *t, uint32_t i, unsigned k) {
    NT rem;
    FN(digits)(t->NT_THR[t->cells[i].thr], t->NT_CS, k, &rem);
    return rem;
}

/* Lumbroso's Bernoulli continued from the remainder rem < cs: at each step,
 * rem = 2 rem; the digit is rem >= cs; a fresh coin equal to 1 stops with that digit. */
static uint32_t FN(lumbroso)(NT rem, NT cs) {
    if (rem == 0) {
        return 0;
    }
    for (;;) {
        rem <<= 1;
        if (rem == cs) {
            return bs_bit();             /* remaining expansion 0.1000... */
        }
        uint32_t y = (rem > cs);
        if (y) {
            rem -= cs;
        }
        if (bs_bit()) {
            return y;
        }
    }
}

#undef FN
#undef CAT
#undef CAT2
