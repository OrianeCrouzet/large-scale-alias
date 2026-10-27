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
 * build.c -- construction of the alias table (Algorithm 1) for the three
 * backends, and the exact slow paths of the sampling step.
 *
 * The three backends build the same table (same cells, in the same order);
 * they differ only in the type used for the weights and in the storage of the
 * exact thresholds needed by the slow paths.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "alias_internal.h"
#include "bitstream.h"

/* ------------------------------------------------------------------------- */
/* Conversions between GMP and machine integers                              */
/* ------------------------------------------------------------------------- */

static uint64_t u64_from_mpz(const mpz_t x) {
    uint64_t v = 0;
    mpz_export(&v, NULL, -1, sizeof(v), 0, 0, x);     /* x < 2^64 */
    return v;
}

static void u64_to_mpz(mpz_t out, uint64_t v) {
    mpz_import(out, 1, -1, sizeof(v), 0, 0, &v);
}

static u128 u128_from_mpz(const mpz_t x) {
    uint64_t limbs[2] = {0, 0};
    mpz_export(limbs, NULL, -1, sizeof(uint64_t), 0, 0, x);     /* x < 2^128 */
    return ((u128)limbs[1] << 64) | limbs[0];
}

static void u128_to_mpz(mpz_t out, u128 v) {
    uint64_t limbs[2] = {(uint64_t)v, (uint64_t)(v >> 64)};
    mpz_import(out, 2, -1, sizeof(uint64_t), 0, 0, limbs);
}

/* ------------------------------------------------------------------------- */
/* Native backends                                                           */
/* ------------------------------------------------------------------------- */

#define NT uint64_t
#define NT_SUFFIX u64
#define NT_FROM_MPZ u64_from_mpz
#define NT_TO_MPZ u64_to_mpz
#define NT_CS cs64
#define NT_THR thr64
#include "build_native.h"
#undef NT
#undef NT_SUFFIX
#undef NT_FROM_MPZ
#undef NT_TO_MPZ
#undef NT_CS
#undef NT_THR

#define NT u128
#define NT_SUFFIX u128
#define NT_FROM_MPZ u128_from_mpz
#define NT_TO_MPZ u128_to_mpz
#define NT_CS cs128
#define NT_THR thr128
#include "build_native.h"
#undef NT
#undef NT_SUFFIX
#undef NT_FROM_MPZ
#undef NT_TO_MPZ
#undef NT_CS
#undef NT_THR

/* ------------------------------------------------------------------------- */
/* GMP backend                                                               */
/* ------------------------------------------------------------------------- */

/* Exact threshold of cell i (cells with an alias), from its significant bytes. */
static void gmp_threshold(const alias_table *t, uint32_t i, mpz_t out) {
    uint32_t start = t->cells[i].thr;
    uint32_t end = (i + 1 < t->m) ? t->cells[i + 1].thr : t->thr_size;
    mpz_import(out, end - start, -1, 1, 0, 0, t->thr_bytes + start);
}

static void gmp_set_alone(alias_table *t, uint32_t *back, uint32_t x) {
    alias_cell *c = &t->cells[--*back];
    c->primary = x;
    c->alias = x;
    c->digits = (t->mode == ALIAS_SPEED) ? UINT32_MAX : 0;
    c->thr = 0;
}

static void gmp_set_cell(alias_table *t, uint32_t *back, size_t *capacity, const mpz_t tw,
                         uint32_t x2, uint32_t x, mpz_t tmp, mpz_t rem) {
    if (mpz_sgn(tw) == 0) {
        gmp_set_alone(t, back, x);
        return;
    }
    alias_cell *c = &t->cells[t->m];
    c->primary = x2;
    c->alias = x;
    size_t len = (mpz_sizeinbase(tw, 2) + 7) / 8;
    if (t->thr_size + len > *capacity) {
        *capacity = 2 * (t->thr_size + len);
        t->thr_bytes = realloc(t->thr_bytes, *capacity);
    }
    mpz_export(t->thr_bytes + t->thr_size, &len, -1, 1, 0, 0, tw);
    c->thr = t->thr_size;
    t->thr_size += (uint32_t)len;
    unsigned k = (t->mode == ALIAS_SPEED) ? SPEED_DIGITS : ALEA_DIGITS;
    mpz_mul_2exp(tmp, tw, k);
    mpz_fdiv_qr(tmp, rem, tmp, t->cs_mpz);       /* digits = floor(tw 2^k / cs) */
    uint32_t digits = (uint32_t)mpz_get_ui(tmp);
    if (t->mode == ALIAS_SPEED) {
        c->digits = digits;
    } else {
        c->digits = (digits << (32 - ALEA_DIGITS)) | (mpz_sgn(rem) == 0);
    }
    t->m++;
}

static int gmp_build(const alias_dist *d, alias_table *t) {
    uint32_t n = d->n;
    mpz_t *D = malloc(((size_t)n + 1) * sizeof(mpz_t));
    for (uint32_t i = 0; i < n; i++) {
        mpz_init_set(D[i], d->weights[i]);
    }
    mpz_t q, r, tmp, rem;
    mpz_inits(q, r, tmp, rem, NULL);
    mpz_fdiv_q_ui(t->cs_mpz, d->total, n);
    mpz_fdiv_qr(q, r, d->total, t->cs_mpz);
    uint32_t nobj = n;
    t->virtual_obj = UINT32_MAX;
    if (mpz_sgn(r) > 0) {
        mpz_init(D[n]);
        mpz_sub(D[n], t->cs_mpz, r);
        t->virtual_obj = n;
        nobj = n + 1;
        mpz_add_ui(q, q, 1);
    }
    int status = -1;
    if (!mpz_fits_uint_p(q) || mpz_get_ui(q) >= UINT32_MAX) {
        goto done;
    }
    t->q = (uint32_t)mpz_get_ui(q);
    t->cells = malloc((size_t)t->q * sizeof(alias_cell));
    size_t capacity = (size_t)nobj * ((mpz_sizeinbase(t->cs_mpz, 2) + 7) / 8) + 1;
    t->thr_bytes = malloc(capacity);
    t->thr_size = 0;
    t->m = 0;
    uint32_t back = t->q;

    uint32_t *small = malloc((size_t)nobj * sizeof(uint32_t));
    uint32_t *large = malloc((size_t)nobj * sizeof(uint32_t));
    uint32_t n_small = 0, n_large = 0;
    for (uint32_t i = 0; i < nobj; i++) {
        if (mpz_cmp(D[i], t->cs_mpz) > 0) {
            large[n_large++] = i;
        } else {
            small[n_small++] = i;
        }
    }
    while (n_large > 0) {
        uint32_t x = large[--n_large];
        if (n_small > 0) {
            uint32_t x2 = small[--n_small];
            if (mpz_cmp(D[x2], t->cs_mpz) != 0) {
                gmp_set_cell(t, &back, &capacity, D[x2], x2, x, tmp, rem);
                mpz_sub(D[x], D[x], t->cs_mpz);     /* D[x] -= cs - D[x2] */
                mpz_add(D[x], D[x], D[x2]);
            } else {
                gmp_set_alone(t, &back, x2);
            }
        } else {
            gmp_set_alone(t, &back, x);
            mpz_sub(D[x], D[x], t->cs_mpz);
        }
        if (mpz_cmp(D[x], t->cs_mpz) > 0) {
            large[n_large++] = x;
        } else {
            small[n_small++] = x;
        }
    }
    while (n_small > 0) {
        gmp_set_alone(t, &back, small[--n_small]);
    }
    free(small);
    free(large);
    if (back == t->m) {
        t->thr_bytes = realloc(t->thr_bytes, (size_t)t->thr_size + 1);
        if (mpz_sgn(r) > 0) {
            mpz_sub(tmp, t->cs_mpz, r);
            mpz_mul_ui(rem, t->cs_mpz, t->q);
            t->reject_prob = mpz_get_d(tmp) / mpz_get_d(rem);
        }
        status = 0;
    }
done:
    for (uint32_t i = 0; i < nobj; i++) {
        mpz_clear(D[i]);
    }
    free(D);
    mpz_clears(q, r, tmp, rem, NULL);
    return status;
}

static uint32_t gmp_lumbroso(mpz_t rem, const mpz_t cs) {
    if (mpz_sgn(rem) == 0) {
        return 0;
    }
    for (;;) {
        mpz_mul_2exp(rem, rem, 1);
        int c = mpz_cmp(rem, cs);
        if (c == 0) {
            return bs_bit();
        }
        uint32_t y = (c > 0);
        if (y) {
            mpz_sub(rem, rem, cs);
        }
        if (bs_bit()) {
            return y;
        }
    }
}

/* ------------------------------------------------------------------------- */
/* Public construction interface                                             */
/* ------------------------------------------------------------------------- */

static uint32_t bit_length(uint32_t x) {
    return x ? 32 - (uint32_t)__builtin_clz(x) : 0;
}

alias_table *alias_build(const alias_dist *d, alias_mode mode, alias_backend backend) {
    alias_backend min = alias_min_backend(d);
    if (backend == ALIAS_BACKEND_AUTO) {
        backend = min;
    }
    if (backend < min) {
        fprintf(stderr, "alias_build: backend %s too small (total weight of %zu bits needs %s)\n",
                alias_backend_name(backend), mpz_sizeinbase(d->total, 2), alias_backend_name(min));
        return NULL;
    }
    alias_table *t = calloc(1, sizeof(alias_table));
    t->mode = mode;
    t->backend = backend;
    t->n = d->n;
    mpz_init(t->cs_mpz);

    int status;
    switch (backend) {
    case ALIAS_BACKEND_U64:  status = build_u64(d, t); break;
    case ALIAS_BACKEND_U128: status = build_u128(d, t); break;
    default:                 status = gmp_build(d, t); break;
    }
    if (status != 0) {
        fprintf(stderr, "alias_build: table too large\n");
        alias_table_free(t);
        return NULL;
    }
    t->q_nbits = (t->q > 1) ? bit_length(t->q - 1) : 0;
    t->q_reject = (uint32_t)(-t->q) % t->q;
    return t;
}

void alias_table_free(alias_table *t) {
    if (t == NULL) {
        return;
    }
    free(t->cells);
    free(t->thr64);
    free(t->thr128);
    free(t->thr_bytes);
    mpz_clear(t->cs_mpz);
    free(t);
}

void alias_info(const alias_table *t, alias_table_info *info) {
    info->mode = t->mode;
    info->backend = t->backend;
    info->n = t->n;
    info->q = t->q;
    info->aliased = t->m;
    info->has_virtual = (t->virtual_obj != UINT32_MAX);
    info->reject_prob = t->reject_prob;
    size_t bytes = sizeof(*t) + (size_t)t->q * sizeof(alias_cell);
    switch (t->backend) {
    case ALIAS_BACKEND_U64:  bytes += (size_t)t->m * sizeof(uint64_t); break;
    case ALIAS_BACKEND_U128: bytes += (size_t)t->m * sizeof(u128); break;
    default:
        bytes += t->thr_size + (mpz_sizeinbase(t->cs_mpz, 2) + 7) / 8;
        break;
    }
    info->bytes = bytes;
}

/* ------------------------------------------------------------------------- */
/* Slow paths                                                                */
/* ------------------------------------------------------------------------- */

static uint32_t slow_bernoulli(const alias_table *t, uint32_t i, unsigned k) {
    switch (t->backend) {
    case ALIAS_BACKEND_U64:
        return lumbroso_u64(remainder_u64(t, i, k), t->cs64);
    case ALIAS_BACKEND_U128:
        return lumbroso_u128(remainder_u128(t, i, k), t->cs128);
    default: {
        mpz_t rem;
        mpz_init(rem);
        gmp_threshold(t, i, rem);
        mpz_mul_2exp(rem, rem, k);
        mpz_fdiv_r(rem, rem, t->cs_mpz);
        uint32_t y = gmp_lumbroso(rem, t->cs_mpz);
        mpz_clear(rem);
        return y;
    }
    }
}

uint32_t alias_slow_tail(const alias_table *t, uint32_t i) {
    return slow_bernoulli(t, i, ALEA_DIGITS);
}

uint32_t alias_slow_tie(const alias_table *t, uint32_t i) {
    const alias_cell *c = &t->cells[i];
    if (c->alias == c->primary) {
        return 1;                        /* cell without alias */
    }
    return slow_bernoulli(t, i, SPEED_DIGITS);
}

void alias_threshold_remainder(const alias_table *t, uint32_t i, unsigned k, mpz_t out) {
    switch (t->backend) {
    case ALIAS_BACKEND_U64:  u64_to_mpz(out, t->thr64[t->cells[i].thr]); break;
    case ALIAS_BACKEND_U128: u128_to_mpz(out, t->thr128[t->cells[i].thr]); break;
    default:                 gmp_threshold(t, i, out); break;
    }
    mpz_mul_2exp(out, out, k);
    mpz_fdiv_r(out, out, t->cs_mpz);
}

/* ------------------------------------------------------------------------- */
/* Random-bit stream                                                         */
/* ------------------------------------------------------------------------- */

void alias_rng_start(void) {
    bs_start();
}

void alias_rng_stop(void) {
    bs_stop();
}

uint64_t alias_rng_bits(void) {
    return bs_bits_consumed();
}
