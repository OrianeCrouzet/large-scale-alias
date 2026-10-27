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
 * dist.c -- reading and analysis of a .dist file.
 *
 * Format (whitespace-separated decimal integers):
 *     Z                 declared total weight (kept for information only)
 *     n w_1 ... w_n     number of objects, then their non-negative weights
 *     [H]               optional: entropy of the distribution
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "alias.h"

static int fail(char *err, size_t errlen, FILE *fp, const char *fmt, const char *path) {
    if (err != NULL && errlen > 0) {
        snprintf(err, errlen, fmt, path);
    }
    if (fp != NULL) {
        fclose(fp);
    }
    return -1;
}

int alias_dist_read(const char *path, alias_dist *d, char *err, size_t errlen) {
    memset(d, 0, sizeof(*d));
    FILE *fp = fopen(path, "r");
    if (fp == NULL) {
        return fail(err, errlen, NULL, "cannot open %s", path);
    }

    mpz_inits(d->declared_total, d->total, d->gcd, NULL);
    mpz_t v;
    mpz_init(v);
    if (mpz_inp_str(d->declared_total, fp, 10) == 0) {
        mpz_clear(v);
        mpz_clears(d->declared_total, d->total, d->gcd, NULL);
        return fail(err, errlen, fp, "%s: cannot read the total weight (first field)", path);
    }
    if (mpz_inp_str(v, fp, 10) == 0 || mpz_sgn(v) <= 0 || mpz_cmp_ui(v, UINT32_MAX - 2) > 0) {
        mpz_clear(v);
        mpz_clears(d->declared_total, d->total, d->gcd, NULL);
        return fail(err, errlen, fp, "%s: invalid number of objects (second field)", path);
    }
    d->n = (uint32_t)mpz_get_ui(v);
    mpz_clear(v);

    d->weights = malloc((size_t)d->n * sizeof(mpz_t));
    for (uint32_t i = 0; i < d->n; i++) {
        mpz_init(d->weights[i]);
        if (mpz_inp_str(d->weights[i], fp, 10) == 0 || mpz_sgn(d->weights[i]) < 0) {
            for (uint32_t j = 0; j <= i; j++) {
                mpz_clear(d->weights[j]);
            }
            free(d->weights);
            mpz_clears(d->declared_total, d->total, d->gcd, NULL);
            return fail(err, errlen, fp, "%s: missing or negative weight", path);
        }
        mpz_gcd(d->gcd, d->gcd, d->weights[i]);
    }
    d->has_entropy = (fscanf(fp, "%lf", &d->declared_entropy) == 1);
    fclose(fp);

    if (mpz_sgn(d->gcd) == 0) {
        alias_dist_free(d);
        return fail(err, errlen, NULL, "%s: all the weights are 0", path);
    }
    /* Dividing by the gcd leaves the distribution unchanged and may allow a smaller backend. */
    for (uint32_t i = 0; i < d->n; i++) {
        mpz_divexact(d->weights[i], d->weights[i], d->gcd);
        mpz_add(d->total, d->total, d->weights[i]);
    }
    /* The construction needs cs = floor(w / n) >= 1. */
    d->scale = 1;
    if (mpz_cmp_ui(d->total, d->n) < 0) {
        d->scale = d->n;
        for (uint32_t i = 0; i < d->n; i++) {
            mpz_mul_ui(d->weights[i], d->weights[i], d->n);
        }
        mpz_mul_ui(d->total, d->total, d->n);
    }
    return 0;
}

void alias_dist_free(alias_dist *d) {
    if (d->weights != NULL) {
        for (uint32_t i = 0; i < d->n; i++) {
            mpz_clear(d->weights[i]);
        }
        free(d->weights);
        d->weights = NULL;
    }
    mpz_clears(d->declared_total, d->total, d->gcd, NULL);
}

alias_backend alias_min_backend(const alias_dist *d) {
    /* Every value handled by the construction and the slow paths is at most the
     * total weight w (see TECHNICAL.md, "Backend selection"). */
    size_t bits = mpz_sizeinbase(d->total, 2);
    if (bits <= 64) {
        return ALIAS_BACKEND_U64;
    }
    if (bits <= 128) {
        return ALIAS_BACKEND_U128;
    }
    return ALIAS_BACKEND_GMP;
}

/* log2 of a positive integer, in double precision. */
static double log2_mpz(const mpz_t x) {
    long e;
    double m = mpz_get_d_2exp(&e, x);
    return log2(m) + (double)e;
}

void alias_analyze(const alias_dist *d, alias_report *r) {
    memset(r, 0, sizeof(*r));
    r->n = d->n;
    r->total_bits = mpz_sizeinbase(d->total, 2);
    r->gcd_bits = mpz_sizeinbase(d->gcd, 2);

    mpz_t cs, t;
    mpz_inits(cs, t, NULL);
    mpz_fdiv_q_ui(cs, d->total, d->n);
    r->cs_bits = mpz_sgn(cs) ? mpz_sizeinbase(cs, 2) : 0;

    r->scale = d->scale;
    mpz_mul(t, d->total, d->gcd);
    mpz_divexact_ui(t, t, d->scale);                /* sum of the original weights */
    r->declared_total_matches = (mpz_cmp(t, d->declared_total) == 0);

    double log2_w = log2_mpz(d->total);
    double h = 0;
    for (uint32_t i = 0; i < d->n; i++) {
        const mpz_t *wi = &d->weights[i];
        if (mpz_sgn(*wi) == 0) {
            r->n_zero++;
            continue;
        }
        size_t b = mpz_sizeinbase(*wi, 2);
        if (b > r->max_weight_bits) {
            r->max_weight_bits = b;
        }
        double l = log2_mpz(*wi) - log2_w;          /* log2 p_i */
        h -= exp2(l) * l;
        /* n w_i / w < 2^-53  <=>  n w_i 2^53 < w */
        mpz_mul_ui(t, *wi, d->n);
        mpz_mul_2exp(t, t, 53);
        if (mpz_cmp(t, d->total) < 0) {
            r->below_double++;
        }
    }
    r->entropy = h;
    r->backend = alias_min_backend(d);
    mpz_clears(cs, t, NULL);
}

const char *alias_backend_name(alias_backend b) {
    switch (b) {
    case ALIAS_BACKEND_U64:  return "u64";
    case ALIAS_BACKEND_U128: return "u128";
    case ALIAS_BACKEND_GMP:  return "gmp";
    default:                 return "auto";
    }
}

const char *alias_mode_name(alias_mode m) {
    return m == ALIAS_SPEED ? "speed" : "alea";
}
