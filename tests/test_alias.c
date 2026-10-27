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
 * test_alias.c -- tests of the Integer-Alias Sampler.
 *
 *   test_alias DIR    runs every test on the .dist files of DIR
 *
 * 1. The backends build the same table (cells, thresholds, virtual object).
 * 2. Chi-square goodness of fit of the samples, for each mode and backend.
 *    Objects of weight 0 and the virtual object are never returned.
 * 3. The slow paths (probability <= 2^-31 in normal use) are called directly
 *    and their frequency of 1 is compared with (T_i 2^k mod cs) / cs.
 */

#include <dirent.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "alias.h"
#include "alias_internal.h"

static int failures = 0;

#define CHECK(cond, ...) do {                                   \
    if (!(cond)) {                                              \
        failures++;                                             \
        printf("  FAIL: ");                                     \
        printf(__VA_ARGS__);                                    \
        printf("\n");                                           \
    }                                                           \
} while (0)

/* 1. Same table for every backend at least as large as the minimal one. */
static void test_backends_agree(const alias_dist *d, alias_mode mode) {
    alias_backend min = alias_min_backend(d);
    alias_table *ref = alias_build(d, mode, ALIAS_BACKEND_GMP);
    CHECK(ref != NULL, "GMP build failed");
    if (ref == NULL) {
        return;
    }
    mpz_t a, b;
    mpz_inits(a, b, NULL);
    for (alias_backend be = min; be < ALIAS_BACKEND_GMP; be++) {
        alias_table *t = alias_build(d, mode, be);
        CHECK(t != NULL, "%s build failed", alias_backend_name(be));
        if (t == NULL) {
            continue;
        }
        int same = t->q == ref->q && t->m == ref->m && t->virtual_obj == ref->virtual_obj
                   && mpz_cmp(t->cs_mpz, ref->cs_mpz) == 0;
        for (uint32_t i = 0; same && i < t->q; i++) {
            const alias_cell *x = &t->cells[i], *y = &ref->cells[i];
            same = x->primary == y->primary && x->alias == y->alias && x->digits == y->digits;
            if (same && i < t->m) {
                /* k = 0: the threshold itself (T_i < cs). */
                alias_threshold_remainder(t, i, 0, a);
                alias_threshold_remainder(ref, i, 0, b);
                same = mpz_cmp(a, b) == 0;
            }
        }
        CHECK(same, "%s and gmp tables differ (%s)", alias_backend_name(be), alias_mode_name(mode));
        printf("  tables %-4s = gmp (%s): %s\n", alias_backend_name(be), alias_mode_name(mode),
               same ? "identical" : "DIFFERENT");
        alias_table_free(t);
    }
    mpz_clears(a, b, NULL);
    alias_table_free(ref);
}

/* Probability of object i, in double precision. */
static double probability(const alias_dist *d, uint32_t i) {
    if (mpz_sgn(d->weights[i]) == 0) {
        return 0.0;
    }
    long ei, ew;
    double mi = mpz_get_d_2exp(&ei, d->weights[i]);
    double mw = mpz_get_d_2exp(&ew, d->total);
    return ldexp(mi / mw, (int)(ei - ew));
}

/* 2. Chi-square test, objects of expected count < 5 pooled into one class.
 * The p-value is approximated by the Wilson-Hilferty transform; the test fails
 * only for z > 4.5 (p < 4e-6), so that it is stable across runs. */
static void test_chi2(const alias_dist *d, alias_mode mode, alias_backend be, uint64_t samples) {
    alias_table *t = alias_build(d, mode, be);
    CHECK(t != NULL, "build failed");
    if (t == NULL) {
        return;
    }
    uint64_t *count = calloc(d->n, sizeof(uint64_t));
    int bad_index = 0;
    uint64_t bits0 = alias_rng_bits();
    for (uint64_t k = 0; k < samples; k++) {
        uint32_t x = alias_sample(t);
        if (x >= d->n) {
            bad_index = 1;
            break;
        }
        count[x]++;
    }
    double bits = (double)(alias_rng_bits() - bits0) / (double)samples;
    CHECK(!bad_index, "index >= n returned (virtual object?)");

    double chi2 = 0, pooled_obs = 0, pooled_exp = 0;
    unsigned classes = 0;
    int zero_drawn = 0;
    for (uint32_t i = 0; i < d->n; i++) {
        double p = probability(d, i);
        if (p == 0.0) {
            zero_drawn |= (count[i] != 0);
            continue;
        }
        double e = p * (double)samples;
        if (e < 5) {
            pooled_obs += (double)count[i];
            pooled_exp += e;
            continue;
        }
        chi2 += ((double)count[i] - e) * ((double)count[i] - e) / e;
        classes++;
    }
    if (pooled_exp > 0) {
        chi2 += (pooled_obs - pooled_exp) * (pooled_obs - pooled_exp) / pooled_exp;
        classes++;
    }
    CHECK(!zero_drawn, "an object of weight 0 was drawn");
    double df = classes > 1 ? classes - 1 : 1;
    double z = (cbrt(chi2 / df) - (1 - 2 / (9 * df))) / sqrt(2 / (9 * df));
    CHECK(z < 4.5, "chi-square too large: chi2 = %.1f, df = %.0f, z = %.2f", chi2, df, z);
    printf("  chi2 %-5s %-4s  %8.1f (df %6.0f, z %+5.2f)  %6.2f bits/sample\n",
           alias_mode_name(mode), alias_backend_name(be), chi2, df, z, bits);
    free(count);
    alias_table_free(t);
}

/* 3. Slow paths: frequency of 1 against (T_i 2^k mod cs) / cs. */
static void test_slow_paths(const alias_dist *d, alias_backend be) {
    for (int mi = 0; mi < 2; mi++) {
        alias_mode mode = mi ? ALIAS_SPEED : ALIAS_ALEA;
        unsigned k = mi ? SPEED_DIGITS : ALEA_DIGITS;
        alias_table *t = alias_build(d, mode, be);
        if (t == NULL || t->m == 0) {
            alias_table_free(t);
            continue;
        }
        mpz_t rem;
        mpz_init(rem);
        double worst = 0;
        const int cells = 5, calls = 100000;
        for (int c = 0; c < cells; c++) {
            uint32_t i = (uint32_t)((uint64_t)c * t->m / cells);
            alias_threshold_remainder(t, i, k, rem);
            double p = mpz_get_d(rem) / mpz_get_d(t->cs_mpz);
            long ones = 0;
            for (int j = 0; j < calls; j++) {
                ones += mi ? alias_slow_tie(t, i) : alias_slow_tail(t, i);
            }
            double f = (double)ones / calls;
            double sd = sqrt(p * (1 - p) / calls);
            double z = sd > 0 ? fabs(f - p) / sd : (f == p ? 0 : INFINITY);
            if (z > worst) {
                worst = z;
            }
        }
        CHECK(worst < 5, "slow path %s/%s: |z| = %.2f", alias_mode_name(mode), alias_backend_name(be), worst);
        printf("  slow path %-5s %-4s  max |z| = %.2f over %d cells\n",
               alias_mode_name(mode), alias_backend_name(be), worst, cells);
        mpz_clear(rem);
        alias_table_free(t);
    }
}

int main(int argc, char **argv) {
    const char *dir = argc > 1 ? argv[1] : "examples";
    DIR *dp = opendir(dir);
    if (dp == NULL) {
        perror(dir);
        return EXIT_FAILURE;
    }
    char *names[256];
    int nfiles = 0;
    struct dirent *e;
    while ((e = readdir(dp)) != NULL && nfiles < 256) {
        size_t len = strlen(e->d_name);
        if (len > 5 && strcmp(e->d_name + len - 5, ".dist") == 0) {
            names[nfiles++] = strdup(e->d_name);
        }
    }
    closedir(dp);

    alias_rng_start();
    for (int f = 0; f < nfiles; f++) {
        char path[1024];
        snprintf(path, sizeof path, "%s/%s", dir, names[f]);
        alias_dist d;
        char err[512];
        if (alias_dist_read(path, &d, err, sizeof err) != 0) {
            CHECK(0, "%s", err);
            continue;
        }
        alias_backend min = alias_min_backend(&d);
        printf("%s (n = %u, %zu-bit total, backend %s)\n", names[f], d.n,
               mpz_sizeinbase(d.total, 2), alias_backend_name(min));
        test_backends_agree(&d, ALIAS_ALEA);
        test_backends_agree(&d, ALIAS_SPEED);
        uint64_t samples = d.n > 5000 ? 5000000 : 2000000;
        for (alias_backend be = min; be <= ALIAS_BACKEND_GMP; be++) {
            test_chi2(&d, ALIAS_ALEA, be, samples);
            test_chi2(&d, ALIAS_SPEED, be, samples);
        }
        test_slow_paths(&d, min);
        if (min != ALIAS_BACKEND_GMP) {
            test_slow_paths(&d, ALIAS_BACKEND_GMP);
        }
        alias_dist_free(&d);
        free(names[f]);
    }
    alias_rng_stop();

    printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "all tests passed", failures,
           failures == 1 ? "" : "s");
    return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
