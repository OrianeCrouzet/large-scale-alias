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
 * bench_backends.c -- compares the backends on one distribution.
 *
 *   bench_backends FILE.dist [SAMPLES]
 *
 * For every backend at least as large as the minimal one, and for both modes:
 * construction time, table memory, sampling time and random bits per sample.
 * Times are the best of 3 runs (CLOCK_MONOTONIC_RAW).
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "alias.h"

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC_RAW, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s FILE.dist [SAMPLES]\n", argv[0]);
        return EXIT_FAILURE;
    }
    long samples = argc > 2 ? atol(argv[2]) : 10000000;
    alias_dist d;
    char err[256];
    if (alias_dist_read(argv[1], &d, err, sizeof err) != 0) {
        fprintf(stderr, "%s\n", err);
        return EXIT_FAILURE;
    }
    alias_rng_start();
    for (alias_backend be = alias_min_backend(&d); be <= ALIAS_BACKEND_GMP; be++) {
        double build[2] = {1e9, 1e9}, sample[2] = {1e9, 1e9}, bits[2] = {0, 0};
        size_t bytes = 0;
        for (int m = 0; m < 2; m++) {
            alias_mode mode = m ? ALIAS_SPEED : ALIAS_ALEA;
            for (int rep = 0; rep < 3; rep++) {
                double t0 = now();
                alias_table *t = alias_build(&d, mode, be);
                double tb = now() - t0;
                if (t == NULL) {
                    return EXIT_FAILURE;
                }
                if (tb < build[m]) {
                    build[m] = tb;
                }
                alias_table_info info;
                alias_info(t, &info);
                bytes = info.bytes;
                uint64_t b0 = alias_rng_bits();
                volatile uint32_t sink = 0;
                t0 = now();
                for (long k = 0; k < samples; k++) {
                    sink += alias_sample(t);
                }
                double ts = (now() - t0) / (double)samples;
                if (ts < sample[m]) {
                    sample[m] = ts;
                }
                bits[m] = (double)(alias_rng_bits() - b0) / (double)samples;
                alias_table_free(t);
            }
        }
        printf("%-4s  build %9.3f ms  memory %10.1f KB  alea %6.2f ns %6.2f bits"
               "  speed %6.2f ns %6.2f bits\n",
               alias_backend_name(be), 1e3 * build[0], bytes / 1024.0,
               1e9 * sample[0], bits[0], 1e9 * sample[1], bits[1]);
    }
    alias_rng_stop();
    alias_dist_free(&d);
    return EXIT_SUCCESS;
}
