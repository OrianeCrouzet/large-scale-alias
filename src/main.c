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
 * main.c -- command-line sampler.
 *
 *   alias_sampler [options] FILE.dist
 *
 * Reads the distribution, analyzes it (size of the weights, entropy, backend),
 * builds the alias table, then draws the requested number of objects.
 * The analysis and the statistics are written on stderr, the samples on stdout.
 */

#include <getopt.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "alias.h"

typedef enum { OUT_SAMPLES, OUT_HISTOGRAM, OUT_NONE } output_kind;

static void usage(FILE *f, const char *prog) {
    fprintf(f,
        "usage: %s [options] FILE.dist\n"
        "\n"
        "Exact sampling from the distribution given by the integer weights of FILE.dist.\n"
        "\n"
        "options:\n"
        "  -m, --mode MODE        speed (default): fastest sampling, 64 random bits per draw\n"
        "                         alea: fewest random bits (Lumbroso's optimal uniform and Bernoulli)\n"
        "  -b, --backend BACKEND  auto (default), u64, u128 or gmp\n"
        "  -n, --samples N        number of objects to draw (default 10)\n"
        "  -o, --output KIND      samples (default, one index per line), histogram, or none\n"
        "  -a, --analyze-only     analyze the distribution and stop\n"
        "  -q, --quiet            no analysis or statistics on stderr\n"
        "  -h, --help             this message\n", prog);
}

static double now_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static void print_report(const alias_report *r, const alias_dist *d) {
    fprintf(stderr, "distribution\n");
    fprintf(stderr, "  objects                 %" PRIu32 " (%" PRIu32 " of weight 0)\n", r->n, r->n_zero);
    fprintf(stderr, "  total weight            %zu bits", r->total_bits);
    if (r->gcd_bits > 1) {
        fprintf(stderr, " (after division by a gcd of %zu bits)", r->gcd_bits);
    }
    if (r->scale > 1) {
        fprintf(stderr, " (weights multiplied by n, as their sum was < n)");
    }
    fprintf(stderr, "\n  largest weight          %zu bits\n", r->max_weight_bits);
    fprintf(stderr, "  cell size cs            %zu bits\n", r->cs_bits);
    fprintf(stderr, "  entropy                 %.5f bits (max log2 n = %.5f)\n", r->entropy, log2((double)r->n));
    if (d->has_entropy && fabs(d->declared_entropy - r->entropy) > 1e-4) {
        fprintf(stderr, "  note: the file declares an entropy of %.5f\n", d->declared_entropy);
    }
    if (!r->declared_total_matches) {
        fprintf(stderr, "  note: the first line of the file differs from the sum of the weights;"
                        " the weights are used\n");
    }
    fprintf(stderr, "  below double resolution %" PRIu32 " objects (n w_i / w < 2^-53: never drawn"
                    " by a double-precision Alias sampler)\n", r->below_double);
    fprintf(stderr, "  backend                 %s (smallest exact representation)\n",
            alias_backend_name(r->backend));
}

int main(int argc, char **argv) {
    alias_mode mode = ALIAS_SPEED;
    alias_backend backend = ALIAS_BACKEND_AUTO;
    unsigned long long samples = 10;
    output_kind output = OUT_SAMPLES;
    int analyze_only = 0, quiet = 0;

    static const struct option opts[] = {
        {"mode", required_argument, NULL, 'm'},
        {"backend", required_argument, NULL, 'b'},
        {"samples", required_argument, NULL, 'n'},
        {"output", required_argument, NULL, 'o'},
        {"analyze-only", no_argument, NULL, 'a'},
        {"quiet", no_argument, NULL, 'q'},
        {"help", no_argument, NULL, 'h'},
        {NULL, 0, NULL, 0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "m:b:n:o:aqh", opts, NULL)) != -1) {
        switch (c) {
        case 'm':
            if (strcmp(optarg, "speed") == 0) {
                mode = ALIAS_SPEED;
            } else if (strcmp(optarg, "alea") == 0) {
                mode = ALIAS_ALEA;
            } else {
                fprintf(stderr, "unknown mode '%s' (speed or alea)\n", optarg);
                return EXIT_FAILURE;
            }
            break;
        case 'b':
            if (strcmp(optarg, "auto") == 0) {
                backend = ALIAS_BACKEND_AUTO;
            } else if (strcmp(optarg, "u64") == 0) {
                backend = ALIAS_BACKEND_U64;
            } else if (strcmp(optarg, "u128") == 0) {
                backend = ALIAS_BACKEND_U128;
            } else if (strcmp(optarg, "gmp") == 0) {
                backend = ALIAS_BACKEND_GMP;
            } else {
                fprintf(stderr, "unknown backend '%s' (auto, u64, u128 or gmp)\n", optarg);
                return EXIT_FAILURE;
            }
            break;
        case 'n': {
            char *end;
            samples = strtoull(optarg, &end, 10);
            if (*end != '\0') {
                fprintf(stderr, "invalid number of samples '%s'\n", optarg);
                return EXIT_FAILURE;
            }
            break;
        }
        case 'o':
            if (strcmp(optarg, "samples") == 0) {
                output = OUT_SAMPLES;
            } else if (strcmp(optarg, "histogram") == 0) {
                output = OUT_HISTOGRAM;
            } else if (strcmp(optarg, "none") == 0) {
                output = OUT_NONE;
            } else {
                fprintf(stderr, "unknown output '%s' (samples, histogram or none)\n", optarg);
                return EXIT_FAILURE;
            }
            break;
        case 'a': analyze_only = 1; break;
        case 'q': quiet = 1; break;
        case 'h': usage(stdout, argv[0]); return EXIT_SUCCESS;
        default:  usage(stderr, argv[0]); return EXIT_FAILURE;
        }
    }
    if (optind != argc - 1) {
        usage(stderr, argv[0]);
        return EXIT_FAILURE;
    }

    /* 1. Reading and analysis of the distribution. */
    alias_dist d;
    char err[512];
    if (alias_dist_read(argv[optind], &d, err, sizeof err) != 0) {
        fprintf(stderr, "%s\n", err);
        return EXIT_FAILURE;
    }
    alias_report report;
    alias_analyze(&d, &report);
    if (!quiet) {
        print_report(&report, &d);
    }
    if (analyze_only) {
        alias_dist_free(&d);
        return EXIT_SUCCESS;
    }
    if (backend == ALIAS_BACKEND_AUTO) {
        backend = report.backend;
    }

    /* 2. Preprocessing. */
    double t0 = now_seconds();
    alias_table *t = alias_build(&d, mode, backend);
    double t_build = now_seconds() - t0;
    if (t == NULL) {
        alias_dist_free(&d);
        return EXIT_FAILURE;
    }
    alias_table_info info;
    alias_info(t, &info);
    if (!quiet) {
        fprintf(stderr, "table\n");
        fprintf(stderr, "  mode, backend           %s, %s\n", alias_mode_name(mode), alias_backend_name(backend));
        fprintf(stderr, "  cells                   %" PRIu32 " (%" PRIu32 " with an alias%s)\n",
                info.q, info.aliased, info.has_virtual ? ", one virtual object" : "");
        fprintf(stderr, "  rejection probability   %.3g per attempt\n", info.reject_prob);
        fprintf(stderr, "  memory                  %zu bytes\n", info.bytes);
        fprintf(stderr, "  preprocessing time      %.3f ms\n", 1e3 * t_build);
    }

    /* 3. Sampling. */
    uint64_t *histogram = NULL;
    if (output == OUT_HISTOGRAM) {
        histogram = calloc(d.n, sizeof(uint64_t));
    }
    alias_rng_start();
    uint64_t sink = 0;
    t0 = now_seconds();
    for (unsigned long long k = 0; k < samples; k++) {
        uint32_t x = alias_sample(t);
        switch (output) {
        case OUT_SAMPLES:   printf("%" PRIu32 "\n", x); break;
        case OUT_HISTOGRAM: histogram[x]++; break;
        default:            sink += x; break;
        }
    }
    double t_sample = now_seconds() - t0;
    uint64_t bits = alias_rng_bits();
    alias_rng_stop();

    if (output == OUT_HISTOGRAM) {
        for (uint32_t i = 0; i < d.n; i++) {
            printf("%" PRIu32 " %" PRIu64 "\n", i, histogram[i]);
        }
        free(histogram);
    }
    if (!quiet && samples > 0) {
        fprintf(stderr, "sampling\n");
        fprintf(stderr, "  samples                 %llu\n", samples);
        fprintf(stderr, "  time per sample         %.2f ns%s\n", 1e9 * t_sample / (double)samples,
                output == OUT_SAMPLES ? " (including the output)" : "");
        fprintf(stderr, "  random bits per sample  %.3f\n", (double)bits / (double)samples);
    }
    if (sink == 1) {
        fprintf(stderr, " ");                /* keeps the loop from being optimized away */
    }
    alias_table_free(t);
    alias_dist_free(&d);
    return EXIT_SUCCESS;
}
