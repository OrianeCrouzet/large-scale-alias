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
 * bitstream.c -- producer thread of the random-bit stream (see bitstream.h).
 */

#include <pthread.h>
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#if defined(__linux__)
#include <sys/random.h>
#endif

#include "bitstream.h"

bs_queue_t bs_q;
static pthread_t bs_thread;

/* Fills buf with len cryptographically secure random bytes. */
static void fill_random(void *buf, size_t len) {
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    arc4random_buf(buf, len);
#elif defined(__linux__)
    unsigned char *p = buf;
    while (len > 0) {
        ssize_t got = getrandom(p, len, 0);
        if (got < 0) {
            perror("getrandom");
            exit(EXIT_FAILURE);
        }
        p += got;
        len -= (size_t)got;
    }
#else
#error "no cryptographic random generator available on this platform"
#endif
}

static void *producer(void *arg) {
    (void)arg;
    while (bs_q.running) {
        uint64_t t = bs_q.tail;
        uint64_t h = __atomic_load_n(&bs_q.head, __ATOMIC_ACQUIRE);
        if (BS_WORDS - (t - h) >= BS_CHUNK) {
            /* BS_WORDS is a multiple of BS_CHUNK: a chunk never wraps around. */
            fill_random(&bs_q.data[t & BS_MASK], BS_CHUNK * sizeof(uint64_t));
            __atomic_store_n(&bs_q.tail, t + BS_CHUNK, __ATOMIC_RELEASE);
        } else {
            /* Full buffer (about two million bits): sleep instead of spinning. */
            usleep(50);
        }
    }
    return NULL;
}

void bs_start(void) {
    memset(&bs_q, 0, sizeof(bs_q));
    bs_q.running = 1;
    if (pthread_create(&bs_thread, NULL, producer, NULL) != 0) {
        perror("bs_start: pthread_create");
        exit(EXIT_FAILURE);
    }
}

void bs_stop(void) {
    bs_q.running = 0;
    pthread_join(bs_thread, NULL);
}

/* Called when every known word has been consumed. */
void bs_wait_slow(void) {
    __atomic_store_n(&bs_q.head, bs_q.c_head, __ATOMIC_RELEASE);
    bs_q.c_publish = 0;
    for (;;) {
        bs_q.c_tail = __atomic_load_n(&bs_q.tail, __ATOMIC_ACQUIRE);
        if (bs_q.c_head != bs_q.c_tail) {
            return;
        }
        sched_yield();
    }
}
