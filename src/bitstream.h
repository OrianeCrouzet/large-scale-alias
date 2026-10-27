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
 * bitstream.h -- stream of cryptographic random bits.
 *
 * A background thread fills a ring buffer of 64-bit words with a
 * cryptographically secure generator (arc4random_buf, or getrandom on
 * Linux). The sampling thread reads bits from the current word with
 * inline functions; it only synchronizes with the producer when it has
 * consumed every word it knows to be available.
 *
 * Bits are consumed from the least significant bit of each word.
 * There is a single stream per process (single consumer thread).
 */

#ifndef ALIAS_BITSTREAM_H
#define ALIAS_BITSTREAM_H

#include <stdint.h>

#define BS_WORDS (1u << 15)     /* 32768 words of 64 bits = 256 KiB */
#define BS_MASK (BS_WORDS - 1)
#define BS_CHUNK 512            /* words written per call to the generator */

typedef struct {
    _Alignas(128) uint64_t tail;     /* written by the producer (release) */
    _Alignas(128) uint64_t head;     /* published by the consumer (release) */
    _Alignas(128) uint64_t c_head;   /* consumer: index of the current word */
    uint64_t c_tail;                 /* consumer: last tail seen */
    uint32_t c_off;                  /* consumer: bits already read in the current word */
    uint32_t c_publish;              /* consumer: words consumed since the last publication */
    _Alignas(128) volatile int running;
    _Alignas(128) uint64_t data[BS_WORDS];
} bs_queue_t;

extern bs_queue_t bs_q;

void bs_start(void);           /* starts the producer thread */
void bs_stop(void);            /* stops and joins it */
void bs_wait_slow(void);       /* waits for new words (out of line) */

/* Number of random bits consumed since bs_start. */
static inline uint64_t bs_bits_consumed(void) {
    return 64 * bs_q.c_head + bs_q.c_off;
}

static inline void bs_ready(void) {
    if (__builtin_expect(bs_q.c_head == bs_q.c_tail, 0)) {
        bs_wait_slow();
    }
}

/* Current word shifted so that bit 0 is the next unread bit; *avail = unread bits. */
static inline uint64_t bs_peek(uint32_t *avail) {
    bs_ready();
    *avail = 64 - bs_q.c_off;
    return bs_q.data[bs_q.c_head & BS_MASK] >> bs_q.c_off;
}

/* Consumes nb <= avail bits of the current word. */
static inline void bs_skip(uint32_t nb) {
    bs_q.c_off += nb;
    if (bs_q.c_off == 64) {
        bs_q.c_off = 0;
        bs_q.c_head++;
        if (++bs_q.c_publish == 64) {
            bs_q.c_publish = 0;
            __atomic_store_n(&bs_q.head, bs_q.c_head, __ATOMIC_RELEASE);
        }
    }
}

/* nb random bits, 1 <= nb <= 64. */
static inline uint64_t bs_bits(uint32_t nb) {
    uint32_t avail;
    uint64_t w = bs_peek(&avail);
    if (nb <= avail) {
        bs_skip(nb);
        return (nb == 64) ? w : (w & ((1ull << nb) - 1));
    }
    uint64_t lo = w;                 /* avail < 64 bits */
    bs_skip(avail);
    uint64_t hi = bs_bits(nb - avail);
    return lo | (hi << avail);
}

static inline uint32_t bs_bit(void) {
    uint32_t avail;
    uint64_t w = bs_peek(&avail);
    bs_skip(1);
    return (uint32_t)(w & 1);
}

/* Next 64-bit word; fast when the stream is aligned on a word boundary. */
static inline uint64_t bs_word(void) {
    uint32_t avail;
    uint64_t w = bs_peek(&avail);
    if (__builtin_expect(avail == 64, 1)) {
        bs_skip(64);
        return w;
    }
    return bs_bits(64);
}

/* Reads coins until the first 1, at most cap coins.
 * Returns the position (1..cap) of the first 1, or 0 if the cap coins are all 0. */
static inline uint32_t bs_coins_until_one(uint32_t cap) {
    uint32_t consumed = 0;
    for (;;) {
        uint32_t avail;
        uint64_t w = bs_peek(&avail);
        uint32_t look = avail;
        if (look > cap - consumed) {
            look = cap - consumed;
        }
        if (w != 0) {
            uint32_t tz = (uint32_t)__builtin_ctzll(w);
            if (tz < look) {
                bs_skip(tz + 1);
                return consumed + tz + 1;
            }
        }
        bs_skip(look);
        consumed += look;
        if (consumed == cap) {
            return 0;
        }
    }
}

#endif
