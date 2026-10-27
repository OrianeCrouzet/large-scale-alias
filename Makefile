# Copyright (C) 2026 Oriane Crouzet, Antoine Genitrini, Mehdi Naima
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Integer-Alias Sampler
#
#   make            builds the library (libalias.a) and the command-line sampler
#   make test       builds and runs the tests
#   make bench      compares the backends on random distributions
#   make clean

CC      ?= cc
CFLAGS  ?= -O3
CFLAGS  += -std=c11 -Wall -Wextra -Iinclude -Isrc
GMP_CFLAGS := $(shell pkg-config --cflags gmp 2>/dev/null)
GMP_LIBS   := $(shell pkg-config --libs gmp 2>/dev/null || echo -lgmp)
CFLAGS  += $(GMP_CFLAGS)
LDLIBS  += $(GMP_LIBS) -lpthread -lm

LIB_SRC := src/bitstream.c src/dist.c src/build.c src/sample.c
LIB_OBJ := $(LIB_SRC:.c=.o)

all: libalias.a alias_sampler

libalias.a: $(LIB_OBJ)
	ar rcs $@ $^

src/build.o: src/build.c src/build_native.h src/alias_internal.h src/bitstream.h include/alias.h
src/sample.o: src/sample.c src/alias_internal.h src/bitstream.h include/alias.h
src/bitstream.o: src/bitstream.c src/bitstream.h
src/dist.o: src/dist.c include/alias.h

alias_sampler: src/main.c libalias.a
	$(CC) $(CFLAGS) -o $@ src/main.c libalias.a $(LDLIBS)

tests/test_alias: tests/test_alias.c libalias.a
	$(CC) $(CFLAGS) -o $@ tests/test_alias.c libalias.a $(LDLIBS)

test: tests/test_alias alias_sampler
	./tests/test_alias examples

tools/bench_backends: tools/bench_backends.c libalias.a
	$(CC) $(CFLAGS) -o $@ tools/bench_backends.c libalias.a $(LDLIBS)

# Comparison of the backends on random distributions (about 2 minutes, 100 MB of files).
bench: tools/bench_backends
	python3 tools/make_bench_dists.py bench_dists
	for f in bench_dists/*.dist; do echo "== $$f"; ./tools/bench_backends $$f; done

clean:
	rm -f $(LIB_OBJ) libalias.a alias_sampler tests/test_alias tools/bench_backends

.PHONY: all test bench clean
