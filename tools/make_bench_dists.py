#!/usr/bin/env python3
# Copyright (C) 2026 Oriane Crouzet, Antoine Genitrini, Mehdi Naima
# SPDX-License-Identifier: GPL-3.0-or-later
"""Writes the random distributions used by `make bench` (deterministic).

For n in {10^3, 10^5, 10^6} and weights of 30 or 90 bits:
  unif_*: weights uniform on [1, 2^bits);
  skew_*: weights uniform on [1, 2^b), with b uniform on [1, bits) for each object.
30-bit weights fit every backend (total < 2^64), 90-bit weights need u128 or gmp.
"""

import os
import random
import sys


def write(path, weights):
    with open(path, 'w') as f:
        f.write('%d\n%d %s\n' % (sum(weights), len(weights), ' '.join(map(str, weights))))


def main(out):
    rng = random.Random(7)
    os.makedirs(out, exist_ok=True)
    for n in (1000, 100000, 1000000):
        for bits in (30, 90):
            write(os.path.join(out, 'unif_n%d_b%d.dist' % (n, bits)),
                  [rng.randrange(1, 1 << bits) for _ in range(n)])
            write(os.path.join(out, 'skew_n%d_b%d.dist' % (n, bits)),
                  [rng.randrange(1, 1 << rng.randrange(1, bits)) for _ in range(n)])


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'bench_dists')
