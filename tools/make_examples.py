#!/usr/bin/env python3
# Copyright (C) 2026 Oriane Crouzet, Antoine Genitrini, Mehdi Naima
# SPDX-License-Identifier: GPL-3.0-or-later
"""Writes the example distributions of examples/ (deterministic).

Format of a .dist file: the total weight, then n followed by the n weights,
then the entropy of the distribution.
"""

import math
import os
import random
import sys


def entropy(weights):
    total = sum(weights)
    h = 0.0
    for w in weights:
        if w:
            p = w / total
            h -= p * math.log2(p)
    return h


def write(path, weights):
    with open(path, 'w') as f:
        f.write('%d\n' % sum(weights))
        f.write('%d %s\n' % (len(weights), ' '.join(str(w) for w in weights)))
        f.write('%.5f\n' % entropy(weights))
    print(path, 'n=%d' % len(weights), 'total=%d bits' % sum(weights).bit_length())


def main(out):
    rng = random.Random(2027)
    os.makedirs(out, exist_ok=True)
    # Small distribution, fits in 64 bits.
    write(os.path.join(out, 'small.dist'), [3, 1, 4, 1, 5, 9, 2, 6, 5, 3])
    # 1000 weights below 2^50: total below 2^64.
    write(os.path.join(out, 'u64.dist'), [rng.randrange(1, 1 << 50) for _ in range(1000)])
    # 1000 weights below 2^100, skewed: total below 2^128.
    write(os.path.join(out, 'u128.dist'),
          [rng.randrange(1, 1 << rng.randrange(1, 100)) for _ in range(1000)])
    # Weights 0 and a total smaller than n (the weights are multiplied by n).
    write(os.path.join(out, 'zeros.dist'), [0, 1, 0, 2, 0, 0])


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'examples')
