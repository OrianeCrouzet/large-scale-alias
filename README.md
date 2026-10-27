# Integer-Alias Sampler

Exact sampling from a discrete distribution given by **integer weights of
arbitrary size**, with the alias method, in constant time per sample.

No floating-point number is used to choose the sampled object: the output
follows the distribution `w_i / Σ w_j` exactly, even when the weights have
hundreds of bits — a regime where double-precision alias samplers silently
never draw some objects.

This is the implementation accompanying the paper *Exact Integer Alias Method:
Robust and Error-Free Large-Scale Sampling* (O. Crouzet, A. Genitrini,
M. Naima).

## Features

* **Two sampling modes**, same exact distribution:
  * `speed`: fastest; one 64-bit random word per sample (Lemire's uniform,
    Bernoulli on 32 precomputed digits);
  * `alea`: fewest random bits (Lumbroso's optimal uniform and Bernoulli),
    `log2 n + O(1)` bits per sample.
* **Automatic choice of the arithmetic**: the distribution is analyzed before
  the construction, and the table is built with `uint64_t`, `unsigned __int128`
  or GMP integers, whichever is the smallest exact representation.
* Linear-time construction (Algorithm 1 of the paper), 16-byte cells.
* Cryptographic random bits (`arc4random` / `getrandom`), produced by a
  background thread.

## Build

Requirements: a C11 compiler with `unsigned __int128` (GCC or Clang),
[GMP](https://gmplib.org/), POSIX threads.

```sh
make            # libalias.a and alias_sampler
make test       # runs the tests on examples/
```

## Usage

```sh
./alias_sampler -m speed -n 5 examples/small.dist
./alias_sampler -m alea -n 1000000 -o histogram examples/dirichlet_20484.dist > counts.txt
./alias_sampler -a examples/u128.dist          # analysis only
```

Example of analysis (on stderr):

```
distribution
  objects                 20484 (0 of weight 0)
  total weight            169 bits
  largest weight          164 bits
  cell size cs            154 bits
  entropy                 7.08217 bits (max log2 n = 14.32221)
  below double resolution 17525 objects (n w_i / w < 2^-53: never drawn by a double-precision Alias sampler)
  backend                 gmp (smallest exact representation)
```

Options:

| Option | Values | Default |
|---|---|---|
| `-m, --mode` | `speed`, `alea` | `speed` |
| `-b, --backend` | `auto`, `u64`, `u128`, `gmp` | `auto` |
| `-n, --samples` | number of samples | 10 |
| `-o, --output` | `samples` (one index per line), `histogram`, `none` | `samples` |
| `-a, --analyze-only` | | |
| `-q, --quiet` | | |

## Input format

A `.dist` file holds whitespace-separated decimal integers: a declared total
weight (informative), then `n` followed by the `n` weights, then optionally the
entropy:

```
39
10 3 1 4 1 5 9 2 6 5 3
3.06058
```

Objects are numbered from 0 in the order of the file. Weights may be 0.
`tools/make_examples.py` regenerates the small examples of `examples/`.

## Library

```c
#include "alias.h"

alias_dist d;
char err[256];
alias_dist_read("weights.dist", &d, err, sizeof err);
alias_report r;
alias_analyze(&d, &r);
alias_table *t = alias_build(&d, ALIAS_SPEED, r.backend);
alias_rng_start();
uint32_t x = alias_sample(t);
alias_rng_stop();
alias_table_free(t);
alias_dist_free(&d);
```

Link with `libalias.a -lgmp -lpthread -lm`.

## Documentation

[TECHNICAL.md](TECHNICAL.md) describes the algorithms, the data structures,
the backend selection and the tests in detail.

## Reference

If you use this code, please cite the paper describing the method:

> Oriane Crouzet, Antoine Genitrini, Mehdi Naima.
> *Exact Integer Alias Method: Robust and Error-Free Large-Scale Sampling*.
> 2026. Link to the archived version: to be added.

```bibtex
@misc{CGN26,
  author = {Crouzet, Oriane and Genitrini, Antoine and Naima, Mehdi},
  title  = {Exact Integer Alias Method: Robust and Error-Free Large-Scale Sampling},
  year   = {2026},
  note   = {Preprint}
}
```

## License

Copyright (C) 2026 Oriane Crouzet, Antoine Genitrini, Mehdi Naima.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version. See [LICENSE](LICENSE).
