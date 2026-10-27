# Integer-Alias Sampler — technical description

This document describes how the sampler works, step by step, and the data
structures behind it. The algorithms are those of the paper *Exact Integer
Alias Method: Robust and Error-Free Large-Scale Sampling* (Crouzet, Genitrini,
Naima); this document focuses on their implementation.

## Contents

1. [Overview](#1-overview)
2. [Input: the `.dist` file](#2-input-the-dist-file)
3. [Analysis and backend selection](#3-analysis-and-backend-selection)
4. [Construction of the alias table](#4-construction-of-the-alias-table)
5. [Data structures](#5-data-structures)
6. [Sampling step](#6-sampling-step)
7. [Random-bit stream](#7-random-bit-stream)
8. [Exactness and tests](#8-exactness-and-tests)
9. [Performance](#9-performance)
10. [Limits](#10-limits)
11. [Interfaces](#11-interfaces)

---

## 1. Overview

The program samples objects `0 .. n-1` with probability `w_i / w`, where the
`w_i` are non-negative integers of arbitrary size and `w = Σ w_i`. Sampling is
**exact**: no floating-point number is ever used to decide which object is
returned, so the output distribution is exactly the one given by the weights,
whatever their size.

A run has four phases:

```
 .dist file ──► reading ──► analysis ──► construction ──► sampling
                (GMP)       backend      alias table      mode speed | alea
                            u64/u128/gmp (Algorithm 1)    (v2)         (v1)
```

* **Analysis** (§3) measures the size of the weights and chooses the smallest
  exact integer representation (*backend*): `uint64_t`, `unsigned __int128`,
  or GMP.
* **Construction** (§4) builds the alias table once, in linear time.
* **Sampling** (§6) draws objects in constant time. Two modes are available:
  * `speed` (v2) minimizes the running time and reads one 64-bit random word
    per attempt;
  * `alea` (v1) minimizes the number of random bits, with Lumbroso's optimal
    uniform and Bernoulli generators.

Both modes sample the exact same distribution and use the same table layout.

Source files:

| File | Role |
|---|---|
| `include/alias.h` | public interface |
| `src/dist.c` | reading and analysis of the `.dist` file, backend selection |
| `src/build.c` | construction for the three backends, slow paths, table information |
| `src/build_native.h` | construction on a machine integer type (included twice, for `uint64_t` and `unsigned __int128`) |
| `src/sample.c` | sampling step, both modes |
| `src/bitstream.{h,c}` | stream of cryptographic random bits |
| `src/alias_internal.h` | table and cell structures |
| `src/main.c` | command-line program `alias_sampler` |
| `tests/test_alias.c` | tests (§8) |
| `tools/bench_backends.c`, `tools/make_bench_dists.py` | comparison of the backends (§9, `make bench`) |
| `tools/make_examples.py` | example distributions of `examples/` |

---

## 2. Input: the `.dist` file

A `.dist` file contains whitespace-separated decimal integers:

```
Z                       declared total weight (informative only)
n w_1 w_2 ... w_n       number of objects, then the n weights (any size, >= 0)
H                       optional: entropy of the distribution
```

This is the format written by the experiment scripts of the paper. The first
line is **not** used for sampling: the weights themselves define the
distribution. If `Z` differs from `Σ w_i`, the analysis reports it.

Reading (`alias_dist_read`) parses the weights with GMP, whatever their size,
and normalizes them in two ways that leave the distribution unchanged:

1. **Division by the gcd** of the weights. Smaller integers may allow a
   smaller backend.
2. **Multiplication by `n`** if `w < n`. The construction needs a cell size
   `cs = ⌊w / n⌋ ≥ 1`. After the multiplication, `w ≥ n`.

Objects of weight 0 are kept (indices stay those of the file) and are never
returned (§4.4).

Limits: `1 ≤ n ≤ 2^32 − 3`.

---

## 3. Analysis and backend selection

`alias_analyze` runs before the construction and reports:

| Field | Meaning |
|---|---|
| `n`, `n_zero` | number of objects, objects of weight 0 |
| `total_bits` | bit length of `w` (after normalization) |
| `max_weight_bits` | bit length of the largest weight |
| `gcd_bits`, `scale` | normalization applied (§2) |
| `cs_bits` | bit length of the cell size `cs = ⌊w / n⌋` |
| `entropy` | Shannon entropy `−Σ p_i log2 p_i`, computed from the weights |
| `declared_total_matches` | first line of the file equals the sum of the original weights |
| `below_double` | number of objects with `n·w_i / w < 2^-53` |
| `backend` | smallest exact backend |

**`below_double`.** In an alias table, object `i` occupies a share
`n·w_i / w` of a cell, and a floating-point sampler decides between the two
objects of a cell by comparing a uniform double with a threshold. A double has
53 significant bits, so an object whose share is below `2^-53` is never drawn
by such a sampler. This count shows how far a double-precision Alias sampler
would be from the target distribution; it is not used by the program itself.

### Backend selection

| Backend | Type | Condition |
|---|---|---|
| `u64` | `uint64_t` | `w < 2^64` |
| `u128` | `unsigned __int128` | `w < 2^128` |
| `gmp` | `mpz_t` | always |

The condition is on the total weight `w` only, because every value handled by
the construction and by the slow paths fits in the type as soon as `w` does:

* the weights, their sum `w`, `cs = ⌊w/n⌋`, `r = w mod cs`, and the weight
  `cs − r` of the virtual object are at most `w`;
* during the construction, the residual weight of an object only decreases:
  `D[x] ← D[x] − (cs − D[x2])` with `D[x] > cs ≥ cs − D[x2]`, and
  `D[x] ← D[x] − cs` with `D[x] > cs`; no intermediate sum is formed;
* the binary long divisions (§4.3) and the slow paths (§6.4) double a
  remainder `t < cs`, and `2t < 2cs ≤ 2w/n ≤ w` for `n ≥ 2` (for `n = 1`, the
  table has a single cell without alias and these steps never run).

The backend can be forced with `--backend` (or the `backend` argument of
`alias_build`), as long as it is not smaller than the minimal one: forcing
`gmp` on small weights is allowed, forcing `u64` on 100-bit weights is refused.

**Effect of the backend.** The backend only changes how the table is built and
stored. The three backends produce **the same table**, cell by cell (this is
tested, §8), and the sampling step never touches the integers except in slow
paths of probability at most `2^-31` per attempt. So the sampling time does not
depend on the backend; the construction time and the memory do (§9).

---

## 4. Construction of the alias table

### 4.1 Cell size and virtual object

The construction (Algorithm 1 of the paper) first computes

```
cs = ⌊w / n⌋          cell size (integer)
q  = ⌊w / cs⌋         number of full cells
r  = w mod cs
```

If `r > 0`, a **virtual object** of weight `cs − r` (index `n`) is added, so
that the total weight becomes `w̄ = w + cs − r = (q + 1)·cs`, and `q` is
incremented. Every cell then has exactly the integer capacity `cs`. The table
has `q` cells, with `n ≤ q ≤ n + ⌈(n − 1)/cs⌉` (tight bounds of the paper).

When a sampling attempt returns the virtual object, it is rejected and the
attempt is repeated (§6.1). Its probability per attempt is
`ρ = (cs − r) / (q·cs) < 1/q`.

### 4.2 Filling the cells

Objects are split in two stacks (LIFO), by their residual weight `D[x]`:

* `L`: objects with `D[x] > cs` (large);
* `S`: objects with `D[x] ≤ cs` (small).

```
while L is not empty:
    x = pop(L)
    if S is not empty:
        x2 = pop(S)
        if D[x2] < cs:
            cell (x2, alias x, threshold T = D[x2])   # x2 w.p. T/cs, x otherwise
            D[x] = D[x] − (cs − D[x2])
        else:                                        # D[x2] == cs
            cell (x2) without alias
    else:
        cell (x) without alias
        D[x] = D[x] − cs
    push x on L if D[x] > cs, on S otherwise
while S is not empty:
    cell (pop(S)) without alias
```

Since the weights are integers and every cell has the integer capacity `cs`,
all operations are exact; the paper proves that the procedure fills exactly
`q` cells (Theorem 1). The implementation checks it.

### 4.3 Precomputed digits

For each cell with an alias, the construction stores the first binary digits
of the threshold `T / cs` (with `0 < T < cs`), obtained by binary long
division: `t ← 2t; if t ≥ cs then (t ← t − cs, digit 1) else digit 0`.

* `alea` mode: the 31 first digits `d_1 … d_31`, plus a flag set when the
  remainder after 31 digits is 0 (finite expansion).
* `speed` mode: the 32 first digits, i.e. `D = ⌊T·2^32 / cs⌋`.

The native backends run the long division on the machine type; the GMP
backend computes `⌊T·2^k / cs⌋` with one GMP division.

### 4.4 Special cases

* An object of weight 0 is small; when it is popped against a large object, it
  gives a cell with threshold 0, which is stored as a cell **without alias**
  holding the large object. A weight-0 object therefore never appears in the
  table and is never returned.
* A small object of weight exactly `cs` gives a cell without alias.
* Cells with an alias are written from the front of the cell array, cells
  without alias from the back. The order of the cells does not matter (a cell
  is drawn uniformly), and this keeps the exact thresholds contiguous (§5.3).
  The resulting table is identical to the one of the paper's implementation
  (`alias.integers_opt`), which was verified on the experimental distributions.

Complexity: `O(n)` operations on the backend integers, plus the digit
computation (32 steps per cell, or one division per cell for GMP).

---

## 5. Data structures

### 5.1 Cell (16 bytes)

```c
typedef struct {
    uint32_t digits;    /* precomputed digits of T / cs (see below)   */
    uint32_t primary;   /* object returned with probability T / cs    */
    uint32_t alias;     /* object returned otherwise; == primary if no alias */
    uint32_t thr;       /* position of the exact threshold T (cells with an alias) */
} alias_cell;
```

| Field | `alea` mode (v1) | `speed` mode (v2) |
|---|---|---|
| `digits` | bits 31…1: `d_1 … d_31` (`d_j` at bit `32 − j`); bit 0: finite-expansion flag | `D = ⌊T·2^32 / cs⌋`; `UINT32_MAX` for a cell without alias |
| `primary`, `alias` | object indices; `alias == primary` marks a cell without alias | same |
| `thr` | index (u64, u128) or byte offset (gmp) of `T` | same |

The cell is 16 bytes, so four cells share a 64-byte cache line, and one cell
access is enough for all but the slow paths.

### 5.2 Table

```c
struct alias_table {
    alias_mode mode;  alias_backend backend;
    uint32_t n;          /* objects */
    uint32_t q;          /* cells */
    uint32_t q_nbits;    /* bit length of q − 1 (Fast Dice Roller) */
    uint32_t q_reject;   /* (2^32 − q) mod q (Lemire's rejection threshold) */
    uint32_t virtual_obj;/* n if there is a virtual object, UINT32_MAX otherwise */
    uint32_t m;          /* cells with an alias: cells 0 .. m−1 */
    alias_cell *cells;   /* q cells */
    double reject_prob;  /* ρ = (cs − r) / (q·cs) */
    /* exact thresholds of cells 0 .. m−1 and cs, in the backend's representation */
    uint64_t cs64;  uint64_t *thr64;          /* u64  */
    u128     cs128; u128     *thr128;         /* u128 */
    mpz_t cs_mpz;   unsigned char *thr_bytes; uint32_t thr_size;   /* gmp */
};
```

### 5.3 Exact thresholds

The exact thresholds are only read by the slow paths. Their storage depends on
the backend:

| Backend | Storage of `T` for cell `i < m` | Size |
|---|---|---|
| `u64` | `thr64[cells[i].thr]` | 8 bytes per aliased cell |
| `u128` | `thr128[cells[i].thr]` | 16 bytes per aliased cell |
| `gmp` | significant bytes of `T`, little-endian, at `thr_bytes + cells[i].thr`; the length is the offset of cell `i + 1` (or `thr_size`) minus that of cell `i` | `⌈bits(T) / 8⌉` bytes per aliased cell |

Memory of a table: `q·16` bytes of cells, plus the thresholds, plus a constant.
The GMP storage only keeps the significant bytes, so it can be smaller than
the fixed-width native storage when many thresholds are small.

---

## 6. Sampling step

### 6.1 Attempts and virtual-object rejection

```
repeat:
    i      = uniform cell in [0, q)
    b      = Bernoulli(T_i / cs)        (1 for a cell without alias)
    result = b ? primary_i : alias_i
until result != virtual object
return result
```

The probability of object `x` in one attempt is `w_x / w̄`; conditioned on
not being the virtual object, it is `w_x / w`. The expected number of attempts
is `1 / (1 − ρ) = w̄ / w ≤ 1 + 1/n` (Lemma 3 of the paper).

### 6.2 Mode `alea` (v1): fewest random bits

**Uniform cell — Fast Dice Roller** (Lumbroso, 2013). The state is a pair
`(x, bound)` with `x` uniform on `[0, bound)`. It starts by reading
`⌈log2 q⌉` bits at once (`bound = 2^⌈log2 q⌉`), then:

```
while bound < q or x >= q:
    if bound >= q: bound −= q; x −= q        # reject, keep the leftover
    bound = 2·bound; x = 2·x + bit()
return x
```

For an isolated draw, its expected number of bits equals the optimum of Knuth
and Yao for the uniform distribution on `q` values (between `log2 q` and
`log2 q + 2`).

**Bernoulli — Lumbroso's algorithm on precomputed digits.** Lumbroso's
Bernoulli(`p`) reads coins `c_1, c_2, …` and returns the digit `d_k` of `p` at
the first coin `c_k = 1`; if the binary expansion of `p` ends at digit `L`, it
returns `c_L` at step `L`. Hence `P(1) = Σ d_k 2^-k = p`, with 2 bits on
average (the Knuth–Yao optimum for a non-dyadic `p`).

The implementation finds the first coin equal to 1 with a single `ctz`
instruction on the current 64-bit word (`bs_coins_until_one`), then reads
`d_k` in `cell.digits`:

* finite expansion ending at digit `L ≤ 31` (flag set): at most `L` coins are
  read; the first 1 at position `k ≤ L` gives `d_k` (for `k = L` this is
  `c_L = 1`), and `L` zero coins give 0;
* otherwise, the first 1 at position `k ≤ 31` gives `d_k`;
* if the 31 coins are 0 (probability `2^-31`), the **slow path**
  `alias_slow_tail` continues Lumbroso's algorithm exactly from the remainder
  `T·2^31 mod cs` (§6.4).

The coins are read in the same order and number as in Lumbroso's algorithm, so
the consumption is exactly that of the paper (Proposition 4: at most
`log2 q + 3` bits per attempt on average).

### 6.3 Mode `speed` (v2): fewest instructions

Each attempt reads **one** 64-bit word `W` (`bs_word`, aligned reads):

**Uniform cell — Lemire's multiply-shift** (Lemire, 2019) on the high half
`h = W >> 32`:

```
m = h·q                       (64-bit product)
i = m >> 32                   (the cell)
if (m mod 2^32) < q:          (probability < q / 2^32)
    while (m mod 2^32) < (2^32 − q) mod q:   m = bits(32)·q
    i = m >> 32
```

This is exactly uniform on `[0, q)`: the rejection removes the
`(2^32 − q) mod q` values of `h` that would bias the result.

**Bernoulli — comparison with 32 digits.** The low half `ℓ = W mod 2^32` is
compared with `D_i = ⌊T_i·2^32 / cs⌋`. Let `V = (ℓ + Θ) / 2^32` with `Θ`
uniform on `[0, 1)`; `V` is uniform on `[0, 1)` and the attempt must return
the primary object iff `V < T_i / cs`:

* `ℓ < D_i` implies `V < (D_i)/2^32 ≤ T_i/cs`: primary object;
* `ℓ > D_i` implies `V ≥ (D_i + 1)/2^32 > T_i/cs`: alias;
* `ℓ = D_i` (probability `2^-32`): the **slow path** `alias_slow_tie` returns
  Bernoulli(`(T_i·2^32 mod cs) / cs`) with fresh bits, the exact probability
  that `V < T_i / cs` given `ℓ = D_i`.

The low half is independent of the high half, including after a Lemire
rejection (which only redraws the high part). A cell without alias has
`D = UINT32_MAX` and `alias == primary`, so the selection is branch-free.

Consumption: 64 bits per attempt (plus rare redraws), about 3.6 times more
than `alea` on the experimental distributions, for about half the time (§9).

### 6.4 Slow paths

Both slow paths reduce to Lumbroso's Bernoulli continued from a remainder
`rem = T·2^k mod cs` (`k = 31` for `alea`, `k = 32` for `speed`):

```
if rem == 0: return 0
loop:
    rem = 2·rem
    if rem == cs: return bit()          # remaining expansion is 0.1000…
    y = (rem > cs); if y: rem −= cs     # next digit
    if bit() == 1: return y
```

The remainder is computed on demand from the exact threshold, in the
backend's arithmetic (`remainder_u64`, `remainder_u128`, or GMP). These paths
run with probability at most `2^-31` per attempt; they are tested directly
(§8).

---

## 7. Random-bit stream

The random bits come from a cryptographically secure generator
(`arc4random_buf` on macOS and BSD, `getrandom` on Linux), through a ring
buffer filled by a background thread.

```c
typedef struct {
    uint64_t tail;        /* producer: words written      (own cache line) */
    uint64_t head;        /* consumer: words released     (own cache line) */
    uint64_t c_head;      /* consumer: current word index                  */
    uint64_t c_tail;      /* consumer: last tail seen                      */
    uint32_t c_off;       /* consumer: bits read in the current word       */
    uint32_t c_publish;   /* consumer: words since the last release        */
    volatile int running;
    uint64_t data[32768]; /* 256 KiB ring buffer                           */
} bs_queue_t;
```

* **Producer.** While at least 512 words (4 KiB) are free, it fills them with
  one call to the generator and publishes `tail` (release). When the buffer is
  full (about two million bits ahead), it sleeps 50 µs instead of spinning.
* **Consumer.** It reads bits from the current word with inline functions
  (`bs_bits`, `bs_bit`, `bs_word`, `bs_coins_until_one`), least significant bit
  first. It reads `tail` (acquire) only when it has consumed every word it
  knows of, and releases consumed words every 64 words.
* **Accounting.** `bs_bits_consumed() = 64·c_head + c_off` is the exact
  number of bits consumed; `alias_rng_bits()` exposes it.

Generating the bits costs about 0.8 ns per 32 bits on an Apple M4 Pro; this
time is spent on the producer thread, in parallel with sampling.

There is one stream per process and one consumer thread: `alias_sample` is not
thread-safe.

---

## 8. Exactness and tests

`make test` runs `tests/test_alias` on the distributions of `examples/`:

1. **Backends agree.** For each distribution and mode, every backend at least
   as large as the minimal one builds the same table as GMP: same number of
   cells, same virtual object, same `cs`, and for each cell the same primary
   object, alias, digits and exact threshold.
2. **Goodness of fit.** For each mode and backend, 2·10^6 samples (5·10^6 for
   the large distribution) are compared with the exact probabilities by a
   chi-square test (objects with an expected count below 5 pooled into one
   class; Wilson–Hilferty approximation, failure for `z > 4.5`). The test also
   checks that no index `≥ n` (virtual object) and no object of weight 0 is
   ever returned.
3. **Slow paths.** They are too rare to be reached by sampling, so
   `alias_slow_tail` and `alias_slow_tie` are called directly, 10^5 times on
   five cells, and their frequency of 1 is compared with
   `(T·2^k mod cs) / cs` (failure for `|z| > 5`).

The examples cover the three backends, weights of 0, a total smaller than `n`
(scaling), and a distribution of the paper's protocol (`n = 20 484`,
`w ≈ 4.7·10^50`). The tables were also checked to be identical, cell by cell,
to those of the paper's implementations (`alias.integers_opt` for `alea`,
`alias.integers_v2` for `speed`) on distributions of the experimental protocol.

---

## 9. Performance

All measurements: Apple M4 Pro, one process, best of 3 runs, 10^7 samples.
Sampling times vary by up to about ±30 % between sessions; bit counts are
exact. `make bench` reproduces the comparison of the backends
(`tools/make_bench_dists.py` writes the distributions, `tools/bench_backends`
measures them).

### 9.1 Comparison of the backends

Twelve random distributions: `n ∈ {10^3, 10^5, 10^6}`, weights of 30 bits
(total < 2^64: every backend applies) or 90 bits (u128 and gmp only), either
uniform on `[1, 2^bits)` (*uniform*) or uniform on `[1, 2^b)` with a random
`b < bits` for each object (*skewed*).

**Construction time** (Algorithm 1 and the digits):

| n | u64 | u128 | gmp |
|---|---|---|---|
| 10^3 | 0.04 ms | 0.07 ms | 0.10–0.14 ms |
| 10^5 | 3.7–3.8 ms | 6.9–7.3 ms | 9.3–13.4 ms |
| 10^6 | 37–38 ms | 70–72 ms | 92–131 ms |

`u64` builds about 2.5 times faster than `gmp`, `u128` 1.3 to 1.8 times
faster than `gmp`; `u128` costs about twice `u64`.

**Table memory** (n = 10^6):

| Weights | u64 | u128 | gmp |
|---|---|---|---|
| 30 bits, uniform | 23.4 MB | 31.3 MB | 19.5 MB |
| 30 bits, skewed | 23.4 MB | 31.3 MB | 17.7 MB |
| 90 bits, uniform | — | 31.3 MB | 27.0 MB |
| 90 bits, skewed | — | 31.3 MB | 21.4 MB |

The cells take 16 bytes in every case; the difference comes from the exact
thresholds (§5.3): 8 or 16 bytes each for `u64` and `u128`, only their
significant bytes for `gmp`, which is therefore 15 to 45 % smaller.

**Sampling.** The sampling time is the same for the three backends, within
±2 %, as expected: the sampling code is common and only the slow paths
(probability ≤ 2^-31) touch the integers.

| n | `alea`, uniform weights | `alea`, skewed weights | `speed` | bits: `alea` / `speed` |
|---|---|---|---|---|
| 10^3 | ≈ 16 ns | ≈ 10–11 ns | ≈ 6.6 ns | 12.2 / 64 |
| 10^5 | ≈ 22 ns | ≈ 14–16 ns | ≈ 8.3 ns | 19.7 / 64 |
| 10^6 | ≈ 22 ns | ≈ 13–15 ns | ≈ 11 ns | 22.3 / 64 |

* `alea` is about 1.5 times slower on uniform weights: the thresholds are then
  spread over `(0, 1)` rather than close to 0 or 1, and the outcome of the
  Bernoulli becomes unpredictable for the branch predictor. `speed` does not
  depend on the shape of the distribution.
* Both modes slow down as `n` grows, because of the cache: at `n = 10^6`, the
  cells alone take 16 MB.
* `alea` consumes the Fast Dice Roller's cost (between `log2 q` and
  `log2 q + 2` bits) plus about 2 bits for the Bernoulli; `speed` consumes 64.

**Choice of the backend.** The automatic choice (smallest exact
representation) gives the fastest construction. `gmp` would only be preferable
when the memory of the table matters more than its construction time.

### 9.2 Distribution of the paper's protocol

`examples/dirichlet_20484.dist` (`n = 20 484`, 169-bit total, backend `gmp`):

| Mode | ns / sample | bits / sample |
|---|---|---|
| `alea` | ≈ 10–15 | 17.60 |
| `speed` | ≈ 4.5–7.6 | 64.00 |

Construction: 1.7 ms; table memory: 0.42 MB.

---

## 10. Limits

* `n ≤ 2^32 − 3` and `q < 2^32 − 1` (indices on 32 bits).
* One consumer thread per process (§7); for parallel sampling, the stream
  would have to be made per-thread.
* The random bits are cryptographic and cannot be seeded: runs are not
  reproducible bit for bit. Exactness does not depend on it.
* The native backends store the exact thresholds at fixed width (8 or 16
  bytes), the GMP backend at variable width: `gmp` uses 15 to 45 % less memory
  than `u64` and `u128` on the distributions of §9.1. Storing the native
  thresholds by significant bytes too would remove this difference without
  affecting the sampling time.

---

## 11. Interfaces

### Command line

```
alias_sampler [options] FILE.dist
  -m, --mode speed|alea         default speed
  -b, --backend auto|u64|u128|gmp   default auto (from the analysis)
  -n, --samples N               default 10
  -o, --output samples|histogram|none
  -a, --analyze-only            print the analysis and stop
  -q, --quiet                   nothing on stderr
```

The analysis, the table information and the sampling statistics (time per
sample, random bits per sample) are written on stderr; the samples (one index
per line, from 0) or the histogram (`index count` lines) on stdout.

### Library

```c
#include "alias.h"

alias_dist d;
char err[256];
if (alias_dist_read("weights.dist", &d, err, sizeof err) != 0) { /* err */ }

alias_report r;
alias_analyze(&d, &r);                               /* r.backend: u64, u128 or gmp */

alias_table *t = alias_build(&d, ALIAS_ALEA, r.backend);
alias_rng_start();
for (int k = 0; k < 1000; k++) {
    uint32_t x = alias_sample(t);                    /* 0 <= x < d.n */
}
alias_rng_stop();
alias_table_free(t);
alias_dist_free(&d);
```

Link with `libalias.a -lgmp -lpthread -lm`.

### References

* D. Lemire. Fast Random Integer Generation in an Interval. *ACM TOMACS* 29(1), 2019.
* J. Lumbroso. Optimal Discrete Uniform Generation from Coin Flips, and Applications. arXiv:1304.1916, 2013.
* M. D. Vose. A Linear Algorithm for Generating Random Numbers with a Given Distribution. *IEEE TSE* 17(9), 1991.
* A. J. Walker. An Efficient Method for Generating Discrete Random Variables with General Distributions. *ACM TOMS* 3(3), 1977.
* D. E. Knuth, A. C. Yao. The Complexity of Nonuniform Random Number Generation. 1976.
