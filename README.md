# dyadic-examples

Examples and extension headers for [**dyadic**](https://github.com/StevenDBennett/dyadic) — the Six-Axiom 2-Adic Operator Calculus library.

```bash
cmake -B build && cmake --build build && cmake --build build --target run
```

---

## What This Is

[dyadic](https://github.com/StevenDBennett/dyadic) is a single-header C++20 library for arithmetic and calculus over the 2-adic integers ℤ₂ and the ring of formal power series ℤ₂[[t]]. This companion project provides **31 runnable examples** covering the full dyadic API, plus **extension headers** that add functionality best kept separate from the core library.

## Extensions

Four standalone headers, drop-in ready — just `#include <dynamic_polynomial.h>` (or whichever) in your project:

| Header | What It Adds | Depends On |
|--------|-------------|------------|
| [`dynamic_polynomial.h`](#dynamic_polynomial) | Heap-allocated polynomial with runtime-variable degree | dyadic.h |
| [`pade.h`](#pade) | Padé approximant construction in ℤ₂[[t]] | dyadic.h |
| [`continued_fractions.h`](#continued_fractions) | Continued fraction expansion and convergent evaluation in ℤ₂[[t]] | dyadic.h, dynamic_polynomial.h |
| [`matrix.h`](#matrix) | Generic M×N matrix over ℤ/2^Wℤ with linear algebra | dyadic.h |

---

## `dynamic_polynomial.h`

Runtime-degree polynomial — the same `eval`, formal derivative, forward difference, and basis-conversion operations as the compile-time `Polynomial<N,W,Basis>`, but with heap-allocated storage so the degree is determined at runtime.

```cpp
DynamicPolynomial<uint32_t> p({1, 2, 3});          // 1 + 2x + 3x²
auto df = formal_derivative(p);                      // 2 + 6x
auto taylor = change_basis<TaylorBasis>(p);

// Convert to/from compile-time Polynomial
Polynomial<4, uint32_t> static_p{{1, 2, 3, 4}};
auto dyn = to_dynamic(static_p);
auto back = to_static<4>(dyn);
```

Arithmetic (`+`, `−`, `*`) works between `DynamicPolynomial` values of any degree; mixed-degree operations auto-resize the result.

---

## `pade.h`

Constructs the [m/n] Padé approximant `P(t)/Q(t)` from the first `m+n+1` terms of a formal power series `A(t) = Σ a_i t^i`, matching `A(t)` to order `O(t^{m+n})`. Uses Gaussian elimination over ℤ₂ to solve for Q's coefficients (requires `a_m` odd for a unique normal solution with `Q(0)=1`).

```cpp
// Geometric series: 1 + t + t²
Polynomial<3, uint32_t> geom{{1, 1, 1}};
auto [P, Q] = pade_approximant<1, 1>(geom);
// P = 1, Q = 1 − t   →  1/(1-t) recovered
```

Works at compile time — all operations are `constexpr`.

---

## `continued_fractions.h`

Expands a power series `A(t)` into the continued fraction form

```
A(t) = c₀ + t / (c₁ + t / (c₂ + t / (...)))
```

by repeated constant-term extraction and degree reduction. The k-th convergent `P_k/Q_k` is computed via the recurrence:

```
P₋₁ = 1,  P₀ = c₀,  Q₋₁ = 0,  Q₀ = 1
P_k = c_k · P_{k-1} + t · P_{k-2}
Q_k = c_k · Q_{k-1} + t · Q_{k-2}
```

```cpp
DynamicPolynomial<uint32_t> geom({1, 1, 1, 1, 1, 1, 1, 1});
auto c = cf_expand(geom, 6);                     // {1, 1, 1, 1, 1, 1}
auto [P, Q] = cf_convergent(c, 3);               // 3rd convergent
W val = P.eval(2) * modinv_odd(Q.eval(2));       // evaluate at t=2
```

Rational functions produce finitely terminating CF expansions; transcendental series give infinite ones.

---

## `matrix.h`

Generic `Matrix<M, N, W>` over ℤ/2^Wℤ with Gaussian elimination using odd-pivot selection (invertible elements in ℤ₂ are exactly the odd integers).

```cpp
Matrix<3, 3, uint32_t> A{{{ {1, 2, 3}, {0, 1, 4}, {5, 6, 0} }}};
auto det = A.determinant();
auto inv = A.inverse();    // M×M, Gauss-Jordan
auto x   = A.solve({7, 11, 13});  // linear system
int r    = A.rank();       // Gaussian elimination mod 2
```

All operations (`+`, `−`, `*`, determinant, rank, inverse, solve, transpose, trace, rref) are `constexpr`.

---

## The 31 Examples

Each example is a self-contained `.cpp` file in `examples/`. Build and run them all with `cmake --build build --target run`.

| # | Example | What It Demonstrates |
|---|---------|---------------------|
| 01 | `basis_conversion` | Monomial ↔ FallingFactorial ↔ Taylor roundtrip |
| 02 | `formal_derivative` | D = d/dt, Dⁿ annihilates deg N-1 polynomials |
| 03 | `forward_difference` | Δ = e^D − I, identity ΔP = P(t+1) − P(t) |
| 04 | `taylor_shift` | P(t) → P(t + δ) in monomial and FF bases |
| 05 | `indefinite_sum` | Σ = Δ⁻¹, exact inverse in FF basis |
| 06 | `witt_vectors` | Witt addition/multiplication, ghost map, Frobenius, Verschiebung |
| 07 | `adams_teichmuller` | Adams operations ψⁿ, Teichmüller lifts |
| 08 | `compose_reversion` | Power series composition P(Q(t)), Lagrange inversion |
| 09 | `carry_chain` | C = (I−N)⁻¹, one-pass carry propagation |
| 10 | `2adic_primitives` | v₂, modinv_odd, exact_divide, Artin-Schreier ℘(x) |
| 11 | `stirling_numbers` | S₂(n,k), s₁(n,k), |s₁|(n,k) tables |
| 12 | `polynomial_arithmetic` | +, −, ×, eval in all three bases |
| 13 | `log_exp_series` | Binomial series (1+t)^α, geometric series |
| 14 | `newton_iteration` | Newton's method for sqrt, cuberoot in ℤ₂ |
| 15 | `pade_approximants` | Rational evaluation, geometric series in ℤ₂[[t]] |
| 16 | `bernoulli_numbers` | Bₙ via generating function, Σ C(n+1,k)B_k = 0 |
| 17 | `eulerian_numbers` | A(n,k) table, Worpitzky identity |
| 18 | `bell_numbers` | Bₙ via Stirling sum, recurrence B_{n+1} = Σ C(n,k)B_k |
| 19 | `hensel_lifting` | Root lifting mod 2^k for polynomials over ℤ₂ |
| 20 | `witt_division` | Witt vector inverse a⁻¹, ghost inverse verification |
| 21 | `resultant_discriminant` | Resultant, discriminant via Sylvester matrix |
| 22 | `continued_fractions` | CF expansion, convergent P/Q, finite CF for rationals |
| 23 | `benchmark_precision` | Timing across word sizes (uint8–uint64), precision analysis |
| 24 | `compiletime_benchmarks` | Compile-time vs runtime operation timing |
| 25 | `csv_output` | Generates CSV files for Stirling/Witt/precision visualization |
| 26 | `witt_log_exp` | Witt inverse, exp/log homomorphisms, Adams compose |
| 27 | `dynamic_polynomial` | Runtime-degree polynomial, basis conversion, D, Δ |
| 28 | `hardware_arithmetic` | adc, add_overflow, mul_overflow, carry chain integration |
| 29 | `pade_approximant` | Padé [1/1], [2/2], rational recovery, constexpr usage |
| 30 | `matrix` | Matrix operations, determinant, inverse, rank, solve |
| 31 | `discrete_hedging` | Δ/Σ operators in finance: daily returns, discrete gamma, cumulative P&L, continuous vs discrete hedge ratios |

## Build

### Dependencies

- **C++20** compiler (GCC 12+, Clang 17+)
- [dyadic](https://github.com/StevenDBennett/dyadic) — expected at `../dyadic/dyadic.h` relative to this repo

### Quick start

```bash
cmake -B build
cmake --build build

# Run all 31 examples
cmake --build build --target run

# Run a single example
./build/01_basis_conversion
```

The `run` target runs every example in sequence and prints PASS/FAIL results.

### Project integration

If you already have a CMake project, link against the `dyadic` interface target:

```cmake
add_subdirectory(path/to/dyadic-examples)
target_link_libraries(my_app PRIVATE dyadic)
```

This gives you the `-I` paths for both `<dyadic.h>` and the extension headers in a single step.

## Project Structure

```
dyadic-examples/
├── CMakeLists.txt           # Build system for all 31 examples
├── cmake/
│   └── run_all.sh           # Script invoked by `--target run`
├── examples/
│   ├── 01_basis_conversion.cpp
│   ├── 02_formal_derivative.cpp
│   ├── ...                  # (31 .cpp files total)
│   ├── 30_matrix.cpp
│   └── 31_discrete_hedging.cpp
├── dynamic_polynomial.h     # Extension: runtime-degree polynomials
├── pade.h                   # Extension: Padé approximants
├── continued_fractions.h    # Extension: continued fractions in ℤ₂[[t]]
├── matrix.h                 # Extension: matrix / linear algebra
├── LICENSE                  # MIT
└── README.md
```

## Related

- [dyadic](https://github.com/StevenDBennett/dyadic) — the core library
- [dyadic_verify.h](https://github.com/StevenDBennett/dyadic/blob/main/dyadic_verify.h) — compile-time proofs

## License

MIT — see [LICENSE](LICENSE).
