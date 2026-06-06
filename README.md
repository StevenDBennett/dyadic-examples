# dyadic-examples — Examples and Extensions for the dyadic Library

[dyadic](https://github.com/StevenDBennett/dyadic) is a single-header C++20 library for the Six-Axiom 2-Adic Operator Calculus. This project provides 30 runnable examples and optional extension headers.

## Extensions

| Header | Description |
|--------|-------------|
| `dynamic_polynomial.h` | Runtime-degree polynomial with basis conversion, D, Δ |
| `pade.h` | Padé approximant construction in ℤ₂[[t]] |
| `continued_fractions.h` | Continued fraction expansion and convergent recurrence |
| `matrix.h` | Matrix/linear algebra over ℤ/2^Wℤ |

## Build

```bash
cmake -B build
cmake --build build

# Run all 30 examples
cmake --build build --target run
```

## Requirements

- C++20 compiler (GCC 12+, Clang 17+)
- [dyadic](https://github.com/StevenDBennett/dyadic) (cloned alongside as `../dyadic/`)

## License

MIT — see LICENSE.
