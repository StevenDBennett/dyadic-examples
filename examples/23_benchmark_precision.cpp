#include "dyadic.h"
#include <cstdint>
#include <cstdio>
#include <chrono>

using namespace dyadic;

// Benchmark the performance of key operations across different word sizes

template<std::unsigned_integral W>
void benchmark_operations() {
    constexpr int N = 10;
    using namespace std::chrono;

    Polynomial<N, W, MonomialBasis> a, b;
    for (int i = 0; i < N; ++i) {
        a[i] = static_cast<W>(i + 1);
        b[i] = static_cast<W>(N - i);
    }

    auto start = steady_clock::now();

    // Addition
    [[maybe_unused]] auto sum = a + b;

    // Multiplication
    [[maybe_unused]] auto prod = a * b;

    // Formal derivative
    [[maybe_unused]] auto da = formal_derivative(a);

    // Forward difference
    [[maybe_unused]] auto delta_a = forward_difference(a);

    // Change of basis
    [[maybe_unused]] auto ff = change_basis<FallingFactorialBasis>(a);
    [[maybe_unused]] auto back = change_basis<MonomialBasis>(ff);

    auto end = steady_clock::now();
    auto duration = duration_cast<nanoseconds>(end - start);

    std::printf("  Word size: %lu bits, Time: %lu ns\n",
           (unsigned long)(8 * sizeof(W)), (unsigned long)duration.count());
}

template<typename T>
void run_benchmark() {
    benchmark_operations<T>();
}

int main() {
    std::printf("=== 23 — Precision and Performance Benchmarks ===\n\n");

    std::printf("Benchmarking operations across word sizes:\n\n");

    std::printf("uint8_t (8 bits):\n");
    run_benchmark<uint8_t>();

    std::printf("uint16_t (16 bits):\n");
    run_benchmark<uint16_t>();

    std::printf("uint32_t (32 bits):\n");
    run_benchmark<uint32_t>();

    std::printf("uint64_t (64 bits):\n");
    run_benchmark<uint64_t>();

    std::printf("\nPrecision analysis:\n");
    std::printf("  The precision window for operations depends on the word size.\n");
    std::printf("  - Taylor basis: T_k = k! * FF_k must avoid overflow\n");
    std::printf("  - Witt vectors: ghost recovery needs r_j < 2^(word_bits - j)\n");
    std::printf("  - Higher word sizes = larger precision windows\n\n");

    std::printf("  Word Size | Max Exact Taylor Order | Max Witt Length\n");
    std::printf("  ----------|----------------------|----------------\n");
    std::printf("  8 bits    | ~5                   | ~8\n");
    std::printf("  16 bits   | ~8                   | ~16\n");
    std::printf("  32 bits   | ~12                  | ~32\n");
    std::printf("  64 bits   | ~20                  | ~64\n\n");

    std::printf("Note: Exact limits depend on coefficient magnitudes.\n");
    std::printf("      Use dyadic::check_taylor_roundtrip_precision() and\n");
    std::printf("      dyadic::check_witt_recovery_precision() to verify.\n\n");

    return 0;
}
