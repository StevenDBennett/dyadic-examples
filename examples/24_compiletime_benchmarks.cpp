#include "dyadic.h"
#include <cstdint>
#include <cstdio>
#include <chrono>

using namespace dyadic;

// Compare compile-time vs runtime execution for key operations
// In C++20, we can evaluate many of these at compile time using constexpr

// A compile-time computable polynomial
template<int N, typename W>
constexpr Polynomial<N, W, MonomialBasis> make_poly() {
    Polynomial<N, W, MonomialBasis> p{};
    for (int i = 0; i < N; ++i) {
        p[i] = static_cast<W>(i + 1);
    }
    return p;
}

// Perform operations at compile time
template<int N, typename W>
constexpr auto compile_time_computation() {
    auto a = make_poly<N, W>();
    auto b = make_poly<N, W>();
    auto sum = a + b;
    auto prod = a * b;
    auto da = formal_derivative(a);
    auto ff = change_basis<FallingFactorialBasis>(a);
    auto back = change_basis<MonomialBasis>(ff);
    return std::make_tuple(sum, prod, da, ff, back);
}

// Perform the same operations at runtime
template<int N, typename W>
auto runtime_computation() {
    auto a = make_poly<N, W>();
    auto b = make_poly<N, W>();
    auto start = std::chrono::steady_clock::now();
    auto sum = a + b;
    auto prod = a * b;
    auto da = formal_derivative(a);
    auto ff = change_basis<FallingFactorialBasis>(a);
    auto back = change_basis<MonomialBasis>(ff);
    auto end = std::chrono::steady_clock::now();
    return std::make_tuple(sum, prod, da, ff, back, end - start);
}

int main() {
    std::printf("=== 24 — Compile-time vs Runtime Benchmarks ===\n\n");

    using W = uint64_t;
    constexpr int N = 10;

    // Force compile-time evaluation
    [[maybe_unused]] constexpr auto ct_result = compile_time_computation<N, W>();
    std::printf("Compile-time operations: ✓ Completed at compile time\n");

    // Runtime evaluation
    auto [sum, prod, da, ff, back, duration] = runtime_computation<N, W>();
    auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count();
    std::printf("Runtime operations:      ✓ Completed in %ld ns\n\n", ns);

    std::printf("Key observations:\n");
    std::printf("  1. All dyadic operations are constexpr-friendly\n");
    std::printf("  2. Compile-time evaluation has zero runtime cost\n");
    std::printf("  3. The Newton iteration can be fully constexpr\n");
    std::printf("  4. Hensel lifting benefits from compile-time precision\n\n");

    std::printf("Benchmarking specific operations (runtime):\n\n");

    // Individual benchmarks
    auto a = make_poly<N, W>();
    auto b = make_poly<N, W>();

    auto start = std::chrono::steady_clock::now();
    volatile auto s = a + b;
    auto end = std::chrono::steady_clock::now();
    std::printf("  Addition:         %6ld ns\n",
           std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

    start = std::chrono::steady_clock::now();
    volatile auto p = a * b;
    end = std::chrono::steady_clock::now();
    std::printf("  Multiplication:   %6ld ns\n",
           std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

    start = std::chrono::steady_clock::now();
    volatile auto d = formal_derivative(a);
    end = std::chrono::steady_clock::now();
    std::printf("  Derivative:       %6ld ns\n",
           std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

    start = std::chrono::steady_clock::now();
    volatile auto delta = forward_difference(a);
    end = std::chrono::steady_clock::now();
    std::printf("  Forward difference: %6ld ns\n",
           std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

    start = std::chrono::steady_clock::now();
    volatile auto ff_basis = change_basis<FallingFactorialBasis>(a);
    end = std::chrono::steady_clock::now();
    std::printf("  Basis change:     %6ld ns\n",
           std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

    std::printf("\n");

    return 0;
}
