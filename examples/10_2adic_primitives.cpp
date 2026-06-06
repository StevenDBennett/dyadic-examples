#include "dyadic.h"
#include <cstdio>

int main() {
    using namespace dyadic;

    std::printf("=== 10 — 2-Adic Primitives ===\n\n");

    std::printf("v2(x) — 2-adic valuation (number of trailing zeros):\n");
    for (uint64_t x : {uint64_t(0), uint64_t(1), uint64_t(2), uint64_t(12),
                       uint64_t(48), uint64_t(1024), uint64_t(1) << 63}) {
        std::printf("  v2(%20lu) = %d\n", x, v2(x));
    }
    std::printf("  (v2(0) = word width by convention: %d)\n\n", v2(uint64_t(0)));

    std::printf("modinv_odd(a) — modular inverse of odd a mod 2^W:\n");
    for (uint64_t a : {uint64_t(1), uint64_t(3), uint64_t(5), uint64_t(7), uint64_t(9),
                       uint64_t(0xFFFFFFFFFFFFFFFFULL)}) {
        if (a % 2 == 0) continue;
        auto inv = modinv_odd(a);
        auto prod = a * inv;
        std::printf("  %lu⁻¹ mod 2⁶⁴ = %lu  (verify: %lu·%lu = %lu)\n", a, inv, a, inv, prod);
    }
    std::printf("\n");

    std::printf("exact_divide(x, d) — exact 2-adic division x/d:\n");
    struct { uint64_t x, d; } div_cases[] = {{42, 3}, {100, 5}, {0xFFFF, 3}, {1ULL << 63, 3}};
    for (auto [x, d] : div_cases) {
        auto q = exact_divide(x, d);
        auto check = q * d;
        std::printf("  %lu / %lu = %lu  (verify: %lu·%lu = %lu) %s\n",
            x, d, q, q, d, check, check == x ? "✓" : "✗");
    }
    std::printf("\n");

    std::printf("Artin-Schreier operator ℘(x) = x² − x:\n");
    std::printf("  Kernel: ℘(0) = %lu, ℘(1) = %lu\n", artin_schreier(uint64_t(0)), artin_schreier(uint64_t(1)));
    std::printf("  Symmetry: ℘(x) = ℘(1−x):\n");
    for (uint64_t x : {uint64_t(2), uint64_t(5), uint64_t(10), uint64_t(100), uint64_t(0xDEADBEEF)}) {
        auto p1 = artin_schreier(x);
        auto p2 = artin_schreier(uint64_t(1) - x);
        std::printf("    ℘(0x%lx) = 0x%lx, ℘(1−x) = 0x%lx  %s\n", x, p1, p2, p1 == p2 ? "✓" : "✗");
    }

    return 0;
}
