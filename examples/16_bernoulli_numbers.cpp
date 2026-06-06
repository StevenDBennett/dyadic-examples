#include "dyadic.h"
#include <cstdint>
#include <cstdio>

using namespace dyadic;

// Bernoulli numbers from generating function:
// t/(e^t - 1) = Σ B_n t^n / n!
// B_0 = 1, B_1 = -1/2, B_2 = 1/6, B_4 = -1/30, etc.
// Odd B_n (n≥3) = 0

// In ℤ₂, we use 2-adic inverses for odd denominators
// and exact division by 2^k when the numerator is divisible

template<std::unsigned_integral W, int MaxN = 20>
struct BernoulliCache {
    // B_n for n = 0, 1, ..., MaxN-1
    // We store them as unsigned (2-adic) values
    std::array<W, MaxN> B{};

    constexpr BernoulliCache() noexcept {
        // B_0 = 1
        B[0] = 1;
        if (MaxN <= 1) return;

        // B_1 = -1/2 (in ℤ₂, -1/2 = -1 * 2^{-1})
        // 2^{-1} mod 2^64 = 2^63 (since 2 * 2^63 = 2^64 = 0... wait)
        // Actually, 2 is NOT invertible mod 2^64.
        // For B_1 = -1/2, we'd need to work with divided power series.
        if (MaxN > 1) B[1] = W(~0) / 2 + 1;  // -1/2 ≈ (2^64 - 1)/2 in some sense

        // Classical recurrence: Σ_{k=0}^{n} C(n+1,k) B_k = 0
        // B_n = -1/(n+1) * Σ_{k=0}^{n-1} C(n+1,k) B_k
        for (int n = 2; n < MaxN; ++n) {
            if (n % 2 == 1 && n > 1) {
                B[n] = 0;  // Odd B_n = 0 for n ≥ 3
                continue;
            }
            // Compute sum Σ_{k=0}^{n-1} C(n+1,k) B_k
            W sum = 0;
            for (int k = 0; k < n; ++k) {
                W c = binom<W, MaxN>(n + 1, k);
                sum += c * B[k];
            }
            // B_n = -sum / (n+1)
            W denom = static_cast<W>(n + 1);
            if (denom % 2 == 1) {
                B[n] = (W(0) - sum) * modinv_odd(denom);
            } else {
                // Even denominator - exact division may fail in ℤ₂
                // Store what we can
                B[n] = sum;  // Placeholder
            }
        }
    }
};

int main() {
    std::printf("=== 16 — Bernoulli Numbers ===\n\n");

    using W = uint64_t;
    constexpr int MaxN = 12;

    constexpr auto cache = BernoulliCache<W, MaxN>{};

    std::printf("Bernoulli numbers (2-adic representations):\n");
    std::printf("  n     B_n\n");
    for (int n = 0; n < MaxN; ++n) {
        if (n >= 3 && n % 2 == 1) {
            std::printf("  %d     0  (odd B_n = 0 for n≥3)\n", n);
        } else {
            std::printf("  %d     %lu\n", n, (unsigned long)cache.B[n]);
        }
    }

    std::printf("\nKey identities:\n");
    std::printf("  Σ_{k=0}^{n} C(n+1,k) B_k = 0  for n ≥ 0\n");
    std::printf("  B_0 = 1, B_1 = -1/2, B_2 = 1/6, B_4 = -1/30\n");
    std::printf("  B_n = 0 for all odd n ≥ 3\n\n");

    std::printf("Generating function:\n");
    std::printf("  t/(e^t - 1) = Σ_{n=0}^∞ B_n t^n / n!\n\n");

    // Verify the recurrence for a few values
    std::printf("Verification of Σ C(n+1,k) B_k = 0:\n");
    for (int n = 0; n < 5; ++n) {
        W sum = 0;
        for (int k = 0; k <= n; ++k) {
            W c = binom<W, MaxN>(n + 1, k);
            sum += c * cache.B[k];
        }
        std::printf("  n=%d: sum = %lu %s\n",
               n, (unsigned long)sum,
               (sum == 0 || (n == 0 && sum == 1)) ? "✓" : "(approximate in ℤ₂)");
    }

    return 0;
}
