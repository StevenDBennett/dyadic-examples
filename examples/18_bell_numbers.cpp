#include "dyadic.h"
#include <cstdint>
#include <cstdio>

using namespace dyadic;

// Bell number B_n = number of partitions of a set with n elements.
// Dobinski's formula: B_n = (1/e) * Σ_{k=0}^∞ k^n / k!
// Also: B_n = Σ_{k=0}^n S(n,k) where S(n,k) are Stirling numbers of the 2nd kind.

template<std::unsigned_integral W, int MaxN = 15>
struct BellCache {
    std::array<W, MaxN> B{};

    constexpr BellCache() noexcept {
        // B_0 = 1
        B[0] = 1;

        // Recurrence using Stirling numbers of the 2nd kind
        for (int n = 1; n < MaxN; ++n) {
            W sum = 0;
            for (int k = 1; k <= n; ++k) {
                sum += stirling_2<W, MaxN>(n, k);
            }
            B[n] = sum;
        }
    }
};

int main() {
    std::printf("=== 18 — Bell Numbers ===\n\n");

    using W = uint64_t;
    constexpr int MaxN = 12;

    constexpr auto cache = BellCache<W, MaxN>{};

    std::printf("Bell numbers (number of set partitions):\n");
    for (int n = 0; n < MaxN; ++n) {
        std::printf("  B(%d) = %lu\n", n, (unsigned long)cache.B[n]);
    }

    std::printf("\nKey formulas:\n");
    std::printf("  B_n = Σ_{k=0}^n S(n,k)  (sum of Stirling numbers of 2nd kind)\n");
    std::printf("  B_{n+1} = Σ_{k=0}^n C(n,k) B_k  (explicit recurrence)\n\n");

    // Verify the recurrence
    std::printf("Verification of B_{n+1} = Σ C(n,k) B_k:\n");
    for (int n = 0; n < MaxN - 1; ++n) {
        W sum = 0;
        for (int k = 0; k <= n; ++k) {
            W c = binom<W, MaxN>(n, k);
            sum += c * cache.B[k];
        }
        std::printf("  n=%d: Σ C(n,k) B_k = %lu, B_{%d} = %lu %s\n",
               n, (unsigned long)sum, n+1, (unsigned long)cache.B[n+1],
               (sum == cache.B[n+1]) ? "✓" : "×");
    }

    // Show Stirling contributions
    std::printf("\nStirling number decomposition S(n,k):\n");
    for (int n = 1; n < 8; ++n) {
        std::printf("  B(%d) = ", n);
        W total = 0;
        for (int k = 1; k <= n; ++k) {
            W s = stirling_2<W, MaxN>(n, k);
            total += s;
            std::printf("S(%d,%d)=%lu", n, k, (unsigned long)s);
            if (k < n) std::printf(" + ");
        }
        std::printf(" = %lu\n", (unsigned long)total);
    }

    return 0;
}
