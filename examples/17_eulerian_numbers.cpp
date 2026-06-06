#include "dyadic.h"
#include <cstdint>
#include <cstdio>

using namespace dyadic;

// Eulerian numbers <n, k> count permutations of {1,...,n} with exactly k descents.
// Recurrence: A(n,k) = (k+1) A(n-1,k) + (n-k) A(n-1,k-1)
// With A(0,0) = 1, A(n,0) = 1 for all n, and A(n,k) = 0 for k ≥ n.

// Connection to operator calculus:
// The forward difference operator Δ and the shift operator E = I + Δ
// relate to Eulerian numbers through generating functions.
// Σ A(n,k) x^k = (1-x)^{n+1} * Σ k^n x^k

template<std::unsigned_integral W, int MaxN = 12>
struct EulerianCache {
    std::array<std::array<W, MaxN>, MaxN> A{};

    constexpr EulerianCache() noexcept {
        if (MaxN > 0) {
            A[0][0] = 1;
        }
        for (int n = 1; n < MaxN; ++n) {
            for (int k = 0; k < n; ++k) {
                W term1 = (k + 1 < MaxN) ? A[n-1][k] : W(0);
                W term2 = (k > 0) ? A[n-1][k-1] : W(0);
                A[n][k] = (k + 1) * term1 + (n - k) * term2;
            }
        }
    }
};

int main() {
    std::printf("=== 17 — Eulerian Numbers ===\n\n");

    using W = uint64_t;
    constexpr int MaxN = 10;

    constexpr auto cache = EulerianCache<W, MaxN>{};

    std::printf("Eulerian numbers A(n,k):\n");
    std::printf("  n\\k");
    for (int k = 0; k < MaxN - 1; ++k) {
        std::printf("  %3d", k);
    }
    std::printf("\n");

    for (int n = 0; n < MaxN; ++n) {
        std::printf("  %3d", n);
        for (int k = 0; k < n; ++k) {
            std::printf("  %3lu", (unsigned long)cache.A[n][k]);
        }
        std::printf("\n");
    }

    std::printf("\nProperties:\n");
    std::printf("  Σ_k A(n,k) = n!  (counts all permutations)\n");
    std::printf("  A(n,k) = A(n, n-1-k)  (symmetry)\n\n");

    // Verify sum = n!
    std::printf("Verification that Σ_k A(n,k) = n!:\n");
    W fact = 1;
    for (int n = 1; n < MaxN; ++n) {
        fact *= n;
        W sum = 0;
        for (int k = 0; k < n; ++k) {
            sum += cache.A[n][k];
        }
        std::printf("  n=%d: Σ A(n,k) = %lu, n! = %lu %s\n",
               n, (unsigned long)sum, (unsigned long)fact,
               (sum == fact) ? "✓" : "×");
    }

    // Connection to forward differences
    std::printf("\nConnection to forward difference calculus:\n");
    std::printf("  (1 - Δ)^{-(n+1)} = Σ_{k≥0} A(n+k,k) Δ^k\n");
    std::printf("  Where Δ is the forward difference operator.\n");
    std::printf("  This appears in umbral calculus / finite difference theory.\n\n");

    // Worpitzky's identity
    std::printf("Worpitzky's identity:\n");
    std::printf("  x^n = Σ_k A(n,k) * C(x+k, n)\n");
    std::printf("  where C(x+k, n) is the binomial coefficient.\n");

    return 0;
}
