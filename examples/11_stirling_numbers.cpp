#include <dyadic.h>
#include <cstdio>

int main() {
    using namespace dyadic;

    std::printf("=== 11 — Stirling Numbers ===\n\n");

    constexpr int MAX = 8;

    std::printf("Stirling numbers of the second kind S(n,k):\n");
    std::printf("  n\\k");
    for (int k = 0; k <= MAX; ++k) std::printf("%6d", k);
    std::printf("\n");
    for (int n = 0; n <= MAX; ++n) {
        std::printf("  %2d ", n);
        for (int k = 0; k <= n; ++k) {
            std::printf("%6lu", stirling_2<uint64_t>(n, k));
        }
        std::printf("\n");
    }
    std::printf("\n");

    std::printf("Stirling numbers of the first kind s(n,k) (signed, as 2-adic):\n");
    std::printf("  n\\k");
    for (int k = 0; k <= MAX; ++k) std::printf("%6d", k);
    std::printf("\n");
    for (int n = 0; n <= MAX; ++n) {
        std::printf("  %2d ", n);
        for (int k = 0; k <= n; ++k) {
            uint64_t s = stirling_1<uint64_t>(n, k);
            std::printf("%6lu", s);
        }
        std::printf("\n");
    }
    std::printf("\n");

    std::printf("Key identities verified:\n");
    bool s2_recurrence = true;
    for (int n = 2; n <= 6; ++n) {
        for (int k = 1; k < n; ++k) {
            auto lhs = stirling_2<uint64_t>(n, k);
            auto rhs = stirling_2<uint64_t>(n-1, k-1) + static_cast<uint64_t>(k) * stirling_2<uint64_t>(n-1, k);
            if (lhs != rhs) { s2_recurrence = false; break; }
        }
    }
    std::printf("  S(n,k) = S(n-1,k-1) + k·S(n-1,k): %s\n", s2_recurrence ? "PASS" : "FAIL");

    bool s1_recurrence = true;
    for (int n = 2; n <= 6; ++n) {
        for (int k = 1; k < n; ++k) {
            auto lhs = stirling_1<uint64_t>(n, k);
            auto rhs = stirling_1<uint64_t>(n-1, k-1) - static_cast<uint64_t>(n-1) * stirling_1<uint64_t>(n-1, k);
            if (lhs != rhs) { s1_recurrence = false; break; }
        }
    }
    std::printf("  s(n,k) = s(n-1,k-1) − (n-1)·s(n-1,k): %s\n", s1_recurrence ? "PASS" : "FAIL");

    std::printf("\nUnsigned Stirling numbers of the first kind |s(n,k)|:\n");
    std::printf("  n\\k");
    for (int k = 0; k <= MAX; ++k) std::printf("%6d", k);
    std::printf("\n");
    for (int n = 0; n <= MAX; ++n) {
        std::printf("  %2d ", n);
        for (int k = 0; k <= n; ++k) {
            std::printf("%6lu", stirling_1_unsigned<uint64_t>(n, k));
        }
        std::printf("\n");
    }

    return 0;
}
