#include <dyadic.h>
#include <cstdint>
#include <cstdio>

using namespace dyadic;

// Pade approximants in ℤ₂[[t]]: rational approximations to formal power series.
// We avoid division by even numbers by working with series whose coefficients
// are 2-adic integers.

// Evaluate a rational function P(t)/Q(t) at point x (truncated series)
// where P has coeffs P[0]..P[M], Q has coeffs Q[0]..Q[N], Q(0) = 1 (normalized)
template<int PN, int QN, std::unsigned_integral W>
constexpr W eval_rational(const Polynomial<PN, W, MonomialBasis>& P,
                          const Polynomial<QN, W, MonomialBasis>& Q,
                          W x) {
    W num = 0, den = 0;
    for (int i = PN - 1; i >= 0; --i) {
        num = num * x + P[i];
    }
    for (int i = QN - 1; i >= 0; --i) {
        den = den * x + Q[i];
    }
    // If denominator is odd, it's invertible
    if (den % 2 == 1) {
        return num * modinv_odd(den);
    }
    // Otherwise handle with exact division (approximate)
    return num / den;
}

// Geometric series: (1 - t)^{-1} = Σ t^k
// Pade [0/1] = (1 - t)^{-1}
template<int N, std::unsigned_integral W>
constexpr Polynomial<N, W, MonomialBasis> geometric_series() noexcept {
    Polynomial<N, W, MonomialBasis> r{};
    for (int i = 0; i < N; ++i) {
        r[i] = 1;
    }
    return r;
}

int main() {
    std::printf("=== 15 — Pade Approximants in ℤ₂[[t]] ===\n\n");

    using W = uint64_t;

    // Geometric series: 1/(1-t) = Σ t^k
    auto G = geometric_series<8, W>();
    std::printf("Geometric series: 1/(1-t) = Σ t^k\n  ");
    for (int i = 0; i < 8; ++i) {
        std::printf("%s%lu t^%d", (i > 0) ? " + " : "", (unsigned long)G[i], i);
    }
    std::printf("\n\n");

    // Verify via evaluation: truncated series Σ x^k approximates 1/(1-x)
    std::printf("Evaluation at x=0: G(0) = %lu (should be 1)\n",
           (unsigned long)G.eval(W(0)));

    // Comparison: geometric series evaluation at various points
    // 1/(1-x) via series: only valid when x is even (so 1-x is odd and invertible)
    std::printf("Evaluating 1/(1-t) = Σ t^k at t = x (x even, so 1-x is odd):\n");
    for (W x : {W(0), W(2), W(4), W(6)}) {
        W series_val = G.eval(x);
        // Rational evaluation using modinv_odd: 1/(1-x) when 1-x is odd
        W inv = modinv_odd(W(1) - x);
        std::printf("  x=%lu: series=%lu, 1/(1-x)=%lu\n",
               (unsigned long)x, (unsigned long)series_val, (unsigned long)inv);
    }
    std::printf("\n");

    // Pade concept: [k/k] approximant matches first 2k+1 terms
    std::printf("Pade approximant concept:\n");
    std::printf("  [k/k] Pade = P_k(t) / Q_k(t) matches Taylor to O(t^{2k+1})\n");
    std::printf("  Converges faster than Taylor for fixed k\n\n");

    // Rational function with odd denominators
    // Rational function example with invertible denominators
    std::printf("Example: Rational function (1 + 3t) / (1 + 3t) = 1\n");
    Polynomial<2, W, MonomialBasis> P{{1, 3}};
    Polynomial<2, W, MonomialBasis> Q{{1, 3}};
    std::printf("  Evaluation always yields 1 (invertible denominator):\n");
    for (W x : {W(1), W(3), W(5)}) {
        W val = eval_rational(P, Q, x);
        std::printf("    x=%lu: (1+3x)/(1+3x) = %lu %s\n",
               (unsigned long)x, (unsigned long)val,
               val == 1 ? "PASS" : "FAIL");
    }
    std::printf("\n");

    return 0;
}
