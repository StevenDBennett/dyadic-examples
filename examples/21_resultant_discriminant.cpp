#include "dyadic.h"
#include <cstdint>
#include <cstdio>

using namespace dyadic;

// Resultant of two polynomials: res(P,Q) = Π_{P(a)=0} Q(a)^{mult(a)}
// Using the Sylvester matrix determinant: res(P,Q) = det(Sylvester(P,Q))

template<int N, int M, std::unsigned_integral W>
constexpr W resultant_monomial(const Polynomial<N, W, MonomialBasis>& P,
                                const Polynomial<M, W, MonomialBasis>& Q) {
    // For monomial basis, we use formulas for small cases
    // General Sylvester matrix would be (N+M-2) x (N+M-2)
    // This is a simplified version for demonstration

    if constexpr (N == 2 && M == 2) {
        // P = a0 + a1*x, Q = b0 + b1*x
        // Resultant = a1*b0 - a0*b1 (both linear)
        return P[1] * Q[0] - P[0] * Q[1];
    } else if constexpr (N == 3 && M == 2) {
        // P = a0 + a1*x + a2*x^2, Q = b0 + b1*x
        // Resultant = a2*b0^2 - a1*b0*b1 + a0*b1^2
        return P[2] * Q[0] * Q[0] - P[1] * Q[0] * Q[1] + P[0] * Q[1] * Q[1];
    } else {
        // Fallback: return 0 (not implemented for general degree)
        return W(0);
    }
}

// Discriminant: disc(P) = (-1)^{n(n-1)/2} / leading_coeff * res(P, P')
template<int N, std::unsigned_integral W>
constexpr W discriminant(const Polynomial<N, W, MonomialBasis>& P) {
    if constexpr (N <= 2) return W(0);  // No discriminant for degree < 2

    auto P_prime = formal_derivative(P);
    // For degree 2: P = a + b*x + c*x^2, disc = b^2 - 4ac
    if constexpr (N == 3) {
        return P[1] * P[1] - 4 * P[2] * P[0];
    }
    return W(0);
}

int main() {
    std::printf("=== 21 — Resultant and Discriminant ===\n\n");

    using W = uint64_t;

    // Example 1: Linear polynomials
    // P(x) = 2 + 3x, Q(x) = 5 + 7x
    Polynomial<2, W, MonomialBasis> P1{{2, 3}};
    Polynomial<2, W, MonomialBasis> Q1{{5, 7}};
    W res1 = resultant_monomial(P1, Q1);
    std::printf("Resultant(2+3x, 5+7x) = %lu\n", (unsigned long)res1);
    std::printf("  Expected: 3*5 - 2*7 = 15 - 14 = 1\n");
    std::printf("  Match: %s\n\n", (res1 == 1) ? "PASS" : "FAIL");

    // Example 2: Quadratic and linear
    // P(x) = 1 + 2x + 3x^2, Q(x) = 4 + 5x
    Polynomial<3, W, MonomialBasis> P2{{1, 2, 3}};
    Polynomial<2, W, MonomialBasis> Q2{{4, 5}};
    W res2 = resultant_monomial(P2, Q2);
    std::printf("Resultant(1+2x+3x^2, 4+5x) = %lu\n", (unsigned long)res2);
    std::printf("  Expected: 3*16 - 2*20 + 1*25 = 48 - 40 + 25 = 33\n");
    std::printf("  Match: %s\n\n", (res2 == 33) ? "PASS" : "FAIL");

    // Example 3: Discriminant of a quadratic
    // P(x) = 1 + 2x + 3x^2
    // disc = b^2 - 4ac = 4 - 12 = -8
    W disc = discriminant(P2);
    std::printf("Discriminant(1+2x+3x^2) = %lu\n", (unsigned long)disc);
    std::printf("  Expected: 4 - 12 = -8 (mod 2^64: %lu)\n",
           (unsigned long)(W(0) - 8));
    std::printf("  (In ℤ₂, -8 wraps to a large number)\n\n");

    std::printf("Applications:\n");
    std::printf("  Resultant = 0  ⇔  P and Q have a common root\n");
    std::printf("  Discriminant = 0  ⇔  P has a repeated root\n\n");

    std::printf("Sylvester matrix for P = Σ a_i x^i and Q = Σ b_j x^j:\n");
    std::printf("  [ a_n    a_{n-1}  ...  a_0    0    ...   0  ]\n");
    std::printf("  [  0      a_n     ...  a_1   a_0   ...   0  ]\n");
    std::printf("  [  ...                               ...      ]\n");
    std::printf("  [  0       0      ...   0     b_m   ...   b_0 ]\n");
    std::printf("  Resultant = determinant of this (n+m) × (n+m) matrix\n");

    return 0;
}
