#include "dyadic.h"
#include <cstdio>

template<int N, std::unsigned_integral W>
void print_poly(const char* label, const dyadic::Polynomial<N, W, dyadic::MonomialBasis>& p) {
    std::printf("  %-20s = {", label);
    for (int i = 0; i < N; ++i) std::printf("%s%lu", i ? ", " : "", (unsigned long)p[i]);
    std::printf("}\n");
}

int main() {
    using namespace dyadic;

    std::printf("=== 08 — Power Series Composition & Reversion ===\n\n");

    Polynomial<5, uint64_t, MonomialBasis> P{{0, 1, 1, 0, 0}};
    Polynomial<4, uint64_t, MonomialBasis> Q{{0, 1, 2, 0}};
    print_poly("P(t) = t + t²", P);
    print_poly("Q(t) = t + 2t²", Q);
    std::printf("\n");

    auto PQ = compose(P, Q);
    int PQ_deg = 4;
    print_poly("P(Q(t))", PQ);
    std::printf("  Expected degree: (4)*(2)+1 = 9, computed degree = %d\n\n", PQ_deg);

    std::printf("Evaluating P(Q(3)):\n");
    uint64_t q3 = Q.eval(3);
    uint64_t pq3 = P.eval(q3);
    uint64_t pq_direct = PQ.eval(3);
    std::printf("  Q(3) = %lu, P(Q(3)) = %lu, (P∘Q)(3) = %lu  %s\n\n",
        q3, pq3, pq_direct, pq3 == pq_direct ? "✓" : "✗");

    std::printf("Reversion (Lagrange inversion):\n");
    std::printf("  Find R such that P(R(t)) = t\n\n");

    Polynomial<6, uint64_t, MonomialBasis> S{{0, 1, 1}};
    print_poly("S(t) = t + t²", S);

    auto R = reversion(S);
    print_poly("R(t) (reversion)", R);
    std::printf("  Expected: R(t) = t − t² + 2t³ − 5t⁴ + 14t⁵ − 42t⁶ + …\n");
    std::printf("  (Negative coefficients stored as 2^64 − value in ℤ₂)\n\n");

    auto SR = compose(S, R);
    print_poly("S(R(t))", SR);
    std::printf("  Should be: {0, 1, 0, 0, 0, ...} = t\n");
    bool is_t = (SR[0] == 0 && SR[1] == 1);
    for (int i = 2; i < 6; ++i) if (SR[i] != 0) { is_t = false; break; }
    std::printf("  S(R(t)) == t: %s\n\n", is_t ? "PASS" : "FAIL");

    std::printf("Compose with wider polynomials:\n");
    Polynomial<5, uint64_t, MonomialBasis> A{{0, 2, 3, 4, 0}};
    Polynomial<3, uint64_t, MonomialBasis> B{{0, 1, 1}};
    print_poly("A(t) = 2t+3t²+4t³", A);
    print_poly("B(t) = t+t²", B);
    auto AB = compose(A, B);
    print_poly("A(B(t))", AB);

    return 0;
}
