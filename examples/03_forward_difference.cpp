#include "dyadic.h"
#include <cstdio>

template<int N, std::unsigned_integral W>
void print_poly(const char* label, const dyadic::Polynomial<N, W, dyadic::MonomialBasis>& p) {
    std::printf("%s = {", label);
    for (int i = 0; i < N; ++i) std::printf("%s%lu", i ? ", " : "", (unsigned long)p[i]);
    std::printf("}\n");
}

int main() {
    using namespace dyadic;

    std::printf("=== 03 — Forward Difference Δ ===\n\n");
    std::printf("Δ(P)(t) = P(t+1) − P(t)\n\n");

    Polynomial<5, uint64_t, MonomialBasis> p{{1, 2, 3, 4, 5}};
    print_poly("P", p);

    auto dp = forward_difference(p);
    print_poly("ΔP", dp);
    std::printf("\n  Verify: P(t+1) − P(t):\n");
    std::printf("    P(5)  = %lu,  P(4)  = %lu,  diff = %lu\n", p.eval(5), p.eval(4), p.eval(5) - p.eval(4));
    std::printf("    P(10) = %lu, P(9)  = %lu,  diff = %lu\n", p.eval(10), p.eval(9), p.eval(10) - p.eval(9));
    std::printf("    ΔP(4) = %lu\n\n", dp.eval(4));

    std::printf("Δ in Falling Factorial basis: Δ((t)_n) = n·(t)_{n−1}\n");
    auto ff = change_basis<FallingFactorialBasis>(p);
    auto d_ff = forward_difference(ff);
    std::printf("  Δ(P) in FF basis = {");
    for (int i = 0; i < 4; ++i) std::printf("%s%lu", i ? ", " : "", (unsigned long)d_ff[i]);
    std::printf("}\n\n");

    auto d1 = formal_derivative(p);
    auto d2 = formal_derivative(d1);
    auto d3 = formal_derivative(d2);
    auto d4 = formal_derivative(d3);

    std::printf("Identity: Δ = e^D − I  (truncated series)\n");
    std::printf("  ΔP = DP + D²P/2! + D³P/3! + ...\n\n");
    std::printf("  deg(DP)=%d, deg(D²P)=%d, deg(D³P)=%d, deg(D⁴P)=%d\n",
        d1.max_degree, d2.max_degree, d3.max_degree, d4.max_degree);
    std::printf("  DP           = {%lu, %lu, %lu, %lu}\n", d1[0], d1[1], d1[2], d1[3]);
    std::printf("  D²P/2!       = {%lu, %lu, %lu}\n", d2[0]/2, d2[1]/2, d2[2]/2);
    std::printf("  D³P/3!       = {%lu, %lu}\n", d3[0]/6, d3[1]/6);
    std::printf("  D⁴P/4!       = {%lu}\n", d4[0]/24);
    uint64_t s0 = d1[0] + d2[0]/2 + d3[0]/6 + d4[0]/24;
    uint64_t s1 = d1[1] + d2[1]/2 + d3[1]/6;
    uint64_t s2 = d1[2] + d2[2]/2;
    uint64_t s3 = d1[3];
    std::printf("  Sum (e^D-I)  = {%lu, %lu, %lu, %lu}\n", s0, s1, s2, s3);
    std::printf("  ΔP           = {%lu, %lu, %lu, %lu}\n", dp[0], dp[1], dp[2], dp[3]);
    std::printf("  Match: %s\n\n",
        (dp[0] == s0 && dp[1] == s1 && dp[2] == s2 && dp[3] == s3) ? "yes" : "no");

    return 0;
}
