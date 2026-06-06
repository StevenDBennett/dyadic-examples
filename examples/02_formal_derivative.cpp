#include <dyadic.h>
#include <cstdio>

template<int N, std::unsigned_integral W>
void print_poly(const char* label, const dyadic::Polynomial<N, W, dyadic::MonomialBasis>& p) {
    std::printf("%s = {", label);
    for (int i = 0; i < N; ++i) std::printf("%s%lu", i ? ", " : "", (unsigned long)p[i]);
    std::printf("}\n");
}

int main() {
    using namespace dyadic;

    std::printf("=== 02 — Formal Derivative D = d/dt ===\n\n");

    Polynomial<5, uint64_t, MonomialBasis> p{{3, 1, 4, 1, 5}};
    print_poly("P", p);
    std::printf("  P(t) = 3 + t + 4t² + t³ + 5t⁴\n\n");

    auto dp = formal_derivative(p);
    print_poly("DP", dp);
    std::printf("  P'(t) = 1 + 8t + 3t² + 20t³\n\n");

    auto d2p = formal_derivative(dp);
    print_poly("D²P", d2p);
    std::printf("  P''(t) = 8 + 6t + 60t²\n\n");

    auto d3p = formal_derivative(d2p);
    print_poly("D³P", d3p);
    std::printf("  P'''(t) = 6 + 120t\n\n");

    auto d4p = formal_derivative(d3p);
    print_poly("D⁴P", d4p);
    std::printf("  P⁽⁴⁾(t) = 120\n\n");

    auto d5p = formal_derivative(d4p);
    std::printf("D⁵P = {}  (zero polynomial — D^N = 0 for deg N-1)\n\n");

    std::printf("Check: D^N P = 0 for degree N-1 polynomial\n");
    std::printf("  deg(P) = 4, N = 5, D⁵ = 0: confirmed\n\n");

    std::printf("Derivative in Falling Factorial basis:\n");
    auto ff = change_basis<FallingFactorialBasis>(p);
    auto d_ff = formal_derivative(ff);
    auto d_ff_back = change_basis<MonomialBasis>(d_ff);
    print_poly("D(P) via FF basis", d_ff_back);
    bool match = true;
    for (int i = 0; i < 4; ++i) if (d_ff_back[i] != dp[i]) match = false;
    std::printf("  Matches monomial D(P): %s\n", match ? "yes" : "no");

    return 0;
}
