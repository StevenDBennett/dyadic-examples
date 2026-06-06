#include "dyadic.h"
#include <cstdio>

template<int N, std::unsigned_integral W>
void print_poly(const char* label, const dyadic::Polynomial<N, W, dyadic::MonomialBasis>& p) {
    std::printf("  %-14s = {", label);
    for (int i = 0; i < N; ++i) std::printf("%s%lu", i ? ", " : "", (unsigned long)p[i]);
    std::printf("}\n");
}

int main() {
    using namespace dyadic;

    std::printf("=== 05 — Indefinite Sum Σ = Δ⁻¹ ===\n\n");
    std::printf("Σ is the inverse of forward difference: Σ(ΔP) = P,  Δ(ΣP) = P\n\n");

    Polynomial<4, uint64_t, MonomialBasis> p{{1, 2, 3, 4}};
    print_poly("P (monomial)", p);
    std::printf("\n");

    auto delta_p = forward_difference(p);
    print_poly("ΔP", delta_p);
    std::printf("\n");

    auto sum_delta = indefinite_sum(delta_p);
    print_poly("Σ(ΔP)", sum_delta);
    std::printf("  Σ(ΔP) = P − P(0)  (kernel of Σ = constants)\n");
    bool kernel = (sum_delta[0] == 0);
    for (int i = 1; i < 4; ++i) if (sum_delta[i] != p[i]) kernel = false;
    std::printf("  Σ(ΔP) == P − P(0)?  %s\n\n", kernel ? "yes" : "no");

    auto delta_sum = forward_difference(indefinite_sum(p));
    print_poly("Δ(ΣP)", delta_sum);
    bool id2 = true;
    for (int i = 0; i < 4; ++i) if (delta_sum[i] != p[i]) id2 = false;
    std::printf("  Δ(ΣP) == P?  %s\n", id2 ? "yes (Δ = Σ⁻¹)" : "no (limited by 2-division precision)");
    std::printf("  Note: Exact inverse requires p[i] divisible by 2^{v2(i+1)}\n");
    std::printf("  P_FF[1] = 9 is odd → division by 2 is inexact in ℤ₂\n\n");

    std::printf("In Falling Factorial basis: Σ((t)_{n-1}) = (t)_n / n\n");
    auto ff = change_basis<FallingFactorialBasis>(p);
    std::printf("  P in FF basis = {");
    for (int i = 0; i < 4; ++i) std::printf("%s%lu", i ? ", " : "", (unsigned long)ff[i]);
    std::printf("}\n");
    auto sum_ff = indefinite_sum(ff);
    auto sum_ff_mono = change_basis<MonomialBasis>(sum_ff);
    print_poly("ΣP (via FF)", sum_ff_mono);
    auto sum_mono = indefinite_sum(p);
    print_poly("ΣP (monomial)", sum_mono);
    bool match = true;
    for (int i = 0; i < 5; ++i) if (sum_ff_mono[i] != sum_mono[i]) match = false;
    std::printf("  Both methods match: %s\n", match ? "yes" : "no");

    return 0;
}
