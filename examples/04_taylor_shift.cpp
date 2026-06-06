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

    std::printf("=== 04 — Taylor Shift P(t) → P(t + δ) ===\n\n");

    Polynomial<5, uint64_t, MonomialBasis> p{{1, 2, 3, 4, 5}};
    print_poly("P(t) (monomial)", p);
    std::printf("  P(t) = 1 + 2t + 3t² + 4t³ + 5t⁴\n\n");

    for (uint64_t delta : {uint64_t(1), uint64_t(2), uint64_t(10)}) {
        auto shifted = taylor_shift(p, delta);
        std::printf("Shift by δ = %lu:\n", delta);
        print_poly("P(t+δ)", shifted);

        uint64_t direct = p.eval(static_cast<uint64_t>(7 + delta));
        uint64_t shifted_val = shifted.eval(7);
        std::printf("  P(%lu+%lu) = P(%lu) by eval shifted = %lu\n",
            (unsigned long)7, delta, (unsigned long)(7 + delta), shifted_val);
        std::printf("  P(%lu) direct = %lu  %s\n\n",
            (unsigned long)(7 + delta), direct,
            direct == shifted_val ? "(match)" : "(MISMATCH)");
    }

    std::printf("Shift in Falling Factorial basis:\n");
    auto ff = change_basis<FallingFactorialBasis>(p);
    auto shifted_ff = taylor_shift(ff, uint64_t(3));
    auto shifted_ff_back = change_basis<MonomialBasis>(shifted_ff);
    print_poly("P(t+3) via FF shift", shifted_ff_back);

    auto shifted_mono = taylor_shift(p, uint64_t(3));
    bool match = true;
    for (int i = 0; i < 5; ++i) if (shifted_ff_back[i] != shifted_mono[i]) match = false;
    std::printf("  Matches monomial shift: %s\n", match ? "yes" : "no");

    return 0;
}
