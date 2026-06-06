#include "dyadic.h"
#include <cstdio>

template<int N, std::unsigned_integral W, typename B>
void print_poly_basis(const char* label, const dyadic::Polynomial<N, W, B>& p) {
    std::printf("  %-20s = {", label);
    for (int i = 0; i < N; ++i) std::printf("%s%lu", i ? ", " : "", (unsigned long)p[i]);
    std::printf("}\n");
}

int main() {
    using namespace dyadic;

    std::printf("=== 12 — Polynomial Arithmetic ===\n\n");

    Polynomial<4, uint64_t, MonomialBasis> p{{1, 2, 3, 4}};
    Polynomial<4, uint64_t, MonomialBasis> q{{5, 6, 7, 8}};

    print_poly_basis("P(t)", p);
    print_poly_basis("Q(t)", q);
    std::printf("\n");

    auto sum = p + q;
    print_poly_basis("P + Q", sum);
    std::printf("  Coefficient-wise addition in ℤ₂[[t]]\n\n");

    auto diff = p - q;
    print_poly_basis("P − Q", diff);
    std::printf("  (Wrapping: 1-5 = %lu = 2⁶⁴-4, etc.)\n\n", p[0] - q[0]);

    auto prod = p * q;
    print_poly_basis("P × Q", prod);
    std::printf("  Degree: %d (expected: 3 + 3 = 6)\n\n", prod.degree);

    std::printf("Evaluation in all three bases:\n");
    for (uint64_t t_val : {uint64_t(0), uint64_t(1), uint64_t(2), uint64_t(5), uint64_t(10)}) {
        uint64_t v_mono = p.eval(t_val);
        auto ff = change_basis<FallingFactorialBasis>(p);
        uint64_t v_ff = ff.eval(t_val);
        auto taylor = change_basis<TaylorBasis>(p);
        uint64_t v_tay = taylor.eval(t_val);
        std::printf("  P(%2lu) = %5lu  (mono=%lu, ff=%lu, taylor=%lu)%s\n",
            t_val, v_mono, v_mono, v_ff, v_tay,
            (v_mono == v_ff && v_ff == v_tay) ? "" : " MISMATCH");
    }
    std::printf("\n");

    std::printf("Polynomial multiplication by scalar (via coefficient-wise):\n");
    Polynomial<4, uint64_t, MonomialBasis> scaled;
    for (int i = 0; i < 4; ++i) scaled[i] = p[i] * uint64_t(3);
    print_poly_basis("3·P(t)", scaled);
    std::printf("  Verify: P(2) = %lu, 3·P(2) = %lu\n", p.eval(2), scaled.eval(2));
    std::printf("  3·P(2) == (3·P)(2): %s\n\n", scaled.eval(2) == 3 * p.eval(2) ? "yes" : "no");

    std::printf("Falling Factorial arithmetic:\n");
    auto ff_p = change_basis<FallingFactorialBasis>(p);
    auto ff_q = change_basis<FallingFactorialBasis>(q);
    auto ff_sum = ff_p + ff_q;
    auto ff_sum_back = change_basis<MonomialBasis>(ff_sum);
    print_poly_basis("P+Q (FF→Mono)", ff_sum_back);
    bool match = true;
    for (int i = 0; i < 4; ++i) if (ff_sum_back[i] != sum[i]) match = false;
    std::printf("  Matches monomial P+Q: %s\n", match ? "yes" : "no");

    return 0;
}
