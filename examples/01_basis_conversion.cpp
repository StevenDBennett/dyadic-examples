#include <dyadic.h>
#include <cstdio>

int main() {
    using namespace dyadic;

    std::printf("=== 01 — Basis Conversions ===\n\n");

    Polynomial<5, uint64_t, MonomialBasis> p{{1, 2, 3, 4, 5}};

    std::printf("Monomial basis:\n  P = {");
    for (int i = 0; i < 5; ++i) std::printf("%s%lu", i ? ", " : "", p[i]);
    std::printf("}  →  P(t) = 1 + 2t + 3t² + 4t³ + 5t⁴\n\n");

    auto ff = change_basis<FallingFactorialBasis>(p);
    std::printf("Falling Factorial basis:\n  P = {");
    for (int i = 0; i < 5; ++i) std::printf("%s%lu", i ? ", " : "", ff[i]);
    std::printf("}\n  P(t) = ");
    for (int i = 0; i < 5; ++i)
        if (ff[i]) std::printf("%s%lu·(t)%s%d", i ? " + " : "", ff[i], "", i);
    std::printf("\n\n");

    auto back = change_basis<MonomialBasis>(ff);
    bool ok = true;
    for (int i = 0; i < 5; ++i) if (p[i] != back[i]) ok = false;
    std::printf("Roundtrip Mono → FF → Mono: %s\n\n", ok ? "PASS" : "FAIL");

    auto taylor = change_basis<TaylorBasis>(p);
    std::printf("Taylor basis (T_k = k! · FF_k):\n  T = {");
    for (int i = 0; i < 5; ++i) std::printf("%s%lu", i ? ", " : "", taylor[i]);
    std::printf("}\n  T_k = {1, 2, 6, 24, 120} × FF_k\n");
    std::printf("  Taylor roundtrip check: %s\n\n",
        verify::check_taylor_roundtrip_precision(p) ? "PASS" : "FAIL (precision loss)");

    Polynomial<6, uint8_t, MonomialBasis> big{{0, 0, 0, 0, 0, 255}};
    std::printf("Precision window (uint8_t):\n");
    std::printf("  P(t) = 255 t⁵  →  T_5 = 5!·255 = %d > 256\n", 120 * 255);
    std::printf("  Taylor roundtrip: %s\n\n",
        verify::check_taylor_roundtrip_precision(big) ? "PASS" : "FAIL (expected — coefficient too large)");

    auto p_at_3 = p.eval(3);
    auto ff_at_3 = eval(ff, uint64_t(3));
    auto t_at_3 = eval(taylor, uint64_t(3));
    std::printf("Evaluation at t=3:\n  Monomial:        P(3) = %lu\n", p_at_3);
    std::printf("  FallingFacto:   P(3) = %lu\n", ff_at_3);
    std::printf("  Taylor:         P(3) = %lu\n", t_at_3);
    std::printf("  All match: %s\n", (p_at_3 == ff_at_3 && ff_at_3 == t_at_3) ? "yes" : "no");

    return 0;
}
