#include <dyadic.h>
#include <cstdio>

template<int N, std::unsigned_integral W>
void print_witt(const char* label, const dyadic::WittVector<N, W>& w) {
    std::printf("  %-16s = {", label);
    for (int i = 0; i < N; ++i) std::printf("%s%lu", i ? ", " : "", (unsigned long)w[i]);
    std::printf("}\n");
}

int main() {
    using namespace dyadic;

    std::printf("=== 07 — Adams Operations & Teichmüller Lifts ===\n\n");

    WittVector<3, uint32_t> a{{2, 3, 5}};
    print_witt("a", a);
    std::printf("\n");

    std::printf("Adams operation ψⁿ: ghost_j(ψⁿ(a)) = ghost_j(a)^n\n");
    auto psi1 = adams_operation(a, 1);
    print_witt("ψ¹(a) (identity)", psi1);
    bool id_check = true;
    for (int i = 0; i < 3; ++i) if (psi1[i] != a[i]) id_check = false;
    std::printf("  ψ¹(a) == a: %s\n\n", id_check ? "yes" : "no");

    auto psi2 = adams_operation(a, 2);
    print_witt("ψ²(a)", psi2);
    auto ga = a.ghost_vector();
    auto gpsi2 = psi2.ghost_vector();
    std::printf("  Ghost check: G(ψ²(a))_j == G(a)_j² for all j:\n");
    bool ghost_ok = true;
    for (int j = 0; j < 3; ++j) {
        bool ok = (gpsi2[j] == ga[j] * ga[j]);
        std::printf("    j=%d: %u² = %u  vs  %u  %s\n", j, ga[j], ga[j] * ga[j], gpsi2[j], ok ? "✓" : "✗");
        if (!ok) ghost_ok = false;
    }
    std::printf("  Result: %s\n\n", ghost_ok ? "PASS" : "FAIL");

    auto psi3 = adams_operation(a, 3);
    print_witt("ψ³(a)", psi3);

    std::printf("\nψ³ ∘ ψ² = ψ⁶  (Adams compose):\n");
    auto psi6 = adams_operation(a, 6);
    auto psi3_psi2 = adams_operation(adams_operation(a, 2), 3);
    print_witt("ψ⁶(a)", psi6);
    print_witt("ψ³(ψ²(a))", psi3_psi2);
    bool compose_ok = true;
    for (int i = 0; i < 3; ++i) if (psi6[i] != psi3_psi2[i]) compose_ok = false;
    std::printf("  ψ⁶(a) == ψ³(ψ²(a)): %s\n\n", compose_ok ? "PASS" : "FAIL");

    std::printf("Teichmüller lift: τ(x) = (x, 0, 0, ...)\n");
    auto tau_7 = teichmueller_lift<4, uint32_t>(7);
    print_witt("τ(7)", tau_7);

    auto tau_3 = teichmueller_lift<4, uint32_t>(3);
    auto tau_5 = teichmueller_lift<4, uint32_t>(5);
    auto tau_15 = teichmueller_lift<4, uint32_t>(15);
    std::printf("\nMultiplicative: τ(a) × τ(b) = τ(a·b)\n");
    print_witt("τ(3)", tau_3);
    print_witt("τ(5)", tau_5);
    print_witt("τ(3)×τ(5)", tau_3 * tau_5);
    print_witt("τ(15) (expected)", tau_15);

    bool teich_ok = true;
    for (int i = 0; i < 4; ++i) if ((tau_3 * tau_5)[i] != tau_15[i]) teich_ok = false;
    std::printf("  τ(3)×τ(5) == τ(15): %s\n", teich_ok ? "PASS" : "FAIL");

    return 0;
}
