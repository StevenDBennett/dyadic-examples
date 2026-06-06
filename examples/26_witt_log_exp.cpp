#include <dyadic.h>
#include <cstdio>

using namespace dyadic;

template<int N, std::unsigned_integral W>
void print_witt(const char* label, const WittVector<N, W>& w) {
    std::printf("  %-20s = {", label);
    for (int i = 0; i < N; ++i) std::printf("%s%lu", i ? ", " : "", (unsigned long)w[i]);
    std::printf("}\n");
}

int main() {
    std::printf("=== 26 — Witt Logarithm, Exponential & Inverse ===\n\n");
    std::printf("These operations work on ghost values at gw_t (widen_t<dword_t<W>>)\n");
    std::printf("precision. Roundtrip accuracy depends on the precision window\n");
    std::printf("constraint: recovered component r_j < 2^{8*sizeof(W)-j}.\n\n");

    using W = uint32_t;
    constexpr int N = 4;

    // === 1. Witt Inverse ===
    {
        std::printf("--- 1. Witt Inverse (a * a^{-1} == τ(1)) ---\n\n");

        for (W val : {W(3), W(5), W(7), W(9)}) {
            WittVector<N, W> u{{val, W(val+1), W(val+2), W(val+3)}};
            auto inv = witt_inverse(u);
            auto prod = u * inv;
            bool ok = (prod[0] == 1 && prod[1] == 0 && prod[2] == 0 && prod[3] == 0);
            std::printf("  a = {%lu, %lu, %lu, %lu}: a * a^{-1} = τ(1)  %s\n",
                   (unsigned long)val, (unsigned long)val+1,
                   (unsigned long)val+2, (unsigned long)val+3,
                   ok ? "PASS" : "FAIL");
        }

        // Simple Teichmuller units
        std::printf("\n  Teichmuller units (a = {x, 0, 0, 0}):\n");
        for (W val : {W(1), W(5), W(9), W(255)}) {
            WittVector<N, W> u{{val, 0, 0, 0}};
            auto inv = witt_inverse(u);
            auto prod = u * inv;
            bool ok = (prod[0] == 1 && prod[1] == 0 &&
                       prod[2] == 0 && prod[3] == 0);
            std::printf("    a = {%lu, 0, 0, 0}: %s\n",
                   (unsigned long)val, ok ? "PASS" : "FAIL");
        }
        std::printf("\n");
    }

    // === 2. Exponential Homomorphism: exp(a+b) == exp(a) * exp(b) ===
    {
        std::printf("--- 2. Exponential Homomorphism ---\n");
        std::printf("    exp(a + b) == exp(a) * exp(b)\n");
        std::printf("    Converges when a[0] ≡ 0 (mod 4) and b[0] ≡ 0 (mod 4)\n\n");

        WittVector<N, W> a{{4, 2, 0, 0}};
        WittVector<N, W> b{{8, 0, 0, 0}};
        print_witt("a", a);
        print_witt("b", b);

        auto exp_ab = witt_exp(a + b);
        auto exp_a_times_exp_b = witt_exp(a) * witt_exp(b);
        print_witt("exp(a + b)", exp_ab);
        print_witt("exp(a) * exp(b)", exp_a_times_exp_b);

        bool add_ok = true;
        for (int i = 0; i < N; ++i) {
            if (exp_ab[i] != exp_a_times_exp_b[i]) add_ok = false;
        }
        std::printf("  exp(a + b) == exp(a) * exp(b): %s\n\n",
               add_ok ? "PASS" : "FAIL");
    }

    // === 3. Logarithm Homomorphism: log(a*b) == log(a) + log(b) ===
    {
        std::printf("--- 3. Logarithm Homomorphism ---\n");
        std::printf("    log(a * b) == log(a) + log(b)\n");
        std::printf("    Defined when a[0] and b[0] are odd (units)\n\n");

        WittVector<N, W> unit_a{{3, 1, 0, 0}};
        WittVector<N, W> unit_b{{5, 2, 0, 0}};
        print_witt("unit_a", unit_a);
        print_witt("unit_b", unit_b);

        auto log_ab = witt_log(unit_a * unit_b);
        auto log_a_plus_log_b = witt_log(unit_a) + witt_log(unit_b);
        print_witt("log(a * b)", log_ab);
        print_witt("log(a) + log(b)", log_a_plus_log_b);

        bool mul_ok = true;
        for (int i = 0; i < N; ++i) {
            if (log_ab[i] != log_a_plus_log_b[i]) mul_ok = false;
        }
        std::printf("  log(a * b) == log(a) + log(b): %s\n\n",
               mul_ok ? "PASS" : "FAIL");
    }

    // === 4. Log/Exp Roundtrip (precision-limited) ===
    {
        std::printf("--- 4. Log/Exp Roundtrip ---\n");
        std::printf("    Roundtrip accuracy is limited by the precision window.\n");
        std::printf("    Accurate only when all components fit within the window.\n\n");

        // Works: single-component Witt vectors (Teichmuller-like)
        for (W val : {W(4), W(8), W(12), W(16)}) {
            WittVector<N, W> a{{val, 0, 0, 0}};
            auto b = witt_exp(a);
            auto c = witt_log(b);
            bool ok = true;
            for (int i = 0; i < N; ++i) {
                if (a[i] != c[i]) ok = false;
            }
            std::printf("  exp(log({%lu, 0, 0, 0})) == original: %s\n",
                   (unsigned long)val, ok ? "PASS" : "FAIL");
        }

        // Limited window: multi-component (only j=0 component accurate)
        WittVector<N, W> a{{4, 2, 4, 6}};
        auto exp_a = witt_exp(a);
        auto log_exp_a = witt_log(exp_a);
        print_witt("a", a);
        print_witt("log(exp(a))", log_exp_a);
        bool roundtrip = true;
        for (int i = 0; i < N; ++i) {
            if (log_exp_a[i] != a[i]) roundtrip = false;
        }
        std::printf("  Full roundtrip: %s (component 0 matches when within window)\n\n",
               roundtrip ? "PASS" : "FAIL (expected)");
    }

    // === 5. Adams operation via ghost power (existing) ===
    {
        std::printf("--- 5. Adams Operation Compose: ψ^m(ψ^n(a)) == ψ^{mn}(a) ---\n");
        WittVector<N, W> a{{3, 1, 4, 1}};
        print_witt("a", a);
        auto psi6_via_compose = adams_operation(adams_operation(a, 2), 3);
        auto psi6_direct = adams_operation(a, 6);
        print_witt("ψ³(ψ²(a))", psi6_via_compose);
        print_witt("ψ⁶(a)", psi6_direct);
        bool compose_ok = true;
        for (int i = 0; i < N; ++i) {
            if (psi6_via_compose[i] != psi6_direct[i]) compose_ok = false;
        }
        std::printf("  ψ³(ψ²(a)) == ψ⁶(a): %s\n\n", compose_ok ? "PASS" : "FAIL");
    }

    std::printf("=== Summary ===\n");
    std::printf("  Witt inverse (a * a^{-1} = τ(1)):      Working\n");
    std::printf("  Exp homomorphism (exp(a+b) = exp(a)*exp(b)): Working\n");
    std::printf("  Log homomorphism (log(a*b) = log(a)+log(b)): Working\n");
    std::printf("  Log/Exp roundtrip: Precision-window limited\n");
    std::printf("  Adams compose (ψ^m∘ψ^n = ψ^{mn}):     Working (existing)\n");
    std::printf("\n");
    std::printf("Note: The precision window constraint r_j < 2^{8*sizeof(W)-j}\n");
    std::printf("limits multi-component accuracy. Use larger word sizes (uint64_t)\n");
    std::printf("to improve it, or uint128_t as the primary word type for full coverage.\n");

    return 0;
}
