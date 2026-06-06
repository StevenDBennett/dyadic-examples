#include <dyadic.h>
#include <cstdio>

template<int N, std::unsigned_integral W>
void print_poly(const char* label, const dyadic::Polynomial<N, W, dyadic::MonomialBasis>& p) {
    std::printf("  %-16s = {", label);
    for (int i = 0; i < N; ++i) std::printf("%s%lu", i ? ", " : "", (unsigned long)p[i]);
    std::printf("}\n");
}

int main() {
    using namespace dyadic;

    std::printf("=== 31 — Discrete Hedging (Δ/Σ in Finance) ===\n\n");

    // Stock price path modelled in cents (integer cents = $/100):
    // S(t) = 10000 + 200t + 10t²  (t = 0, 1, 2, ...)
    // S(0) = $100.00, S(1) = $102.10, S(2) = $104.40, ...
    Polynomial<4, uint64_t, MonomialBasis> S{{10000, 200, 10, 0}};
    print_poly("S(t) (cents)", S);

    std::printf("\n--- Daily Returns via Forward Difference ---\n");
    std::printf("  Δ(S)(t) = S(t+1) − S(t)\n");
    auto ret = forward_difference(S);
    print_poly("ΔS (returns)", ret);

    std::printf("\n  S(t)  (cents):\n");
    for (int t = 0; t < 6; ++t) {
        uint64_t s = S.eval(t);
        std::printf("    t=%d  $%lu.%02lu\n", t,
            (unsigned long)(s / 100), (unsigned long)(s % 100));
    }
    std::printf("\n  Daily return R(t) = Δ(S)(t) (cents):\n");
    for (int t = 0; t < 3; ++t)
        std::printf("    R(%d) = S(%d)−S(%d) = %+ld¢  (%+.2f%%)\n",
            t, t+1, t, (long)ret.eval(t), (double)ret.eval(t) / 100.0);

    std::printf("\n--- Discrete Gamma (Curvature) ---\n");
    std::printf("  Γ = Δ²(S) = Δ(Δ(S)): change in daily returns\n");
    auto gamma = forward_difference(ret);
    print_poly("Γ = Δ²S", gamma);
    std::printf("  Γ(0) = %ld¢ = constant second difference\n\n", (long)gamma.eval(0));

    std::printf("--- Cumulative P&L via Indefinite Sum ---\n");
    std::printf("  Σ(R)(t) = cumulative change = S(t) − S(0)\n");
    auto cum = indefinite_sum(ret);
    print_poly("Σ(ΔS)", cum);
    bool ok = true;
    for (int t = 1; t < 4; ++t) {
        uint64_t expect = S.eval(t) - S.eval(0);
        uint64_t actual = cum.eval(t);
        bool match = (expect == actual);
        std::printf("    Σ(R)(%d) = %lu  S(%d)−S(0) = %lu  %s\n",
            t, (unsigned long)actual, t, (unsigned long)expect,
            match ? "✓" : "✗");
        if (!match) ok = false;
    }
    std::printf("  Σ = Δ⁻¹ holds: %s\n\n", ok ? "yes" : "no");

    std::printf("--- Delta Hedge: Continuous vs Discrete ---\n");
    auto hedge = formal_derivative(S);
    auto d2 = formal_derivative(hedge);
    std::printf("  Operator identity Δ = e^D − I ≈ D + D²/2! + ...\n");
    std::printf("  D(S) ≈ Δ(S) − ½Δ²(S)  (invert for continuous delta):\n");
    for (int t = 0; t < 3; ++t) {
        uint64_t delta_exact  = forward_difference(S).eval(t);
        uint64_t deriv        = hedge.eval(t);
        uint64_t correction   = d2.eval(t) / 2;
        uint64_t delta_from_d = deriv + correction;
        std::printf("    t=%d: D(S)=%4lu  "
            "Δ(S)=%4lu  D+½D²=%4lu  match=%s\n",
            t, (unsigned long)deriv, (unsigned long)delta_exact,
            (unsigned long)delta_from_d,
            delta_exact == delta_from_d ? "yes" : "no");
    }

    std::printf("\n--- Key Insight ---\n");
    std::printf("  Δ = discrete derivative (daily returns)\n");
    std::printf("  Σ = discrete integral  (cumulative P&L)\n");
    std::printf("  Δ = e^D − I connects discrete → continuous calculus\n");
    std::printf("  All operators exact in ℤ/2^Wℤ — no rounding\n");

    return 0;
}
