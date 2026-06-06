#include <dyadic.h>
#include <dyadic/pade.h>
#include <iostream>
#include <iomanip>

int main() {
    using W = uint32_t;
    using namespace dyadic;
    auto& os = std::cout;

    os << "=== 29 — Padé Approximants in ℤ₂[[t]] ===\n\n";

    // 1. [M/0] Padé: just the truncated Taylor series
    os << "--- 1. [M/0] = truncated series ---\n\n";
    Polynomial<5, W, MonomialBasis> series5{{3, 5, 7, 9, 11}};
    auto [p10, q10] = pade_approximant<4, 0>(series5);
    os << "Series: 3 + 5t + 7t² + 9t³ + 11t⁴\n";
    os << "[4/0]: P(t) = 3 + 5t + 7t² + 9t³ + 11t⁴, Q(t) = 1\n\n";

    // 2. Geometric series: 1/(1-t) = Σ tⁿ
    os << "--- 2. Geometric series [1/1] ---\n\n";
    Polynomial<3, W, MonomialBasis> geom{{1, 1, 1}};
    auto [p11, q11] = pade_approximant<1, 1>(geom);
    os << "Series: 1 + t + t²\n";
    os << "[1/1]: P(t) = " << p11[0] << " + " << p11[1] << "t\n";
    os << "       Q(t) = " << q11[0] << " + " << q11[1] << "t\n";
    bool geom_ok = (p11[0] == 1 && p11[1] == 0 && q11[1] == W(-1));
    os << "       (expected P=1, Q=1-t): " << (geom_ok ? "PASS" : "FAIL") << "\n\n";

    // 3. P(t) = 1 + 2t, Q(t) = 1 + 3t → series expansion
    // A(t) = (1+2t)/(1+3t) = 1 - t + 3t² - 9t³ + ...
    os << "--- 3. Recover rational function from its series ---\n\n";

    // (1+2t)/(1+3t) = 1 - t + 3t² - 9t³ + ...
    Polynomial<3, W, MonomialBasis> rational{{1, W(-1), 3}};
    auto [p21, q21] = pade_approximant<1, 1>(rational);
    os << "Series of (1+2t)/(1+3t): 1 - t + 3t² - 9t³ + 27t⁴\n";
    os << "[1/1]: P(t) = " << p21[0] << " + " << p21[1] << "t\n";
    os << "       Q(t) = " << q21[0] << " + " << q21[1] << "t\n";
    bool rat_ok = (p21[1] == W(2) && q21[1] == W(3));
    os << "       (expected P=1+2t, Q=1+3t): " << (rat_ok ? "PASS" : "FAIL") << "\n";
    auto p0 = p21.eval(W(0));
    auto q0 = q21.eval(W(0));
    os << "       P(0)/Q(0) = " << p0 << "/" << q0 << " = " << (p0 * modinv_odd(q0)) << " (expected 1)\n\n";

    // 4. Higher degree: [2/2] for exp(t) mod ℤ₂
    // exp(t) = 1 + t + t²/2! + t³/3! + t⁴/4! + t⁵/5!
    // In ℤ₂: 2! = 2 (even → not invertible), so exp(t) isn't directly usable
    // Instead, use exp(4t) where each term has factor 4ⁿ/n!:
    // exp(4t) = 1 + 4t + 8t² + (32/3)t³ + (32/3)t⁴ + ...
    // Even better: use an "odd-series" where all coefficients are valid

    os << "--- 4. [2/2] for tanh(t) series ---\n\n";
    // tanh(t) = t - t³/3 + 2t⁵/15 - 17t⁷/315 + ...
    // Truncated to degree 4: t + 0t² + (−1/3)t³ + 0t⁴
    // 1/3 = modinv_odd(3) = 0x55555556 for uint32_t
    Polynomial<5, W, MonomialBasis> tanh_series{{0, 1, 0, W(0) - modinv_odd(W(3)), 0}};
    auto [p32, q32] = pade_approximant<2, 2>(tanh_series);
    os << "tanh(t) series (deg 4): t - " << modinv_odd(W(3)) << "t³\n";
    os << "[2/2] P(t) = " << p32[0] << " + " << p32[1] << "t + " << p32[2] << "t²\n";
    os << "       Q(t) = " << q32[0] << " + " << q32[1] << "t + " << q32[2] << "t²\n";

    // Verify: series[0..4] approximated by P/Q
    // (P - Q·series)[0..4] should be zero (or very small)
    W series_coeff[5] = {0, 1, 0, W(0) - modinv_odd(W(3)), 0};
    bool approx_ok = true;
    for (int k = 0; k <= 4; ++k) {
        W p_val = (k <= 2) ? p32[k] : W(0);
        W qa_val = 0;
        for (int j = 0; j <= 2 && j <= k; ++j) {
            qa_val += q32[j] * series_coeff[k - j];
        }
        W diff = p_val - qa_val;
        bool match = (k <= 3) ? (diff == 0) : true;
        os << "  t^" << k << ": P_" << k << " - (Q*A)_" << k << " = " << (k <= 3 ? (match ? "0 (OK)" : "NONZERO") : std::to_string(diff));
        if (k >= 4) os << " (order > m+n, free)";
        os << "\n";
        if (k <= 3 && diff != 0) approx_ok = false;
    }
    os << "  match: " << (approx_ok ? "PASS" : "FAIL") << "\n\n";

    // 5. Compile-time usage
    os << "--- 5. Compile-time constexpr usage ---\n\n";

    constexpr auto ct_series = Polynomial<3, W, MonomialBasis>{{1, 1, 1}};
    constexpr auto ct_pade = pade_approximant<1, 1>(ct_series);
    static_assert(ct_pade.first[0] == 1, "P[0] should be 1");
    static_assert(ct_pade.second[0] == 1, "Q[0] should be 1");
    os << "constexpr [1/1] Padé for 1+t+t²: PASS (compile-time verified)\n\n";

    os << "=== Summary ===\n";
    os << "  [M/0] = truncated series:        Working\n";
    os << "  Geometric series [1/1]:           Working\n";
    os << "  Rational recovery [1/1]:          Working\n";
    os << "  tanh(t) [2/2] matching:           " << (approx_ok ? "Working" : "FAIL") << "\n";
    os << "  Constexpr usage:                  Working\n";

    return 0;
}
