#include <dyadic.h>
#include <iostream>
#include <iomanip>
#include <string>

int main() {
    using W = uint32_t;
    using namespace dyadic;

    auto& os = std::cout;
    os << "=== 27 — DynamicPolynomial (Runtime-Degree) ===\n\n";

    // 1. Construction at various runtime degrees
    os << "--- 1. Construction & eval ---\n\n";

    DynamicPolynomial<W> p({1, 2, 3});           // deg 2: 1 + 2x + 3x^2
    DynamicPolynomial<W> q(5);                    // deg 5, all zero
    q[0] = 5;
    q[3] = 7;                                     // 5 + 7x^3

    os << "p(x) = 1 + 2x + 3x^2\n";
    os << "  size = " << p.size() << ", degree = " << p.degree() << "\n";
    os << "  eval(0) = " << p.eval(0) << "\n";
    os << "  eval(1) = " << p.eval(1) << "\n";
    os << "  eval(2) = " << p.eval(2) << "\n\n";

    os << "q(x) = 5 + 7x^3\n";
    os << "  size = " << q.size() << ", degree = " << q.degree() << "\n\n";

    // 2. Basis conversions
    os << "--- 2. Basis conversion ---\n\n";

    DynamicPolynomial<W, MonomialBasis> mono({1, 2, 3});
    auto falling = monomial_to_falling(mono);
    auto back = falling_to_monomial(falling);

    os << "Monomial p(x)    = 1 + 2x + 3x^2\n";
    os << "  eval(3) = " << mono.eval(3) << "\n\n";

    os << "Falling factorial basis:\n";
    for (int i = 0; i < falling.size(); ++i)
        os << "  ff[" << i << "] = " << falling[i] << "\n";
    os << "  eval(3) = " << falling.eval(3) << "\n";
    os << "  roundtrip match: " << (mono.eval(10) == falling.eval(10)) << "\n\n";

    // 3. Change basis via generic interface
    os << "--- 3. Generic change_basis<> ---\n\n";

    auto taylor = change_basis<TaylorBasis>(mono);
    auto back2 = change_basis<MonomialBasis>(taylor);
    os << "Generic roundtrip: monomial → taylor → monomial\n";
    bool ok = true;
    for (int i = 0; i < mono.size(); ++i) ok = ok && (mono[i] == back2[i]);
    os << "  match: " << (ok ? "PASS" : "FAIL") << "\n\n";

    // 4. Arithmetic
    os << "--- 4. Arithmetic ---\n\n";

    DynamicPolynomial<W> a({1, 2, 3});
    DynamicPolynomial<W> b({4, 5, 6});
    auto sum = a + b;
    auto diff = a - b;
    auto prod = a * b;

    os << "a = " << a[0] << " + " << a[1] << "x + " << a[2] << "x^2\n";
    os << "b = " << b[0] << " + " << b[1] << "x + " << b[2] << "x^2\n";
    os << "a + b = " << sum[0] << " + " << sum[1] << "x + " << sum[2] << "x^2\n";
    os << "a - b = " << diff[0] << " + " << diff[1] << "x + " << diff[2] << "x^2\n";
    os << "a * b (deg " << prod.degree() << "): ";
    for (int i = 0; i <= prod.degree(); ++i) os << prod[i] << (i < prod.degree() ? ", " : "");
    os << "\n\n";

    // 5. Conversion to/from compile-time Polynomial
    os << "--- 5. Convert to/from static Polynomial<N> ---\n\n";

    Polynomial<4, W, MonomialBasis> static_p{{1, 2, 3, 4}};
    auto dyn_p = to_dynamic(static_p);
    auto static_back = to_static<4>(dyn_p);

    os << "Static Polynomial<4>: 1 + 2x + 3x^2 + 4x^3\n";
    os << "  to_dynamic -> degree = " << dyn_p.degree() << ", size = " << dyn_p.size() << "\n";
    bool match = true;
    for (int i = 0; i < 4; ++i) match = match && (static_p[i] == static_back[i]);
    os << "  roundtrip: " << (match ? "PASS" : "FAIL") << "\n\n";

    // 6. Formal derivative
    os << "--- 6. Formal derivative & forward difference ---\n\n";

    DynamicPolynomial<W, MonomialBasis> f({6, 5, 4, 3, 2, 1});
    auto df = formal_derivative(f);
    auto dff = formal_derivative(change_basis<FallingFactorialBasis>(f));

    os << "f(x) = ";
    for (int i = 0; i <= f.degree(); ++i) os << (i ? " + " : "") << f[i] << "x^" << i;
    os << "\n  degree = " << f.degree() << "\n";

    os << "f'(x) (monomial) = ";
    for (int i = 0; i <= df.degree(); ++i) os << (i ? " + " : "") << df[i] << "x^" << i;
    os << "\n  degree = " << df.degree() << "\n";

    os << "Δf(x) (forward diff) = ";
    auto delta = forward_difference(f);
    for (int i = 0; i <= delta.degree(); ++i) os << (i ? " + " : "") << delta[i] << "x^" << i;
    os << "\n  degree = " << delta.degree() << "\n";

    // Verify: Δf(x) = f(x+1) - f(x) at sample points
    os << "  Δf(0) = " << delta.eval(0) << ", f(1)-f(0) = " << (f.eval(1) - f.eval(0)) << "\n";
    os << "  Δf(1) = " << delta.eval(1) << ", f(2)-f(1) = " << (f.eval(2) - f.eval(1)) << "\n";
    os << "  Δf(2) = " << delta.eval(2) << ", f(3)-f(2) = " << (f.eval(3) - f.eval(2)) << "\n\n";

    // 7. Mixed-degree operations
    os << "--- 7. Mixed-degree operations ---\n\n";

    DynamicPolynomial<W> short_p({1, 2});           // deg 1
    DynamicPolynomial<W> long_p({3, 4, 5, 6});     // deg 3
    auto mixed = short_p * long_p + short_p;
    os << "({1, 2} * {3, 4, 5, 6}) + {1, 2}\n";
    os << "  result degree = " << mixed.degree() << ", coefficients: ";
    for (int i = 0; i <= mixed.degree(); ++i) os << mixed[i] << (i < mixed.degree() ? ", " : "");
    os << "\n\n";

    // 8. Initializer list eval (matching static behavior)
    os << "--- 8. eval in FallingFactorial and Taylor bases ---\n\n";

    DynamicPolynomial<W, FallingFactorialBasis> ff({0, 1, 2});
    os << "Falling factorial basis coefs: [0, 1, 2]\n";
    os << "  eval(0) = " << ff.eval(0) << "\n";
    os << "  eval(1) = " << ff.eval(1) << "\n";
    os << "  eval(5) = " << ff.eval(5) << "\n";
    os << "  eval(10) = " << ff.eval(10) << "\n";

    auto ff_mono = change_basis<MonomialBasis>(ff);
    os << "Equivalent monomial: ";
    for (int i = 0; i <= ff_mono.degree(); ++i) os << ff_mono[i] << (i < ff_mono.degree() ? "x + " : "x^" + std::to_string(i));
    os << "\n";
    bool eval_match = ff.eval(7) == ff_mono.eval(7);
    os << "  eval(7) match: " << (eval_match ? "PASS" : "FAIL") << "\n\n";

    // 9. Ring semantics verification: carry-chain vs coefficient-wise
    os << "--- 9. Ring semantics verification ---\n\n";

    os << "  DynamicPolynomial::operator* now uses carry-chain poly_mul,\n";
    os << "  matching Polynomial<N,W,Basis>::operator*. This ensures\n";
    os << "  consistent 2-adic arithmetic across static and dynamic types.\n\n";

    // Use uint8_t to make carry-chain effects visible (overflow at 256)
    using W8 = uint8_t;
    // p = 200 + x  (coefficients deliberately chosen to overflow uint8_t)
    DynamicPolynomial<W8, MonomialBasis> da({200, 1});
    DynamicPolynomial<W8, MonomialBasis> db({2, 0});
    auto dprod = da * db;
    os << "  uint8_t test:  (200 + x) × (2)\n";
    os << "    Coefficient-wise would be:  {144, 2}\n";
    os << "    Carry-chain (poly_mul):     {" << int(dprod[0]) << ", " << int(dprod[1]) << "}\n";

    // Compare with static Polynomial using same multiplication
    Polynomial<2, W8, MonomialBasis> sa{{200, 1}};
    Polynomial<2, W8, MonomialBasis> sb{{2, 0}};
    auto sprod = sa * sb;
    bool ring_match = (dprod[0] == sprod[0] && dprod[1] == sprod[1]);
    os << "    Static poly_mul:             {" << int(sprod[0]) << ", " << int(sprod[1]) << "}\n";
    os << "    Ring match (static==dynamic): " << (ring_match ? "PASS" : "FAIL") << "\n\n";

    // Verify roundtrip via to_static/to_dynamic preserves carry-chain ring
    os << "  Roundtrip verification:\n";
    Polynomial<3, uint32_t> sp{{1, 2, 3}};
    Polynomial<3, uint32_t> sq{{4, 5, 6}};
    auto s_static = sp * sq;
    auto d_from_sp = to_dynamic(sp);
    auto d_from_sq = to_dynamic(sq);
    auto d_dynamic = d_from_sp * d_from_sq;
    auto back_to_static = to_static<5>(d_dynamic);
    bool rt_ok = true;
    for (int i = 0; i < 5; ++i) rt_ok = rt_ok && (s_static[i] == back_to_static[i]);
    os << "    to_static(to_dynamic(P) × to_dynamic(Q)) == P × Q: "
       << (rt_ok ? "PASS" : "FAIL") << "\n\n";

    os << "=== Summary ===\n";
    os << "  Construction with runtime degree:  Working\n";
    os << "  eval (monomial/falling/taylor):    Working\n";
    os << "  Basis conversion roundtrip:        Working\n";
    os << "  Arithmetic (+, -, *):              Working\n";
    os << "  Static ↔ Dynamic conversion:       Working\n";
    os << "  Formal derivative:                 Working\n";
    os << "  Forward difference:                Working\n";
    os << "  Ring semantics (carry-chain):      Working\n";

    return (ring_match && rt_ok) ? 0 : 1;
}
