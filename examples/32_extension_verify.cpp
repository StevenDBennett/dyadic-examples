// 32_extension_verify — Conformance tests for dyadic-example extension headers
// Returns 0 on pass, 1 on failure. Run via: cmake --build build --target run

#include <dyadic.h>
#include <dynamic_polynomial.h>
#include <pade.h>
#include <continued_fractions.h>
#include <matrix.h>
#include <cstdio>

using namespace dyadic;

using M22 = Matrix<2, 2, uint32_t>;
using M32 = Matrix<3, 2, uint32_t>;
using M23 = Matrix<2, 3, uint32_t>;

static int failures = 0;
#define CHECK(cond, msg) do { \
    if (!(cond)) { std::printf("  FAIL  %s\n", msg); failures++; } \
    else { std::printf("  PASS  %s\n", msg); } \
} while(0)

// ---------------------------------------------------------------------------
// 1. DynamicPolynomial ring semantics
// ---------------------------------------------------------------------------
static void test_dynamic_ring_semantics() {
    std::printf("\n--- DynamicPolynomial ring semantics ---\n");

    DynamicPolynomial<uint8_t, MonomialBasis> a{{200, 1}};
    DynamicPolynomial<uint8_t, MonomialBasis> b{{2}};
    auto p = a * b;
    CHECK(p[0] == 144 && p[1] == 3, "carry-chain propagation in uint8_t");

    Polynomial<3, uint32_t> sp{{1, 2, 3}};
    Polynomial<3, uint32_t> sq{{4, 5, 6}};
    auto s_prod = sp * sq;
    auto d_prod = to_dynamic(sp) * to_dynamic(sq);
    auto back = to_static<5>(d_prod);
    bool ok = true;
    for (int i = 0; i < 5; ++i) ok = ok && (s_prod[i] == back[i]);
    CHECK(ok, "ring consistency: to_static(to_dynamic(P) × to_dynamic(Q)) == P × Q");
}

static void test_dynamic_empty_and_truncation() {
    std::printf("\n--- DynamicPolynomial edge cases ---\n");

    DynamicPolynomial<uint32_t> empty;
    CHECK(empty.degree() == -1, "empty polynomial degree == -1");
    CHECK(empty.eval(42) == 0, "empty polynomial eval == 0");
    CHECK(empty.size() == 0, "empty polynomial size == 0");

    auto prod = empty * DynamicPolynomial<uint32_t>({1, 2, 3});
    CHECK(prod.size() == 0, "empty * nonempty == empty");

    auto sum = empty + DynamicPolynomial<uint32_t>({5});
    CHECK(sum.size() == 1 && sum[0] == 5, "empty + scalar == scalar");

    // to_static<N> truncation when dynamic poly exceeds static size
    DynamicPolynomial<uint32_t> big({1, 2, 3, 4, 5});
    auto truncated = to_static<3>(big);
    bool ok = (truncated.size() == 3 && truncated[0] == 1 && truncated[1] == 2 && truncated[2] == 3);
    CHECK(ok, "to_static<N> truncates excess coefficients");

    // Mixed-degree multiplication
    DynamicPolynomial<uint32_t> p2({1, 2});
    DynamicPolynomial<uint32_t> p3({3, 4, 5});
    auto m = p2 * p3;
    // (1 + 2x)(3 + 4x + 5x^2) = 3 + 10x + 13x^2 + 10x^3
    CHECK(m.size() == 4 && m[0] == 3 && m[1] == 10 && m[2] == 13 && m[3] == 10,
          "mixed-degree carry-chain multiplication");

    // Mixed-degree addition with auto-resize
    auto s = p2 + p3;
    CHECK(s.size() == 3 && s[0] == 4 && s[1] == 6 && s[2] == 5,
          "mixed-degree addition auto-resize");

    // Degree-0 polynomial
    DynamicPolynomial<uint32_t> deg0({42});
    CHECK(deg0.degree() == 0, "degree-0 polynomial");
    CHECK(deg0.eval(100) == 42, "degree-0 eval constant");
}

// ---------------------------------------------------------------------------
// 2. DynamicPolynomial basis roundtrip
// ---------------------------------------------------------------------------
static void test_dynamic_basis_roundtrip() {
    std::printf("\n--- DynamicPolynomial basis roundtrip ---\n");

    DynamicPolynomial<uint32_t, MonomialBasis> mono({1, 2, 3, 4});
    auto ff = change_basis<FallingFactorialBasis>(mono);
    auto back = change_basis<MonomialBasis>(ff);
    bool ok = true;
    for (int i = 0; i < mono.size(); ++i) ok = ok && (mono[i] == back[i]);
    CHECK(ok, "Monomial → FallingFactorial → Monomial roundtrip");

    DynamicPolynomial<uint32_t, MonomialBasis> small({1, 2, 3});
    auto taylor = change_basis<TaylorBasis>(small);
    auto tback = change_basis<MonomialBasis>(taylor);
    ok = true;
    for (int i = 0; i < small.size(); ++i) ok = ok && (small[i] == tback[i]);
    CHECK(ok, "Monomial → Taylor → Monomial roundtrip (small coefficients)");
}

// ---------------------------------------------------------------------------
// 3. DynamicPolynomial calculus
// ---------------------------------------------------------------------------
static void test_dynamic_calculus() {
    std::printf("\n--- DynamicPolynomial calculus ---\n");

    DynamicPolynomial<uint32_t, MonomialBasis> p({1, 2, 3});
    auto dp = formal_derivative(p);
    CHECK(dp.size() == 2 && dp[0] == 2 && dp[1] == 6, "formal_derivative D");

    auto delta = forward_difference(p);
    CHECK(delta.size() == 2 && delta[0] == 5 && delta[1] == 6, "forward_difference Δ");

    auto d_delta = formal_derivative(forward_difference(p));
    auto delta_d = forward_difference(formal_derivative(p));
    bool ok = true;
    for (int i = 0; i < d_delta.size(); ++i) ok = ok && (d_delta[i] == delta_d[i]);
    CHECK(ok, "D ∘ Δ = Δ ∘ D");
}

// ---------------------------------------------------------------------------
// 4. Padé approximant
// ---------------------------------------------------------------------------
static void test_pade() {
    std::printf("\n--- Padé approximants ---\n");

    Polynomial<3, uint32_t, MonomialBasis> geom{{1, 1, 1}};
    auto [P, Q] = pade_approximant<1, 1>(geom);
    CHECK(P[0] == 1 && P.actual_degree() == 0, "Padé [1/1] P = 1");
    CHECK(Q[0] == 1 && Q[1] == uint32_t(-1), "Padé [1/1] Q = 1 - t");

    // [0/0] Padé: just the constant term
    Polynomial<1, uint32_t> const_series{{7}};
    auto [P00, Q00] = pade_approximant<0, 0>(const_series);
    CHECK(P00[0] == 7 && Q00[0] == 1, "Padé [0/0] P_0 = a_0, Q_0 = 1");

    // [M/0] Padé: truncated series
    Polynomial<5, uint32_t> series5{{1, 2, 3, 4, 5}};
    auto [PM0, QM0] = pade_approximant<4, 0>(series5);
    bool ok = true;
    for (int i = 0; i < 5; ++i) ok = ok && (PM0[i] == series5[i]);
    CHECK(ok && QM0[0] == 1, "Padé [4/0] P = series, Q = 1");

    // Even a_M triggers silent failure (zero return)
    // a_1 = 2 (even) — a_M odd requirement violated
    Polynomial<4, uint32_t> even_series{{1, 2, 4, 8}};
    auto [Pf, Qf] = pade_approximant<2, 1>(even_series);
    CHECK(Pf[0] == 0 && Pf.actual_degree() == -1, "Padé with even a_M returns zero P");
    CHECK(Qf[0] == 0 && Qf.actual_degree() == -1, "Padé with even a_M returns zero Q");

    // [2/2] for a known series
    Polynomial<5, uint32_t> tanh_like{{0, 1, 0, uint32_t(-2863311531u), 0}};
    auto [P22, Q22] = pade_approximant<2, 2>(tanh_like);
    CHECK(P22[0] == 0 && P22[1] == 1 && P22[2] == 0, "Padé [2/2] P = t");
    CHECK(Q22[0] == 1 && Q22[1] == 0, "Padé [2/2] Q_0 = 1, Q_1 = 0");
}

// ---------------------------------------------------------------------------
// 5. Continued fractions
// ---------------------------------------------------------------------------
static void test_continued_fractions() {
    std::printf("\n--- Continued fractions ---\n");

    DynamicPolynomial<uint32_t, MonomialBasis> geom({1, 1, 1, 1, 1, 1});
    auto c = cf_expand(geom, 5);
    bool ok = true;
    for (size_t i = 0; i < c.size(); ++i) ok = ok && (c[i] == 1);
    CHECK(ok && c.size() >= 4, "geometric series → CF coefficients all 1");

    // Single-term series: CF contains just the constant
    auto c1 = cf_expand(DynamicPolynomial<uint32_t>({42}), 10);
    CHECK(c1.size() == 1 && c1[0] == 42, "single-term series → CF [42]");

    // All-zero series: all CF coefficients are 0
    auto cz = cf_expand(DynamicPolynomial<uint32_t>({0, 0, 0, 0}), 4);
    ok = true;
    for (auto v : cz) ok = ok && (v == 0);
    CHECK(ok && cz.size() == 4, "all-zero series → CF coefficients all 0");

    // cf_convergent n=0: returns c_0 / 1
    std::vector<uint32_t> cf = {3, 5, 7};
    auto [P0, Q0] = cf_convergent(cf, 0);
    CHECK(P0.eval(0) == 3 && Q0.eval(0) == 1, "cf_convergent n=0: P = c_0, Q = 1");

    // cf_convergent n=1: P_1 = c_0*c_1 + t, Q_1 = c_1
    auto [P1, Q1] = cf_convergent(cf, 1);
    // P_1 = c_1*P_0 + t*P_{-1} = 5*3 + t*1 = 15 + t
    // Q_1 = c_1*Q_0 + t*Q_{-1} = 5*1 + t*0 = 5
    CHECK(P1.eval(2) == 15 + 2 && Q1.eval(2) == 5, "cf_convergent n=1 recurrence");

    // Polynomial CF terminates within its degree+1 coefficients
    // 1 + 2t + 3t^2 has CF [1, 2, 3] (3 coefficients for degree 2)
    DynamicPolynomial<uint32_t> poly_series({1, 2, 3});
    auto cp = cf_expand(poly_series, 10);
    CHECK(cp.size() <= 3, "polynomial CF terminates within degree+1 terms");

    // cf_eval match for simple series
    auto val = cf_eval(cf, 1, uint32_t(0));
    CHECK(val.size() == 1 && val[0] == 3, "cf_eval at t=0 equals c_0");
}

// ---------------------------------------------------------------------------
// 6. Matrix operations
// ---------------------------------------------------------------------------
static void test_matrix() {
    std::printf("\n--- Matrix ---\n");

    M22 A{{ {{3, 1}, {2, 5}} }};
    M22 B{{ {{7, 2}, {3, 9}} }};

    auto detA = A.determinant();
    auto detB = B.determinant();
    auto detAB = (A * B).determinant();
    CHECK(detAB == detA * detB, "det(AB) = det(A) · det(B)");

    auto inv = A.inverse();
    auto I = A * inv;
    CHECK(I == M22::identity(), "A · A^{-1} = I");

    std::array<uint32_t, 2> b{{10, 19}};
    auto x = A.solve(b);
    uint32_t check0 = A[0][0] * x[0] + A[0][1] * x[1];
    uint32_t check1 = A[1][0] * x[0] + A[1][1] * x[1];
    bool ok = (check0 == b[0] && check1 == b[1]);
    CHECK(ok, "A · solve(A, b) = b");

    // trace
    CHECK(A.trace() == 3 + 5, "Matrix trace = sum of diagonal");

    // singular matrix detection
    M22 singular{{ {{2, 4}, {6, 8}} }};
    CHECK(singular.is_singular(), "is_singular on even-determinant matrix");
    CHECK(singular.determinant() % 2 == 0, "singular matrix has even determinant");

    // singular solve returns zero
    auto xs = singular.solve(b);
    CHECK(xs[0] == 0 && xs[1] == 0, "singular solve returns zero vector");

    // singular inverse returns zero matrix
    auto sinv = singular.inverse();
    CHECK(sinv == M22::zero(), "singular inverse returns zero matrix");

    // rank on rectangular matrices
    M32 tall{{ {{1, 2}, {3, 4}, {5, 6}} }};
    CHECK(tall.rank() >= 1, "tall (M>N) matrix rank > 0");

    M23 wide{{ {{1, 2, 3}, {4, 5, 6}} }};
    CHECK(wide.rank() >= 1, "wide (M<N) matrix rank > 0");

    // operator== and operator!=
    M22 C = A;
    CHECK(A == C, "Matrix operator== on equal matrices");
    CHECK(A != B, "Matrix operator!= on different matrices");
    CHECK(A + B == B + A, "Matrix addition commutativity via operator==");
}

int main() {
    std::printf("=== 32 — Extension Conformance Tests ===\n");

    test_dynamic_ring_semantics();
    test_dynamic_empty_and_truncation();
    test_dynamic_basis_roundtrip();
    test_dynamic_calculus();
    test_pade();
    test_continued_fractions();
    test_matrix();

    std::printf("\n=== %s ===\n", failures ? "SOME TESTS FAILED" : "ALL TESTS PASSED");
    return failures ? 1 : 0;
}
