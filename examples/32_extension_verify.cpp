// 32_extension_verify — Conformance tests for dyadic-example extension headers
// Returns 0 on pass, 1 on failure. Run via: cmake --build build --target run

#include <dyadic.h>
#include <dynamic_polynomial.h>
#include <pade.h>
#include <continued_fractions.h>
#include <matrix.h>
#include <cstdio>

using namespace dyadic;

static int failures = 0;
#define CHECK(cond, msg) do { \
    if (!(cond)) { std::printf("  FAIL  %s\n", msg); failures++; } \
    else { std::printf("  PASS  %s\n", msg); } \
} while(0)

// ---------------------------------------------------------------------------
// 1. DynamicPolynomial ring semantics: operator* must use carry-chain poly_mul
// ---------------------------------------------------------------------------
static void test_dynamic_ring_semantics() {
    std::printf("\n--- DynamicPolynomial ring semantics ---\n");

    // uint8_t: coefficients that overflow word boundary to test carry chain
    DynamicPolynomial<uint8_t, MonomialBasis> a{{200, 1}};  // 200 + x
    DynamicPolynomial<uint8_t, MonomialBasis> b{{2}};       // 2
    auto p = a * b;
    // Carry-chain: (200*2) = 400 → 400 mod 256 = 144, carry 1 → 1*2+1 = 3
    CHECK(p[0] == 144 && p[1] == 3, "carry-chain propagation in uint8_t");

    // Ring consistency: DynamicPolynomial × DynamicPolynomial == to_dynamic(static × static)
    Polynomial<3, uint32_t> sp{{1, 2, 3}};
    Polynomial<3, uint32_t> sq{{4, 5, 6}};
    auto s_prod = sp * sq;
    auto d_prod = to_dynamic(sp) * to_dynamic(sq);
    auto back = to_static<5>(d_prod);
    bool ok = true;
    for (int i = 0; i < 5; ++i) ok = ok && (s_prod[i] == back[i]);
    CHECK(ok, "ring consistency: to_static(to_dynamic(P) × to_dynamic(Q)) == P × Q");
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

    // Taylor roundtrip with small coefficients (precision window safe)
    DynamicPolynomial<uint32_t, MonomialBasis> small({1, 2, 3});
    auto taylor = change_basis<TaylorBasis>(small);
    auto tback = change_basis<MonomialBasis>(taylor);
    ok = true;
    for (int i = 0; i < small.size(); ++i) ok = ok && (small[i] == tback[i]);
    CHECK(ok, "Monomial → Taylor → Monomial roundtrip (small coefficients)");
}

// ---------------------------------------------------------------------------
// 3. DynamicPolynomial formal derivative / forward difference
// ---------------------------------------------------------------------------
static void test_dynamic_calculus() {
    std::printf("\n--- DynamicPolynomial calculus ---\n");

    DynamicPolynomial<uint32_t, MonomialBasis> p({1, 2, 3});
    auto dp = formal_derivative(p);
    // D(1 + 2x + 3x^2) = 2 + 6x
    CHECK(dp.size() == 2 && dp[0] == 2 && dp[1] == 6, "formal_derivative D");

    auto delta = forward_difference(p);
    // Δ(1 + 2x + 3x^2) = P(x+1)-P(x) = (2+3) + (2+6)x + 3x^2 ignoring truncation
    // Closed form: ΔP = 2*C(1,0) + 3*C(2,0) + 3*C(2,1)*x = 5 + 6x
    CHECK(delta.size() == 2 && delta[0] == 5 && delta[1] == 6, "forward_difference Δ");

    // Verify D and Δ commute: D(Δ(P)) = Δ(D(P))
    auto d_delta = formal_derivative(forward_difference(p));
    auto delta_d = forward_difference(formal_derivative(p));
    bool ok = true;
    for (int i = 0; i < d_delta.size(); ++i) ok = ok && (d_delta[i] == delta_d[i]);
    CHECK(ok, "D ∘ Δ = Δ ∘ D");
}

// ---------------------------------------------------------------------------
// 4. Padé approximant: [1/1] of geometric series recovers 1/(1-t)
// ---------------------------------------------------------------------------
static void test_pade() {
    std::printf("\n--- Padé approximants ---\n");

    // Geometric series 1 + t + t^2 → [1/1] = 1/(1-t)
    // Requires M+N+1 = 3 terms
    Polynomial<3, uint32_t, MonomialBasis> geom{{1, 1, 1}};
    auto [P, Q] = pade_approximant<1, 1>(geom);
    // P should be 1 (degree 0), Q should be (1, -1) = {1, W-1}
    CHECK(P[0] == 1 && P.actual_degree() == 0, "Padé [1/1] P = 1");
    CHECK(Q[0] == 1 && Q[1] == uint32_t(-1), "Padé [1/1] Q = 1 - t");
}

// ---------------------------------------------------------------------------
// 5. Continued fraction: rational function produces finite CF
// ---------------------------------------------------------------------------
static void test_continued_fractions() {
    std::printf("\n--- Continued fractions ---\n");

    // For a geometric series 1 + t + t^2 + ..., the CF coefficients are all 1
    DynamicPolynomial<uint32_t, MonomialBasis> geom({1, 1, 1, 1, 1, 1});
    auto c = cf_expand(geom, 5);
    bool ok = true;
    for (size_t i = 0; i < c.size(); ++i) ok = ok && (c[i] == 1);
    CHECK(ok && c.size() >= 4, "geometric series → CF coefficients all 1");
}

// ---------------------------------------------------------------------------
// 6. Matrix: det(AB) = det(A)det(B) and inverse roundtrip
// ---------------------------------------------------------------------------
static void test_matrix() {
    std::printf("\n--- Matrix ---\n");

    // 2×2 matrices with odd diagonal for invertibility
    Matrix<2, 2, uint32_t> A{{ {{3, 1}, {2, 5}} }};
    Matrix<2, 2, uint32_t> B{{ {{7, 2}, {3, 9}} }};

    // det(AB) = det(A)det(B)
    auto detA = A.determinant();
    auto detB = B.determinant();
    auto detAB = (A * B).determinant();
    CHECK(detAB == detA * detB, "det(AB) = det(A) · det(B)");

    // A * A^{-1} = I
    auto inv = A.inverse();
    auto I = A * inv;
    bool ok = true;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            ok = ok && (I[i][j] == (i == j ? uint32_t(1) : uint32_t(0)));
    CHECK(ok, "A · A^{-1} = I");

    // Solve: A · x = b → A · x should recover b
    std::array<uint32_t, 2> b{{10, 19}};
    auto x = A.solve(b);
    uint32_t check0 = A[0][0] * x[0] + A[0][1] * x[1];
    uint32_t check1 = A[1][0] * x[0] + A[1][1] * x[1];
    ok = (check0 == b[0] && check1 == b[1]);
    CHECK(ok, "A · solve(A, b) = b");
}

int main() {
    std::printf("=== 32 — Extension Conformance Tests ===\n");

    test_dynamic_ring_semantics();
    test_dynamic_basis_roundtrip();
    test_dynamic_calculus();
    test_pade();
    test_continued_fractions();
    test_matrix();

    std::printf("\n=== %s ===\n", failures ? "SOME TESTS FAILED" : "ALL TESTS PASSED");
    return failures ? 1 : 0;
}
