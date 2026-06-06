#include <dyadic.h>
#include <pade.h>
#include <cstdint>
#include <cstdio>

using namespace dyadic;

int main() {
    std::printf("=== 15 — Padé Approximants via pade.h ===\n\n");

    using W = uint32_t;

    // [1/1] Padé for geometric series 1/(1-t) = 1 + t + t² + ...
    // P(t) / Q(t) should match first 2·1+1 = 3 series terms
    Polynomial<3, W, MonomialBasis> geom{{1, 1, 1}};
    auto [p, q] = pade_approximant<1, 1>(geom);

    std::printf("Series: 1 + t + t²\n");
    std::printf("[1/1] P(t) = %u + %u t\n", p[0], p[1]);
    std::printf("      Q(t) = %u + %u t\n", q[0], q[1]);
    std::printf("      (expected P=1, Q=1-t)\n\n");

    // Verify: P - Q·series = 0 up to t²
    bool ok = (p[0] == 1 && p[1] == 0 && q[1] == W(-1));
    std::printf("  %s\n\n", ok ? "PASS" : "FAIL");

    // [M/0] is just the truncated series (needs M+1 = 5 terms)
    Polynomial<5, W, MonomialBasis> geom5{{1, 1, 1, 1, 1}};
    auto [pm0, qm0] = pade_approximant<4, 0>(geom5);
    std::printf("[4/0] P(t) = %u + %u t + %u t² + %u t³ + %u t⁴\n",
            pm0[0], pm0[1], pm0[2], pm0[3], pm0[4]);
    std::printf("      Q(t) = %u\n", qm0[0]);

    return ok ? 0 : 1;
}
