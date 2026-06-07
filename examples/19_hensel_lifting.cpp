#include <dyadic.h>
#include <cstdio>

using namespace dyadic;

// Hensel's Lemma: If f(a) ≡ 0 (mod p) and f'(a) not≡ 0 (mod p),
// then there exists a unique root a* in ℤ_p with a* ≡ a (mod p).
// Lift: a_{n+1} = a_n - f(a_n) * inv(f'(a_n))  (mod p^{2n})

// Careful Hensel lift for p=2:
// Given a solution to f(x) ≡ 0 (mod 2^k), lift to (mod 2^{2k})
template<std::unsigned_integral W, typename Func, typename Deriv>
constexpr W hensel_lift(W a_k, Func f, Deriv f_prime, int bit_precision) {
    W a = a_k;
    for (int i = 0; i < bit_precision; ++i) {
        W f_val = f(a);
        W f_p = f_prime(a);
        if (f_val == 0) break;
        // If f'(a) is odd (unit in ℤ₂), we can invert it
        if (f_p % 2 == 1) {
            W inv = modinv_odd(f_p);
            a = a - f_val * inv;
        } else {
            // f'(a) even: need more careful handling (degenerate case)
            break;
        }
    }
    return a;
}

int main() {
    std::printf("=== 19 — Hensel Lifting ===\n\n");

    // Example: Lift root of f(x) = x^2 - 1 = 0 from mod 2 to mod 2^64
    // Roots in ℤ₂: x = 1 and x = -1 (as a power series)
    auto f = [](uint64_t x) -> uint64_t { return x * x - 1; };
    auto f_prime = [](uint64_t x) -> uint64_t { return 2 * x; };

    uint64_t root = hensel_lift<uint64_t>(uint64_t{1}, f, f_prime, 64);
    std::printf("Root of x^2 - 1 = 0:\n");
    std::printf("  Lifted from mod 2: x = %lu\n", root);
    std::printf("  Verify: %lu^2 - 1 = %lu\n", root, f(root));
    std::printf("  f(x) == 0: %s\n\n", (f(root) == 0) ? "PASS" : "FAIL");

    // Example: f(x) = x^3 - x = x(x-1)(x+1) = 0
    // Roots: 0, 1, -1 in ℤ₂
    auto g = [](uint64_t x) -> uint64_t { return x * x * x - x; };
    auto g_prime = [](uint64_t x) -> uint64_t { return 3 * x * x - 1; };

    for (uint64_t guess : {uint64_t{0}, uint64_t{1}, uint64_t{0} - 1}) {
        uint64_t lifted = hensel_lift(guess, g, g_prime, 64);
        std::printf("Root of x^3 - x = 0, starting from x₀ = %lu:\n", guess);
        std::printf("  Lifted: x = %lu\n", lifted);
        std::printf("  Verify: f(%lu) = %lu\n", lifted, g(lifted));
        std::printf("  Match: %s\n\n", (g(lifted) == 0) ? "PASS" : "FAIL");
    }

    // Example: f(x) = x^2 + x + 1 (roots are primitive cube roots of unity)
    // Check if it has roots in ℤ₂ (it doesn't, since discriminant = -3 ≡ 1 mod 8?)
    auto h = [](uint64_t x) -> uint64_t { return x * x + x + 1; };
    auto h_prime = [](uint64_t x) -> uint64_t { return 2 * x + 1; };

    uint64_t h_root = hensel_lift<uint64_t>(uint64_t{1}, h, h_prime, 64);
    std::printf("Root of x^2 + x + 1 = 0, starting from x₀ = 1:\n");
    std::printf("  Lifted: x = %lu\n", h_root);
    std::printf("  Verify: f(%lu) = %lu\n", h_root, h(h_root));
    std::printf("  Note: No root in ℤ₂ since this is irreducible over ℤ₂\n\n");

    return 0;
}
