#include <dyadic.h>
#include <cstdio>

using namespace dyadic;

// Newton-Raphson for sqrt in ℤ₂: x_{n+1} = x_n - (x_n² - a) / (2*x_n)
// In ℤ₂, division by 2 requires exact 2-adic division (diff must be divisible by 2)
// For the Newton method to converge: a ≡ 1 (mod 8) gives a unique square root in ℤ₂
template<std::unsigned_integral W>
constexpr W p_adic_sqrt(W a, int iterations = 12) {
    W x = 1;  // sqrt(1) ≡ 1 (mod 2); also valid when a ≡ 1 (mod 8)
    for (int i = 0; i < iterations; ++i) {
        W f = x * x - a;
        if (f == 0) break;
        // Newton step: x_{n+1} = x_n - f / (2*x_n)
        // Use 2-adic division for the denominator 2
        // f is divisible by increasing powers of 2 as we converge
        W inv_x = modinv_odd(x);  // x must be odd for invertibility
        W correction = f * inv_x;  // f/x
        x = x - div_2k_adic(correction, 1);  // correction/2
    }
    return x;
}

// Cube root: x_{n+1} = x_n - (x_n³ - a) / (3*x_n²)
// For odd x: f'(x) = 3x² is odd → fully invertible in ℤ₂
template<std::unsigned_integral W>
constexpr W p_adic_cuberoot(W a, int iterations = 12) {
    W x = 1;
    for (int i = 0; i < iterations; ++i) {
        W f = x * x * x - a;
        if (f == 0) break;
        W fp = 3 * x * x;  // Derivative = 3x²
        W inv_fp = modinv_odd(fp);  // fp is odd → invertible
        x = x - f * inv_fp;
    }
    return x;
}

// Show convergence by printing successive approximations
template<std::unsigned_integral W>
constexpr void show_convergence(W a, const char* label, W(*method)(W, int), int steps = 6) {
    std::printf("  %s for a=%lu:\n", label, (unsigned long)a);
    W x = 1;
    for (int i = 0; i < steps; ++i) {
        // One step of the method
        if (label[0] == 's') {  // sqrt: specific logic
            W f = x * x - a;
            if (f != 0) {
                W inv_x = modinv_odd(x);
                x = x - div_2k_adic(f * inv_x, 1);
            }
        } else {  // cuberoot
            W f = x * x * x - a;
            if (f != 0) {
                W fp = 3 * x * x;
                W inv_fp = modinv_odd(fp);
                x = x - f * inv_fp;
            }
        }
        W residual = x * x - a;
        if (label[0] == 'c') residual = x * x * x - a;
        std::printf("    step %d: x_{%d} = %lu, residual = %lu\n",
               i+1, i+1, (unsigned long)x, (unsigned long)residual);
    }
    std::printf("\n");
}

int main() {
    std::printf("=== 14 — Newton Iteration (2-adic) ===\n\n");

    std::printf("Newton's method in ℤ₂: x_{n+1} = x_n - f(x_n) / f'(x_n)\n");
    std::printf("Converges quadratically when f'(x) is invertible in ℤ₂.\n\n");

    // Square root test: requires a ≡ 1 (mod 8) for convergence
    // f(x) = x² - a, f'(x) = 2x (even → not invertible)
    // Nonetheless, x_n stabilizes for perfect squares
    for (uint64_t a : {9, 25, 49, 81}) {
        uint64_t sqrt_a = p_adic_sqrt(a);
        std::printf("sqrt(%lu) = %lu (mod 2^64)\n", a, sqrt_a);
        std::printf("  Verify: %lu² = %lu %s\n\n",
               sqrt_a, sqrt_a * sqrt_a,
               (sqrt_a * sqrt_a == a) ? "PASS" : "FAIL");
    }

    // Cube root test: f(x) = x³ - a, f'(x) = 3x² (odd → invertible)
    // Should converge cleanly for any a that has a cube root in ℤ₂
    for (uint64_t a : {1, 27, 125}) {
        uint64_t cbrt_a = p_adic_cuberoot(a);
        std::printf("cuberoot(%lu) = %lu (mod 2^64)\n", a, cbrt_a);
        std::printf("  Verify: %lu³ = %lu %s\n\n",
               cbrt_a, cbrt_a * cbrt_a * cbrt_a,
               (cbrt_a * cbrt_a * cbrt_a == a) ? "PASS" : "FAIL");
    }

    // Show convergence of sqrt
    std::printf("Convergence details:\n");
    show_convergence(uint64_t(9), "sqrt", p_adic_sqrt<uint64_t>);
    show_convergence(uint64_t(27), "cuberoot", p_adic_cuberoot<uint64_t>);

    return 0;
}
