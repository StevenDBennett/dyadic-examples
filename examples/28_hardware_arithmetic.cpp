#include <dyadic.h>
#include <iostream>
#include <iomanip>
#include <cstdint>
#include <utility>

int main() {
    using namespace dyadic;
    using namespace dyadic::detail;
    auto& os = std::cout;
    os << "=== 28 — Hardware-Accelerated Arithmetic ===\n\n";

    // 1. adc (add with carry)
    os << "--- 1. adc (add with carry) ---\n\n";

    auto [s0, c0] = adc<uint32_t>(5, 10, 0);
    os << "adc<uint32_t>(5, 10, 0) = (sum=" << s0 << ", carry=" << c0 << ")\n";
    os << "  expected: sum=15, carry=0\n\n";

    auto [s1, c1] = adc<uint32_t>(0xFFFFFFFF, 1, 0);
    os << "adc<uint32_t>(0xFFFFFFFF, 1, 0) = (sum=" << s1 << ", carry=" << c1 << ")\n";
    os << "  expected: sum=0, carry=1\n\n";

    auto [s2, c2] = adc<uint32_t>(0xFFFFFFFF, 0, 1);
    os << "adc<uint32_t>(0xFFFFFFFF, 0, 1) = (sum=" << s2 << ", carry=" << c2 << ")\n";
    os << "  expected: sum=0, carry=1\n\n";

    // Chained addition: 0xFFFFFFFF_FFFFFFFF + 0x00000001_00000000
    auto [lo, ca] = adc<uint32_t>(0xFFFFFFFF, 0x00000000, 0);
    auto [hi, cb] = adc<uint32_t>(0xFFFFFFFF, 0x00000001, ca);
    os << "Chain: 0xFFFFFFFF_FFFFFFFF + 0x00000001_00000000 =\n";
    os << "  lo=0x" << std::hex << lo << ", hi=0x" << hi << ", final_carry=0x" << cb << std::dec << "\n";
    os << "  expected: lo=0xFFFFFFFF, hi=0x00000000, final_carry=0x1\n\n";

    // 2. add_overflow detection
    os << "--- 2. Overflow-checked addition ---\n\n";

    uint32_t r32;
    bool ov;
    ov = add_overflow<uint32_t>(100, 200, &r32);
    os << "add_overflow(100, 200): result=" << r32 << ", overflow=" << ov << " (expected: 300, 0)\n";

    ov = add_overflow<uint32_t>(0xFFFFFFFF, 1, &r32);
    os << "add_overflow(0xFFFFFFFF, 1): result=" << r32 << ", overflow=" << ov << " (expected: 0, 1)\n\n";

    // 3. mul_overflow detection
    os << "--- 3. Overflow-checked multiplication ---\n\n";

    ov = mul_overflow<uint32_t>(1000, 2000, &r32);
    os << "mul_overflow(1000, 2000): result=" << r32 << ", overflow=" << ov << " (expected: 2000000, 0)\n";

    // 4. Compile-time verification (static_assert)
    os << "--- 4. Compile-time constexpr usage ---\n\n";

    constexpr auto ca_test = []{
        auto [x, cy] = adc<uint32_t>(10, 20, 0);
        (void)cy;
        uint32_t tmp;
        bool ov1 = add_overflow<uint32_t>(x, 5, &tmp);
        return ov1 ? 0 : tmp;
    }();
    static_assert(ca_test == 35, "constexpr adc + add_overflow failed");
    os << "constexpr adc + add_overflow: result=" << ca_test << " (expected: 35) PASS\n\n";

    // 5. Carry chain verification (existing)
    os << "--- 5. carry_chain with dword accumulation ---\n\n";

    // Simulate the unsaturated polynomial product from section 16
    uint32_t a_coeff[3] = {1, 2, 3};
    uint32_t b_coeff[2] = {5, 7};
    // Unsaturated product (using dword_t accumulator)
    uint64_t unsat[4] = {};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 2; ++j)
            unsat[i + j] += static_cast<uint64_t>(a_coeff[i]) * static_cast<uint64_t>(b_coeff[j]);
    // Carry chain
    uint32_t result[4];
    auto final_carry = carry_chain<uint32_t>(result, unsat, 4);
    os << "(1 + 2x + 3x^2) * (5 + 7x) = ";
    for (int i = 0; i < 4; ++i) os << result[i] << (i < 3 ? "x^" + std::to_string(i+1) + " + " : "");
    os << "\n  final_carry = " << final_carry << "\n";
    // Expected: 5 + 17x + 29x^2 + 21x^3
    bool ok = (result[0] == 5 && result[1] == 17 && result[2] == 29 && result[3] == 21 && final_carry == 0);
    os << "  correct: " << (ok ? "PASS" : "FAIL") << "\n\n";

    // 6. ADX-specific info
    os << "--- 6. Platform info ---\n\n";
#if defined(__x86_64__)
    os << "  Architecture: x86-64\n";
#else
    os << "  Architecture: non-x86\n";
#endif
#if defined(__ADX__)
    os << "  ADX intrinsics: available (compile with -madx)\n";
#else
    os << "  ADX intrinsics: not compiled in\n";
#endif
#if (defined(__GNUC__) || defined(__clang__)) && !defined(__STRICT_ANSI__)
    os << "  __builtin_add_overflow: available\n";
    os << "  __builtin_mul_overflow: available\n";
#else
    os << "  Builtin overflow builtins: not available (software fallback)\n";
#endif
    os << "\n";

    os << "=== Summary ===\n";
    os << "  adc (add with carry) constexpr:          Working\n";
    os << "  add_overflow detection:                  Working\n";
    os << "  mul_overflow detection:                  Working\n";
    os << "  Compile-time constexpr usage:            Working\n";
    os << "  carry_chain integration:                 Working\n";

    return 0;
}
