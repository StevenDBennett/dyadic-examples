#include "dyadic.h"
#include <cstdio>
#include <array>

int main() {
    using namespace dyadic;

    std::printf("=== 09 — Carry Chain C = (I−N)⁻¹ ===\n\n");

    std::printf("The carry chain processes an array of dword values into\n");
    std::printf("normalized words, propagating carries in one pass.\n\n");

    {
        std::printf("Simple case: all values fit in a single word\n");
        std::array<uint64_t, 4> input{1, 2, 3, 4};
        std::array<uint64_t, 4> result{};
        dword_t<uint64_t> dword_input[4]{1, 2, 3, 4};
        auto carry = carry_chain(result.data(), dword_input, 4);
        std::printf("  Input:  {1, 2, 3, 4}\n");
        std::printf("  Result: {%lu, %lu, %lu, %lu}\n", result[0], result[1], result[2], result[3]);
        std::printf("  Carry:  %lu (no overflow expected)\n\n", (unsigned long)carry);
    }

    {
        std::printf("Overflow case: large values that generate carries\n");
        dword_t<uint64_t> dword_input[4]{
            dword_t<uint64_t>(1) << 63,
            dword_t<uint64_t>(1) << 63,
            dword_t<uint64_t>(1) << 63,
            dword_t<uint64_t>(1) << 63
        };
        std::array<uint64_t, 4> result{};
        auto carry = carry_chain(result.data(), dword_input, 4);
        std::printf("  Input:  {2⁶³, 2⁶³, 2⁶³, 2⁶³}\n");
        std::printf("  Result: {");
        for (int i = 0; i < 4; ++i) std::printf("%s0x%lx", i ? ", " : "", result[i]);
        std::printf("}\n");
        std::printf("  Carry:  0x%lx\n\n", (unsigned long)carry);
    }

    {
        std::printf("Maximal input: all 0xFFFF...FFFF\n");
        dword_t<uint64_t> dword_input[4]{
            dword_t<uint64_t>::all_ones(),
            dword_t<uint64_t>::all_ones(),
            dword_t<uint64_t>::all_ones(),
            dword_t<uint64_t>::all_ones()
        };
        std::array<uint64_t, 4> result{};
        auto carry = carry_chain(result.data(), dword_input, 4);
        std::printf("  Input:  {2¹²⁸−1, 2¹²⁸−1, 2¹²⁸−1, 2¹²⁸−1}\n");
        std::printf("  Result: {%lu, %lu, %lu, %lu}\n", result[0], result[1], result[2], result[3]);
        std::printf("  Carry:  %lu\n\n", (unsigned long)carry);
    }

    std::printf("Idempotence: C(C(x)) == C(x)\n");
    {
        dword_t<uint64_t> input[3]{
            (dword_t<uint64_t>(1) << 63) + dword_t<uint64_t>(42),
            dword_t<uint64_t>(1) << 62,
            dword_t<uint64_t>::all_ones()
        };
        std::array<uint64_t, 3> first{};
        auto c1 = carry_chain(first.data(), input, 3);

        dword_t<uint64_t> first_dword[3]{first[0], first[1], first[2]};
        std::array<uint64_t, 3> second{};
        auto c2 = carry_chain(second.data(), first_dword, 3);

        bool idempotent = true;
        for (int i = 0; i < 3; ++i) if (first[i] != second[i]) idempotent = false;
        std::printf("  C(C(x)) == C(x): %s (carries: %lu -> %lu)\n",
            idempotent ? "PASS" : "FAIL",
            (unsigned long)c1, (unsigned long)c2);
    }

    return 0;
}
