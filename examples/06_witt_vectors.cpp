#include <dyadic.h>
#include <cstdio>

template<int N, std::unsigned_integral W>
void print_witt(const char* label, const dyadic::WittVector<N, W>& w) {
    std::printf("  %-14s = {", label);
    for (int i = 0; i < N; ++i) std::printf("%s%lu", i ? ", " : "", (unsigned long)w[i]);
    std::printf("}\n");
}

int main() {
    using namespace dyadic;

    std::printf("=== 06 — Witt Vectors ===\n\n");

    WittVector<4, uint32_t> a{{1, 2, 3, 4}};
    WittVector<4, uint32_t> b{{5, 6, 7, 8}};

    print_witt("a", a);
    print_witt("b", b);
    std::printf("\n");

    auto sum = a + b;
    print_witt("a + b (Witt)", sum);
    auto ghost_check_add = [&]() {
        auto ga = a.ghost_vector();
        auto gb = b.ghost_vector();
        auto gs = sum.ghost_vector();
        std::printf("  Ghost add check: G(a+b)_j == G(a)_j + G(b)_j for all j: ");
        bool ok = true;
        for (int j = 0; j < 4; ++j) {
            if (gs[j] != ga[j] + gb[j]) { ok = false; break; }
        }
        std::printf("%s\n", ok ? "PASS" : "FAIL");
    };
    ghost_check_add();

    auto prod = a * b;
    print_witt("a × b (Witt)", prod);
    auto ghost_check_mul = [&]() {
        auto ga = a.ghost_vector();
        auto gb = b.ghost_vector();
        auto gp = prod.ghost_vector();
        std::printf("  Ghost mul check: G(a×b)_j == G(a)_j × G(b)_j for all j: ");
        bool ok = true;
        for (int j = 0; j < 4; ++j) {
            if (gp[j] != ga[j] * gb[j]) { ok = false; break; }
        }
        std::printf("%s\n", ok ? "PASS" : "FAIL");
    };
    ghost_check_mul();

    std::printf("\nFrobenius F(a)_i = a_i²:\n");
    auto fa = a.frobenius();
    print_witt("F(a)", fa);

    std::printf("\nVerschiebung V(a) = (0, a₀, a₁, …):\n");
    auto va = a.verschiebung();
    print_witt("V(a)", va);

    std::printf("\nIdentity: F(V(a)) == V(F(a)):\n");
    auto fv = a.FV();
    auto vf = a.VF();
    print_witt("F(V(a))", fv);
    print_witt("V(F(a))", vf);
    bool fv_eq = true;
    for (int i = 0; i < 4; ++i) if (fv[i] != vf[i]) fv_eq = false;
    std::printf("  FV == VF: %s\n\n", fv_eq ? "PASS" : "FAIL");

    std::printf("Ghost vectors:\n");
    auto ga = a.ghost_vector();
    std::printf("  G(a) = {%u, %u, %u, %u}\n", ga[0], ga[1], ga[2], ga[3]);

    std::printf("\nRing axioms (distributivity):\n");
    WittVector<4, uint32_t> c{{9, 10, 11, 12}};
    auto lhs = a * (b + c);
    auto rhs = a * b + a * c;
    print_witt("a×(b+c)", lhs);
    print_witt("a×b + a×c", rhs);
    bool dist = true;
    for (int i = 0; i < 4; ++i) if (lhs[i] != rhs[i]) dist = false;
    std::printf("  a×(b+c) == a×b + a×c: %s\n", dist ? "PASS" : "FAIL");

    return 0;
}
