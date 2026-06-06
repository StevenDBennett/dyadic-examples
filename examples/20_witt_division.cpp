#include <dyadic.h>
#include <cstdio>

using namespace dyadic;

template<int N, std::unsigned_integral W>
void print_witt(const char* label, const WittVector<N, W>& w) {
    std::printf("  %-14s = {", label);
    for (int i = 0; i < N; ++i) std::printf("%s%lu", i ? ", " : "", (unsigned long)w[i]);
    std::printf("}\n");
}

int main() {
    std::printf("=== 20 — Witt Division ===\n\n");

    using W = uint64_t;
    constexpr int N = 4;

    // Division in the Witt ring via ghost map inversion.
    // Given a Witt vector a (with a[0] odd = unit), a^{-1} is
    // computed by inverting each ghost value then recovering Witt components.

    WittVector<N, W> a{{3, 5, 7, 9}};
    print_witt("a", a);

    // Compute inverse via ghost map
    auto inv = witt_inverse(a);
    print_witt("a^{-1}", inv);

    // Ghost verification
    auto ga = a.ghost_vector();
    auto gi = inv.ghost_vector();
    std::printf("\nGhost values G(a) = {%lu, %lu, %lu, %lu}\n",
           (unsigned long)ga[0], (unsigned long)ga[1],
           (unsigned long)ga[2], (unsigned long)ga[3]);
    std::printf("Ghost values G(a^{-1}) = {%lu, %lu, %lu, %lu}\n",
           (unsigned long)gi[0], (unsigned long)gi[1],
           (unsigned long)gi[2], (unsigned long)gi[3]);

    std::printf("\nVerification: ghost(a^{-1})_j = ghost(a)_j^{-1}:\n");
    bool ghost_ok = true;
    for (int j = 0; j < N; ++j) {
        W prod = ga[j] * gi[j];
        bool ok = (prod == 1);
        std::printf("  G_%d * G_%d^{-1} = %lu * %lu = %lu  %s\n",
               j, j, (unsigned long)ga[j], (unsigned long)gi[j],
               (unsigned long)prod, ok ? "PASS" : "FAIL");
        if (!ok) ghost_ok = false;
    }

    // Verify: a * a^{-1} = τ(1) = {1, 0, 0, 0}
    auto prod = a * inv;
    print_witt("a * a^{-1}", prod);
    bool inv_ok = (prod[0] == 1 && prod[1] == 0 && prod[2] == 0 && prod[3] == 0);
    std::printf("  a * a^{-1} == τ(1): %s\n\n", inv_ok ? "PASS" : "FAIL");

    // Additional examples
    std::printf("More examples:\n");
    for (W unit_val : {W(1), W(5), W(9), W(255)}) {
        WittVector<N, W> u{{unit_val, 0, 0, 0}};
        auto u_inv = witt_inverse(u);
        auto u_prod = u * u_inv;
        bool ok = (u_prod[0] == 1 && u_prod[1] == 0 &&
                   u_prod[2] == 0 && u_prod[3] == 0);
        std::printf("  a = {%lu, 0, 0, 0}: a * a^{-1} = τ(1)  %s\n",
               (unsigned long)unit_val, ok ? "PASS" : "FAIL");
    }

    return 0;
}
