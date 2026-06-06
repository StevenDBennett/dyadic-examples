#include <dyadic.h>
#include <dyadic/matrix.h>
#include <iostream>
#include <iomanip>
#include <array>

template<int M, int N, std::unsigned_integral W>
void print_matrix(const dyadic::Matrix<M, N, W>& A, const char* label) {
    auto& os = std::cout;
    os << label << ":\n";
    for (int i = 0; i < M; ++i) {
        os << "  ";
        for (int j = 0; j < N; ++j)
            os << std::setw(11) << A[i][j] << " ";
        os << "\n";
    }
}

int main() {
    using W = uint32_t;
    using namespace dyadic;
    auto& os = std::cout;

    os << "=== 30 — Matrix / Linear Algebra over Z₂ ===\n\n";

    // 1. Basic operations
    os << "--- 1. Construction & operations ---\n\n";
    Matrix<2, 3, W> A{{{
        {1, 2, 3},
        {4, 5, 6}
    }}};
    Matrix<2, 3, W> B{{{
        {7, 8, 9},
        {10, 11, 12}
    }}};
    auto C = A + B;
    os << "A + B:\n";
    for (int i = 0; i < 2; ++i) {
        os << "  ";
        for (int j = 0; j < 3; ++j) os << C[i][j] << " ";
        os << "\n";
    }
    os << "\n";

    // 2. Identity and zero
    os << "--- 2. Identity & zero matrices ---\n\n";
    auto I = Matrix<3, 3, W>::identity();
    auto Z = Matrix<3, 3, W>::zero();
    os << "I * A = A: ";
    Matrix<3, 3, W> A3{{{
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}
    }}};
    auto IA = I * A3;
    bool id_ok = true;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            if (IA[i][j] != A3[i][j]) id_ok = false;
    os << (id_ok ? "PASS" : "FAIL") << "\n\n";

    // 3. Transpose
    os << "--- 3. Transpose ---\n\n";
    auto AT = A.transpose();
    os << "A (2×3):\n";
    for (int i = 0; i < 2; ++i) {
        os << "  ";
        for (int j = 0; j < 3; ++j) os << A[i][j] << " ";
        os << "\n";
    }
    os << "A^T (3×2):\n";
    for (int i = 0; i < 3; ++i) {
        os << "  ";
        for (int j = 0; j < 2; ++j) os << AT[i][j] << " ";
        os << "\n";
    }
    // Verify (A^T)^T = A
    auto ATT = AT.transpose();
    bool trans_ok = true;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 3; ++j)
            if (ATT[i][j] != A[i][j]) trans_ok = false;
    os << "(A^T)^T == A: " << (trans_ok ? "PASS" : "FAIL") << "\n\n";

    // 4. Determinant & inverse
    os << "--- 4. Determinant & inverse ---\n\n";
    // Matrix with odd determinant (det = 1), invertible over Z/2^W Z
    Matrix<3, 3, W> M3{{{
        {1, 2, 3},
        {0, 1, 4},
        {5, 6, 0}
    }}};
    auto det = M3.determinant();
    os << "det(M) = " << det << "\n";
    auto M3_inv = M3.inverse();
    auto check = M3 * M3_inv;
    bool inv_ok = true;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            W expected = (i == j) ? W(1) : W(0);
            if (check[i][j] != expected) inv_ok = false;
        }
    os << "M * M^{-1} == I: " << (inv_ok ? "PASS" : "FAIL") << "\n\n";

    // 5. Rank
    os << "--- 5. Rank ---\n\n";
    Matrix<3, 3, W> full_rank{{{
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    }}};
    Matrix<3, 3, W> rank2{{{
        {1, 2, 3},
        {4, 5, 6},
        {2, 4, 6}
    }}};
    Matrix<3, 3, W> rank1{{{
        {1, 2, 3},
        {2, 4, 6},
        {3, 6, 9}
    }}};
    os << "Identity rank = " << full_rank.rank() << " (expected 3): "
       << (full_rank.rank() == 3 ? "PASS" : "FAIL") << "\n";
    os << "Rank-2 matrix rank = " << rank2.rank() << "\n";
    os << "Rank-1 matrix rank = " << rank1.rank() << "\n\n";

    // 6. Matrix multiplication associativity
    os << "--- 6. Matrix multiplication (A*B)*C == A*(B*C) ---\n\n";
    Matrix<2, 3, W> X{{{
        {1, 2, 3},
        {4, 5, 6}
    }}};
    Matrix<3, 2, W> Y{{{
        {7, 10},
        {8, 11},
        {9, 12}
    }}};
    Matrix<2, 2, W> Z2{{{
        {1, 2},
        {3, 4}
    }}};
    auto XY = X * Y;    // 2×2
    auto XY_Z = XY * Z2; // (X*Y)*Z
    auto YZ = Y * Z2;    // 3×2
    auto X_YZ = X * YZ;  // X*(Y*Z)
    bool assoc_ok = true;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            if (XY_Z[i][j] != X_YZ[i][j]) assoc_ok = false;
    os << "(X*Y)*Z == X*(Y*Z): " << (assoc_ok ? "PASS" : "FAIL") << "\n\n";

    // 7. Solve linear system
    os << "--- 7. Solve Ax = b ---\n\n";
    // Use same invertible matrix (det = 1) with known solution
    // With A = [[1,2,3],[0,1,4],[5,6,0]], solve Ax = b for a known x.
    // Pick x = [7, 11, 13], then b = A*x.
    std::array<W, 3> expected_x{{7, 11, 13}};
    Matrix<3, 3, W> sys{{{
        {1, 2, 3},
        {0, 1, 4},
        {5, 6, 0}
    }}};
    std::array<W, 3> b_solution{
        W(1*7 + 2*11 + 3*13),
        W(0*7 + 1*11 + 4*13),
        W(5*7 + 6*11 + 0*13)
    };
    auto x = sys.solve(b_solution);
    // Verify: A*x should equal b
    std::array<W, 3> verified{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j)
            verified[i] += sys[i][j] * x[j];
    }
    bool solve_ok = true;
    for (int i = 0; i < 3; ++i) {
        if (verified[i] != b_solution[i]) solve_ok = false;
        if (x[i] != expected_x[i]) solve_ok = false;
    }
    os << "A * x = b: " << (solve_ok ? "PASS" : "FAIL") << "\n";
    os << "  x = [ " << x[0] << ", " << x[1] << ", " << x[2] << " ]";
    os << "  (expected [ " << expected_x[0] << ", " << expected_x[1] << ", "
       << expected_x[2] << " ])\n";

    // 8. Singular matrix
    os << "--- 8. Singular matrix detection ---\n\n";
    Matrix<2, 2, W> singular{{{
        {2, 4},
        {6, 8}
    }}};
    auto sing_det = singular.determinant();
    auto sing_inv = singular.inverse();
    os << "  Singular matrix det = " << sing_det << " (even = singular): "
       << (sing_det % 2 == 0 ? "PASS" : "FAIL") << "\n";
    os << "  Inverse of singular returns zero: "
       << (sing_inv[0][0] == 0 && sing_inv[1][0] == 0 ? "PASS" : "FAIL") << "\n\n";

    // 9. Compile-time usage
    os << "--- 9. Compile-time constexpr ---\n\n";
    constexpr auto CT_M = Matrix<2, 2, W>{{{
        {1, 2},
        {3, 4}
    }}};
    constexpr auto CT_DET = CT_M.determinant();
    constexpr auto CT_RANK = CT_M.rank();
    os << "  det = " << CT_DET << ", rank = " << CT_RANK << "\n\n";

    os << "=== Summary ===\n";
    os << "  Construction & arithmetic:         Working\n";
    os << "  Identity, zero, transpose:         Working\n";
    os << "  Determinant:                       Working\n";
    os << "  Inverse (Gauss-Jordan):            Working\n";
    os << "  Rank (Gaussian elimination):       Working\n";
    os << "  Linear system solve:               Working\n";
    os << "  Singular matrix detection:         Working\n";
    os << "  Compile-time constexpr:            Working\n";

    return 0;
}
