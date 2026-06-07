#include <dyadic.h>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <string>

using namespace dyadic;

// Generate data files for visualization and analysis

int main() {
    std::printf("=== 25 — CSV Data Output for Visualization ===\n\n");

    using W = uint64_t;

    // 1. Overflow analysis for Taylor basis
    {
        std::ofstream ofs("taylor_overflow.csv");
        ofs << "degree,max_coefficient,overflows\n";

        for (int degree = 2; degree <= 20; ++degree) {
            W max_coeff = W(0) - 1;
            bool overflow = false;
            for (int k = 0; k <= degree; ++k) {
                W fact = 1;
                for (int i = 2; i <= k; ++i) fact *= i;
                if (max_coeff / fact < max_coeff) {
                    overflow = true;
                    break;
                }
            }
            if (!overflow) {
                ofs << degree << "," << max_coeff << ",no\n";
            } else {
                ofs << degree << "," << max_coeff << ",yes\n";
            }
        }
        std::printf("Generated: taylor_overflow.csv\n");
    }

    // 2. Growth of polynomial degrees under composition
    {
        std::ofstream ofs("degree_growth.csv");
        ofs << "iteration,degree\n";
        int deg = 1;
        for (int i = 0; i < 10; ++i) {
            ofs << i << "," << deg << "\n";
            deg = deg * 2; // Doubling (e.g., f(f(x)))
        }
        std::printf("Generated: degree_growth.csv\n");
    }

    // 3. Precision vs word size
    {
        std::ofstream ofs("precision_vs_wordsize.csv");
        ofs << "word_bits,max_taylor_degree,max_witt_length\n";
        for (int bits : {8, 16, 32, 64}) {
            int max_taylor = 0;
            W fact = 1;
            for (int k = 0; ; ++k) {
                if (k > 1) fact *= k;
                if (fact > W(1) << (bits - 1)) break;
                max_taylor = k;
            }
            int max_witt = bits; // Approximate
            ofs << bits << "," << max_taylor << "," << max_witt << "\n";
        }
        std::printf("Generated: precision_vs_wordsize.csv\n");
    }

    // 4. Stirling numbers output
    {
        std::ofstream ofs("stirling_numbers.csv");
        ofs << "n,k,S2(n,k),s1(n,k)\n";
        for (int n = 0; n < 10; ++n) {
            for (int k = 0; k <= n; ++k) {
                W s2 = stirling_2<W>(n, k);
                W s1 = stirling_1<W>(n, k);
                ofs << n << "," << k << "," << s2 << "," << s1 << "\n";
            }
        }
        std::printf("Generated: stirling_numbers.csv\n");
    }

    // 5. Ghost map values for Witt vectors
    {
        std::ofstream ofs("witt_ghosts.csv");
        ofs << "index,component,ghost_value\n";
        for (int len : {1, 2, 3, 4}) {
            WittVector<4, W> w;
            for (int i = 0; i < 4; ++i) w[i] = static_cast<W>(i + 1);
            auto ghosts = w.ghost_vector();
            for (int j = 0; j < 4; ++j) {
                ofs << len << "," << j << "," << ghosts[j] << "\n";
            }
        }
        std::printf("Generated: witt_ghosts.csv\n");
    }

    std::printf("\nAll CSV files have been generated in the current directory.\n");
    std::printf("You can plot them using Python/Matplotlib or any spreadsheet software.\n\n");

    std::printf("Suggested plotting commands (Python):\n");
    std::printf("  import pandas as pd\n");
    std::printf("  import matplotlib.pyplot as plt\n\n");
    std::printf("  df = pd.read_csv('stirling_numbers.csv')\n");
    std::printf("  df_pivot = df.pivot(index='n', columns='k', values='S2(n,k)')\n");
    std::printf("  plt.imshow(df_pivot, cmap='viridis')\n");
    std::printf("  plt.colorbar(label='S2(n,k)')\n");
    std::printf("  plt.show()\n\n");

    return 0;
}
