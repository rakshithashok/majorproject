#include "dispersion.hpp"
#include <iostream>
#include <string>

using namespace fmm;

int main() {
    // Input format:
    // n d k
    // x_1 ... x_d
    // ...
    // x_n ... x_d
    //
    // Coordinates are integral (signed 64-bit).
    // Output is the optimum squared minimum distance and one optimal set.

    int n, d, k;
    if (!(std::cin >> n >> d >> k)) {
        std::cerr << "usage: n d k followed by n integral d-dimensional points\n";
        return 1;
    }

    if (n <= 0 || d <= 0 || k < 2 || k > n) {
        std::cerr << "invalid n, d, or k\n";
        return 1;
    }

    std::vector<Point> points(n);
    for (int i = 0; i < n; ++i) {
        points[i].x.resize(d);
        for (int j = 0; j < d; ++j) {
            if (!(std::cin >> points[i].x[j])) {
                std::cerr << "invalid/incomplete point data\n";
                return 1;
            }
        }
    }

    try {
        FMMKDispersion solver(points, k);
        const auto ans = solver.solve();

        if (!ans.feasible) {
            std::cout << "NO SOLUTION\n";
            return 0;
        }

        std::cout << "optimal_squared_distance "
                  << uint128_to_string(ans.squared_distance) << '\n';
        std::cout << "selected_indices";
        for (int v : ans.selected) std::cout << ' ' << v;
        std::cout << '\n';

        std::cout << "selected_points\n";
        for (int v : ans.selected) {
            for (int j = 0; j < d; ++j) std::cout << (j ? " " : "") << points[v].x[j];
            std::cout << '\n';
        }
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
