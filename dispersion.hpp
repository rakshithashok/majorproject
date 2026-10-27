#pragma once

#include "clique.hpp"
#include "matrix.hpp"
#include <cstdint>
#include <vector>
#include <optional>

namespace fmm {

struct Point {
    std::vector<std::int64_t> x;
};

// Sizes of the cliques in the three parts of the auxiliary graph.
// P = M12 * M23 is an n^k1 x n^k2 by n^k2 x n^k3 product.
struct CliqueSplit {
    int k1, k2, k3;
};

CliqueSplit clique_split(int k);

struct DispersionResult {
    bool feasible = false;
    unsigned __int128 squared_distance = 0;
    std::vector<int> selected;
};

class FMMKDispersion {
public:
    FMMKDispersion(std::vector<Point> points, int k, std::size_t gemm_block = 64);

    // Exact k-dispersion for integral coordinates, using squared Euclidean
    // distances. The optimum squared distance is returned.
    DispersionResult solve();

private:
    struct FeasibleResult {
        bool yes = false;
        std::vector<int> witness;
    };

    std::vector<Point> points_;
    int n_ = 0;
    int k_ = 0;
    int d_ = 0;
    std::size_t gemm_block_;
    std::vector<std::vector<unsigned __int128>> dist2_;
    std::vector<unsigned __int128> candidate_thresholds_;

    void build_distances();
    void build_candidates();
    ThresholdGraph build_graph(unsigned __int128 r2) const;

    FeasibleResult fmm_clique(const ThresholdGraph& G) const;

    static void insert_unique_sorted(std::vector<int>& a, int x);
};

std::string uint128_to_string(unsigned __int128 x);

} // namespace fmm
