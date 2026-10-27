# majorproject

Exact k-dispersion (select k of n points maximizing the minimum pairwise
Euclidean distance) via the FMM-based k-clique reduction used in
K. Chen, A. Dumitrescu, *A Couple of Simple Algorithms for k-Dispersion*
(arXiv:2511.00692), Section 3 / Theorem 2.

## Algorithm

1. Compute all squared pairwise distances (integral coordinates, exact).
2. The optimum is one of the pairwise distances. Binary search over the sorted
   distinct values r: build the threshold graph G_r (edge iff dist >= r);
   a k-set with min distance >= r exists iff G_r contains a K_k.
3. K_k detection (Nešetřil–Poljak / Eisenbrand–Grandoni): split
   k = k1 + k2 + k3, let V_i be the k_i-cliques of G_r, connect two clique
   nodes iff their union is a clique. G_r has a K_k iff this 3-partite graph
   has a triangle, found from P = M12 · M23 checked against M13.

The split follows the rectangular products used in the proof of Theorem 2
(`clique_split` in `dispersion.cpp`):

| k | (k1, k2, k3) | exponent e(k) |
|---|--------------|---------------|
| 4 | (1, 1, 2) | ω(1,1,2) < 3.251 |
| 5 | (2, 1, 2) | ω(2,1,2) < 4.086 |
| 6 | (2, 2, 2) | 2ω < 4.744 |
| 7 | (2, 3, 2) | ω(2,3,2) < 5.590 |
| 8 | (2, 2, 4) | 2ω(1,1,2) < 6.502 |

k = 2, 3 use (0,1,1) and (1,1,1). For k ≥ 9 (not covered by the paper) the
balanced split ⌊k/3⌋, ⌈(k−1)/3⌉, ⌈k/3⌉ is used.

## Matrix multiplication backend

`multiply_blocked` in `matrix.cpp` is an exact, cache-blocked classical
product. The reduction is exact, but the asymptotic bounds above only hold if
this routine is replaced by a fast (rectangular) matrix multiplication
algorithm; the paper does not specify one.

## Build and run

```
g++ -O2 -std=c++17 *.cpp -o dispersion
./dispersion < input.txt
```

Input: `n d k` followed by n lines of d integer coordinates.
