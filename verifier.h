#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <limits>

namespace verifier {

    /*@brief  Return floor(sqrt(n)) using 64-bit safe fix-up.
      @param  n  Nonnegative 64-bit integer.
      @return floor(sqrt(n)).
    */
    static inline std::uint64_t isqrt_u64(std::uint64_t n) {
        if (n == 0) return 0;
        long double d = std::sqrt(static_cast<long double>(n));
        std::uint64_t r = static_cast<std::uint64_t>(d);
        while ((r + 1) <= n / (r + 1)) ++r;
        while (r > n / r) --r;
        return r;
    }

    /*@brief  Simple sieve up to N (inclusive) to generate base primes for segmented sieve.
      @param  N        Upper bound (inclusive), typically floor(sqrt(x)).
      @return Vector of all primes <= N.
    */
    static std::vector<int> base_primes_up_to(std::uint64_t N) {
        std::size_t n = static_cast<std::size_t>(N);
        std::vector<uint8_t> isprime(n + 1, 1);
        if (n >= 0) isprime[0] = 0;
        if (n >= 1) isprime[1] = 0;
        std::uint64_t R = static_cast<std::uint64_t>(std::sqrt(static_cast<long double>(N)));
        for (std::uint64_t p = 2; p <= R; ++p) {
            if (!isprime[static_cast<std::size_t>(p)]) continue;
            std::uint64_t start = p * p;
            for (std::uint64_t j = start; j <= N; j += p) {
                isprime[static_cast<std::size_t>(j)] = 0;
            }
        }
        std::vector<int> primes;
        for (std::size_t i = 2; i <= n; ++i) if (isprime[i]) primes.push_back(static_cast<int>(i));
        return primes;
    }

    /*@brief  Count primes in [L, R] inclusive using the provided base primes.
      @param  L            Segment lower bound (inclusive).
      @param  R            Segment upper bound (inclusive). Requires R>=L.
      @param  base_primes  All primes <= floor(sqrt(R)).
      @return Number of primes in [L, R].
    */
    static std::uint64_t count_primes_in_segment(std::uint64_t L,
        std::uint64_t R,
        const std::vector<int>& base_primes) {
        std::size_t len = static_cast<std::size_t>(R - L + 1);
        std::vector<uint8_t> isprime(len, 1);
        if (L == 0) {
            if (len >= 1) isprime[0] = 0;
            if (len >= 2) isprime[1] = 0;
        }
        else if (L == 1) {
            isprime[0] = 0;
        }
        for (int p : base_primes) {
            std::uint64_t pp = static_cast<std::uint64_t>(p);
            std::uint64_t p2 = pp * pp;
            if (p2 > R) break;
            std::uint64_t start = (L + pp - 1) / pp * pp;
            if (start < p2) start = p2;
            for (std::uint64_t j = start; j <= R; j += pp) {
                isprime[static_cast<std::size_t>(j - L)] = 0;
            }
        }
        std::uint64_t cnt = 0;
        for (uint8_t v : isprime) cnt += (v != 0);
        return cnt;
    }

    /*@brief  Segmented sieve to compute π(x) with optional OpenMP parallelism.
      @param  x                 Upper bound x (count primes <= x).
      @param  segment_size      Segment size in number of integers (default: 1<<22).
      @return π(x).
    */
    std::uint64_t prime_pi_sieve(std::uint64_t x) {
        const std::uint64_t segment_size = (1ull << 22);
        if (x < 2) return 0;
        std::uint64_t root = isqrt_u64(x);
        std::vector<int> base = base_primes_up_to(root);

        std::uint64_t seg = std::max<std::uint64_t>(segment_size, 1ull << 20);
        std::uint64_t blocks = (x + seg) / seg;

        std::uint64_t total = 0;

#pragma omp parallel for reduction(+:total) schedule(dynamic)
        for (std::int64_t b = 0; b < static_cast<std::int64_t>(blocks); ++b) {
            std::uint64_t L = static_cast<std::uint64_t>(b) * seg;
            if (L < 2) L = 2;
            std::uint64_t R = std::min<std::uint64_t>(x, L + seg - 1);
            if (L > R) continue;

            std::uint64_t cnt = count_primes_in_segment(L, R, base);
            total += cnt;
        }

        return total;
    }

    

} // namespace verifier