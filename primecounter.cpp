#define _CRT_SECURE_NO_WARNINGS 1
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <limits>
#include <chrono>
#include <iomanip>
#include <random>
#include <cstdint>

#include "lmo_pi.h"
//#include "verifier.h"
//#include "algcmp1.h"

using u64 = std::uint64_t;

#if defined(_MSC_VER) && defined(_M_X64)
#include <intrin.h>
// (a*b) % mod via 128-bit intermediates using MSVC intrinsics
static inline u64 mul_mod_u64(u64 a, u64 b, u64 mod) {
	unsigned __int64 hi;
	unsigned __int64 lo = _umul128(a, b, &hi); // 128 = (hi:lo)
	unsigned __int64 rem;
	_udiv128(hi, lo, mod, &rem);               // rem = (hi:lo) % mod
	return rem;
}
#else
// GCC/Clang: use __int128 directly
static inline u64 mul_mod_u64(u64 a, u64 b, u64 mod) {
	return (u64)((__uint128_t)a * b % mod);
}
#endif

static inline u64 mod_pow_u64(u64 a, u64 d, u64 mod) {
	u64 r = 1 % mod;
	while (d) {
		if (d & 1) r = mul_mod_u64(r, a, mod);
		a = mul_mod_u64(a, a, mod);
		d >>= 1;
	}
	return r;
}

static bool miller_rabin_u64(u64 n) {
	if (n < 2) return false;
	static const u64 small[] = { 2,3,5,7,11,13,17,19,23,29,31,37 };
	for (u64 p : small) { 
		if (n % p == 0) return n == p; 
	}

	// n-1 = d * 2^s
	u64 d = n - 1, s = 0;
	while ((d & 1) == 0) { d >>= 1; ++s; }

	static const u64 bases[] = { 2,3,5,7,11,13,17,19,23 };
	for (u64 a : bases) {
		if (a % n == 0) continue;
		u64 x = mod_pow_u64(a, d, n);
		if (x == 1 || x == n - 1) continue;
		bool witness = true;
		for (u64 r = 1; r < s; ++r) {
			x = mul_mod_u64(x, x, n);
			if (x == n - 1) { witness = false; break; }
		}
		if (witness) return false;
	}
	return true;
}

static std::vector<u64> largest_primes_leq_n_fast(u64 n, int need = 10) {
	std::vector<u64> res;
	if (n < 2 || need <= 0) return res;
	if (n >= 2 && need > 0) {
		u64 x = (n % 2 == 0 ? n - 1 : n);
		while (true) {
			if (x == 1) break;
			if (x % 3 != 0 || x == 3) {
				if (miller_rabin_u64(x)) {
					res.push_back(x);
					if ((int)res.size() >= need) break;
				}
			}
			if (x <= 3) break;
			x -= 2;
		}
		if ((int)res.size() < need && n >= 2) {
			res.push_back(2);
		}
	}
	return res;
}

//template <typename solver1>
//void check_correct(solver1 f1) {
//	const size_t sz = 1e6;
//	for (int i = 1; i < sz; ++i) {
//		auto a1 = f1(i);
//		auto a2 = algcmp1::count_pi(i);
//		if (a1 != a2) {
//			std::cout << "WA!\tx = " << i << '\t' << a1 << '\t' << a2 << '\n';
//		}
//		if (i % 1000 == 0) {
//			std::cout << "i = " << i << std::endl;
//		}
//	}
//	std::cout << "Over." << std::endl;
//}
//
//template <typename solver1, typename solver2>
//void check_time(solver1 f1, solver2 f2, size_t testsize = 500) {
//	static std::random_device rd;
//	static std::mt19937_64 gen(rd());
//	//std::uniform_int_distribution<std::uint64_t> dist(2, (std::uint64_t)(1'000'000'000'000ULL));
//	//std::uniform_int_distribution<std::uint64_t> dist(2, (std::uint64_t)(4'000'000ULL));
//	std::uniform_int_distribution<std::uint64_t> dist((std::uint64_t)(900'000'000'000ULL), (std::uint64_t)(1'100'000'000'000ULL));
//	auto t1 = std::chrono::steady_clock::now();
//	auto t2 = std::chrono::steady_clock::now();
//	auto t3 = std::chrono::steady_clock::now();
//	auto t4 = std::chrono::steady_clock::now();
//	double dt1 = 0;
//	double dt2 = 0;
//
//	std::vector<u64> testx(testsize);
//	std::vector<u64> ans1(testsize), ans2(testsize);
//	bool flag = 0;
//
//	for (int i = 0; i < testsize; ++i) {
//		testx[i] = dist(gen);
//	}
//
//	for (int i = 0; i < testsize; ++i) {
//		t1 = std::chrono::steady_clock::now();
//		ans1[i] = f1(testx[i]);
//		t2 = std::chrono::steady_clock::now();
//		dt1 += std::chrono::duration<double, std::milli>(t2 - t1).count();
//	}
//	std::cout << "Elapsed_1: " << dt1 / testsize << " ms\n";
//
//	for (int i = 0; i < testsize; ++i) {
//		t3 = std::chrono::steady_clock::now();
//		ans2[i] = f2(testx[i]);
//		t4 = std::chrono::steady_clock::now();
//		dt2 += std::chrono::duration<double, std::milli>(t4 - t3).count();
//	}
//	std::cout << "Elapsed_2: " << dt2 / testsize << " ms\n";
//
//	for (int i = 0; i < testsize; ++i) {
//		if (ans1[i] != ans2[i]) {
//			std::cout << "WA!\tx = " << testx[i] << "  \t" << ans1[i] << '\t' << ans2[i] << '\n';
//			flag = 1;
//		}
//	}
//
//	if (flag == 0) {
//		std::cout << "AC." << std::endl;
//	}
//}

int main(int argc, char* argv[]) {

	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr);
	

	if (argc < 2) { // for test
		//omp_set_num_threads(6);
		//getchar(); // input any key to continue
		//check_time(lmo_pi::count_pi,verifier::prime_pi_sieve);
		//check_time(lmo_pi::count_pi, static_cast<int64_t(*)(int64_t)>(primecount::pi));
		//check_time(lmo_pi_modify::count_pi, lmo_pi::count_pi);
		//check_correct(lmo_pi::count_pi);
		//check_time(lmo_pi::count_pi, lmo_pi::count_pi);
	}
	else { // release
		//omp_set_num_threads(8);
		u64 x = std::strtoull(argv[1], nullptr, 10);
		auto t1 = std::chrono::steady_clock::now();
		auto ans = lmo_pi::count_pi(x);
		/*auto t2 = std::chrono::steady_clock::now();
		auto dt = std::chrono::duration<double, std::milli>(t2 - t1).count();
		std::cout << "ans = " << ans << "\nTime: " << dt << " ms" << std::endl;*/
		std::cout << ans << '\n';

		//t1 = std::chrono::steady_clock::now();
		auto top = largest_primes_leq_n_fast(x, 10);
		//t2 = std::chrono::steady_clock::now();

		//auto alldt = dt;
		//dt = std::chrono::duration<double, std::milli>(t2 - t1).count();
		//std::cout << "\n";
		for (size_t i = 0; i < top.size(); ++i) {
			std::cout << top[i] << '\n';
		}

		auto t2 = std::chrono::steady_clock::now();
		auto dt = std::chrono::duration<double, std::milli>(t2 - t1).count();

		//alldt += dt;
		std::cout << "Time: " << dt << " ms" << std::endl;
		//std::cout << "\nAll Time: " << alldt << " ms" << std::endl;
	}

	return 0;
}
