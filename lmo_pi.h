#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <limits>
#include <omp.h>


namespace lmo_pi {

    using ll = long long;

    // 更快地开方
    inline int isqrt(ll n) {
        return static_cast<int>(std::sqrt((double)n));
    }

    // 将奇数映射到下标
    inline ll half(ll n) {
        return (n - 1) >> 1;
    }

    // 更快地除法
    inline ll divide(ll n, ll base) {
        return static_cast<ll>(
            static_cast<double>(n) / static_cast<double>(base)
            );
    }

    // 处理一个素数 p 的整个块: 筛 skips，并更新 roughs/larges/smalls/pCnt
    inline bool process_prime(
        int p,
        int lim,
        ll n,
        int& vsz,
        std::vector<int>& smalls,
        std::vector<int>& roughs,
        std::vector<ll>& larges,
        std::vector<char>& skips,
        int& pCnt
    ) {
        int p2 = p * p;
        // p^2 > sqrt(n) ⇒ p > n^(1/4)
        if (p2 > lim)
            return false;

        skips[p] = 1;
        for (int x = p2; x <= lim; x += (p << 1))
            skips[x] = 1;

        int ns = 0;
        for (int i = 0; i < vsz; ++i) {
            int cur = roughs[i];
            if (skips[cur]) continue;

            ll d = 1LL * cur * p;
            if (d <= lim)
                larges[ns] = larges[i] - larges[smalls[d >> 1] - pCnt] + pCnt;
            else
                larges[ns] = larges[i] - smalls[half(divide(n, d))] + pCnt;

            roughs[ns++] = cur;
        }
        vsz = ns;

        int cx = static_cast<int>(half(lim));
        for (int cy = ((lim / p) - 1) | 1; cy >= p; cy -= 2) {
            int cur = smalls[cy >> 1] - pCnt;
            for (int cz = (cy * p) >> 1; cz <= cx; --cx)
                smalls[cx] -= cur;
        }

        ++pCnt;
        return true;
    }

    inline ll count_pi(ll n) {
        if (n <= 1) return 0LL;
        if (n == 2) return 1LL;

        const int lim = isqrt(n);
        int vsz = (lim + 1) >> 1;
        //std::cout << "vsz = " << vsz << std::endl;

        std::vector<int> smalls(vsz);
        for (int i = 0; i < vsz; ++i)
            smalls[i] = i;

        std::vector<int> roughs(vsz);
        for (int i = 0; i < vsz; ++i)
            roughs[i] = (i << 1) | 1;

        std::vector<ll> larges(vsz);
        for (int i = 0; i < vsz; ++i)
            larges[i] = ((n / ((i << 1) | 1)) - 1) >> 1;

        std::vector<char> skips(lim + 1, 0);
        int pCnt = 0;

        bool stop = false;

        // 先按原逻辑处理 p = 3, 5（如果在范围内且未被标记）
        if (3 <= lim && !skips[3]) {
            if (!process_prime(3, lim, n, vsz, smalls, roughs, larges, skips, pCnt))
                stop = true;
        }
        if (!stop && 5 <= lim && !skips[5]) {
            if (!process_prime(5, lim, n, vsz, smalls, roughs, larges, skips, pCnt))
                stop = true;
        }

        // 然后使用 30-轮跳过 2,3,5 的倍数。
        // 只枚举与 30 互素的候选：7,11,13,17,19,23,29,31,37,...
        if (!stop) {
            static const int wheel_steps[8] = { 4, 2, 4, 2, 4, 6, 2, 6 };
            int step_idx = 0;
            int p = 7;

            while (p <= lim) {
                if (!skips[p]) {
                    if (!process_prime(p, lim, n, vsz, smalls, roughs, larges, skips, pCnt)) {
                        stop = true;
                        break;
                    }
                }
                p += wheel_steps[step_idx];
                ++step_idx;
                if (step_idx == 8)
                    step_idx = 0;
            }
        }

        // ==== 组合公式 ====
        ll result = larges[0] + 1LL * (vsz + ((pCnt - 1) << 1)) * (vsz - 1) / 2;
        for (int i = 1; i < vsz; ++i)
            result -= larges[i];

        // 找出最后需要参与 S2 的 i 上界（仍然是顺序扫描）
        int i_limit = vsz;
        for (int i = 1; i < vsz; ++i) {
            int q = roughs[i];
            ll m = n / q;
            int e = smalls[half(m / q)] - pCnt;
            if (e < i + 1) {
                i_limit = i;
                break;
            }
        }

        // S2 双重和：保持原来的顺序实现和 double divide
#pragma omp parallel for schedule(dynamic) reduction(+:result)
        for (int i = 1; i < i_limit; ++i) {
            int q = roughs[i];
            ll m = n / q;
            int e = smalls[half(m / q)] - pCnt;
            if (e < i + 1) continue;

            ll t = 0;
            for (int j = i + 1; j <= e; ++j)
                t += smalls[half(divide(m, roughs[j]))];

            result += t - 1LL * (e - i) * (pCnt + i - 1);
        }

        return result + 1;
    }

} 
