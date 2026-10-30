#include <bits/stdc++.h>
using namespace std;

using u64 = uint64_t;
using u128 = __uint128_t;

static mt19937_64 rng(
    (u64)chrono::steady_clock::now().time_since_epoch().count()
);

u64 mul_mod(u64 a, u64 b, u64 mod) {
    return (u128)a * b % mod;
}

u64 pow_mod(u64 a, u64 d, u64 mod) {
    u64 r = 1;
    while (d) {
        if (d & 1) r = mul_mod(r, a, mod);
        a = mul_mod(a, a, mod);
        d >>= 1;
    }
    return r;
}

bool isPrime(u64 n) {
    if (n < 2) return false;

    // Small primes first
    for (u64 p : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL, 23ULL, 29ULL, 31ULL, 37ULL}) {
        if (n % p == 0) return n == p;
    }

    u64 d = n - 1, s = 0;
    while ((d & 1) == 0) {
        d >>= 1;
        ++s;
    }

    auto witness = [&](u64 a) -> bool {
        if (a % n == 0) return true;

        u64 x = pow_mod(a, d, n);
        if (x == 1 || x == n - 1) return true;

        for (u64 r = 1; r < s; ++r) {
            x = mul_mod(x, x, n);
            if (x == n - 1) return true;
        }
        return false;
    };

    // Deterministic for all 64-bit integers
    for (u64 a : {2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL}) {
        if (!witness(a)) return false;
    }

    return true;
}

u64 pollard_rho(u64 n) {
    if ((n & 1ULL) == 0) return 2;
    if (n % 3ULL == 0) return 3;

    uniform_int_distribution<u64> dist;

    while (true) {
        u64 c = dist(rng) % (n - 1) + 1;
        u64 x = dist(rng) % (n - 3) + 2;
        u64 y = x;
        u64 d = 1;

        auto f = [&](u64 v) {
            return (mul_mod(v, v, n) + c) % n;
        };

        while (d == 1) {
            x = f(x);
            y = f(f(y));
            u64 diff = (x > y) ? (x - y) : (y - x);
            d = std::gcd(diff, n);
        }

        if (d != n) return d;
    }
}

void factorRec(u64 n, vector<u64>& fac) {
    if (n == 1) return;
    if (isPrime(n)) {
        fac.push_back(n);
        return;
    }

    u64 d = pollard_rho(n);
    factorRec(d, fac);
    factorRec(n / d, fac);
}

u64 abs_u64(long long x) {
    // Safe even for LLONG_MIN
    return x < 0 ? (u64)(-(__int128)x) : (u64)x;
}

vector<pair<long long,int>> factorize(long long x) {
    vector<pair<long long,int>> res;

    if (x == 0) return res;   // 0 has no prime factorization
    if (x == 1) return res;   // empty factorization
    if (x == -1) {
        res.push_back({-1, 1});
        return res;
    }

    bool neg = (x < 0);
    u64 n = abs_u64(x);

    vector<u64> fac;
    factorRec(n, fac);
    sort(fac.begin(), fac.end());

    if (neg) res.push_back({-1, 1});

    for (u64 p : fac) {
        if (!res.empty() && res.back().first == (long long)p) {
            res.back().second++;
        } else {
            res.push_back({(long long)p, 1});
        }
    }

    return res;
}

int getPrimeCount(const vector<pair<long long,int>>& factors, long long p)
{
    for (const auto& [prime, cnt] : factors)
    {
        if (prime == p)
            return cnt;
    }
    return 0;
}

int mainPrime() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    long long x = 72;

    auto f = factorize(x);

    if (x == 0) {
        cout << "0 has no prime factorization\n";
        return 0;
    }

    if (f.empty()) {
        cout << x << "\n";
        return 0;
    }

    bool first = true;
    for (auto [p, cnt] : f) {
        if (!first) cout << " * ";
        first = false;

        cout << p;
        if (cnt > 1) cout << "^" << cnt;
    }
    cout << "\n";

    cout << getPrimeCount(f, 3);
    return 0;
}