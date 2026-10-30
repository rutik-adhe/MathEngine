#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <iostream>
#include <mutex>
#include <random>
#include <thread>
#include <vector>

using boost::multiprecision::cpp_int;

static cpp_int abs_cpp(cpp_int x) {
    return x < 0 ? -x : x;
}

static cpp_int gcd_cpp(cpp_int a, cpp_int b) {
    a = abs_cpp(a);
    b = abs_cpp(b);
    while (b != 0) {
        cpp_int t = a % b;
        a = b;
        b = t;
    }
    return a;
}

static cpp_int mod_pow(cpp_int a, cpp_int d, const cpp_int& mod) {
    cpp_int r = 1;
    a %= mod;
    while (d > 0) {
        if ((d & 1) != 0) r = (r * a) % mod;
        a = (a * a) % mod;
        d >>= 1;
    }
    return r;
}

static bool is_probable_prime(const cpp_int& n) {
    if (n < 2) return false;

    static const int small_primes[] = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37
    };

    for (int p : small_primes) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }

    cpp_int d = n - 1;
    unsigned s = 0;
    while ((d & 1) == 0) {
        d >>= 1;
        ++s;
    }

    auto witness = [&](cpp_int a) -> bool {
        if (a % n == 0) return true;
        cpp_int x = mod_pow(a, d, n);
        if (x == 1 || x == n - 1) return true;

        for (unsigned r = 1; r < s; ++r) {
            x = (x * x) % n;
            if (x == n - 1) return true;
        }
        return false;
    };

    for (int p : small_primes) {
        if (!witness(p)) return false;
    }

    return true;
}

static cpp_int pollard_rho_worker(const cpp_int& n, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<uint64_t> dist;

    while (true) {
        cpp_int c = 1 + cpp_int(dist(rng)) % (n - 1);
        cpp_int x = 2 + cpp_int(dist(rng)) % (n - 3);
        cpp_int y = x;
        cpp_int d = 1;

        auto f = [&](const cpp_int& v) {
            return (v * v + c) % n;
        };

        while (d == 1) {
            x = f(x);
            y = f(f(y));
            cpp_int diff = (x > y) ? (x - y) : (y - x);
            d = gcd_cpp(diff, n);
        }

        if (d != n) return d;
    }
}

static cpp_int pollard_rho_parallel(const cpp_int& n,
                                    unsigned thread_count = std::thread::hardware_concurrency()) {
    if ((n & 1) == 0) return 2;
    if (n % 3 == 0) return 3;

    if (thread_count == 0) thread_count = 4;

    std::atomic<bool> found(false);
    std::mutex out_mutex;
    cpp_int answer = 0;

    auto now_seed = []() -> uint64_t {
        return (uint64_t)std::chrono::steady_clock::now().time_since_epoch().count();
    };

    auto worker = [&](uint64_t seed) {
        std::mt19937_64 local_rng(seed);

        while (!found.load(std::memory_order_relaxed)) {
            cpp_int d = pollard_rho_worker(n, local_rng());
            if (d != 1 && d != n) {
                std::lock_guard<std::mutex> lock(out_mutex);
                if (!found.exchange(true)) {
                    answer = d;
                }
                return;
            }
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(thread_count);

    for (unsigned i = 0; i < thread_count; ++i) {
        threads.emplace_back(worker, now_seed() ^ (0x9e3779b97f4a7c15ULL + i));
    }

    for (auto& t : threads) t.join();
    return answer;
}

static void factor_rec(const cpp_int& n, std::vector<cpp_int>& fac) {
    if (n == 1) return;

    if (is_probable_prime(n)) {
        fac.push_back(n);
        return;
    }

    cpp_int d = pollard_rho_parallel(n);
    factor_rec(d, fac);
    factor_rec(n / d, fac);
}

std::vector<std::pair<cpp_int, int>> factorize(cpp_int n) {
    std::vector<std::pair<cpp_int, int>> res;

    if (n == 0) {
        throw std::invalid_argument("0 has no prime factorization");
    }

    if (n < 0) {
        res.push_back({-1, 1});
        n = -n;
    }

    if (n == 1) return res;

    std::vector<cpp_int> fac;
    factor_rec(n, fac);
    std::sort(fac.begin(), fac.end());

    for (const auto& p : fac) {
        if (!res.empty() && res.back().first == p) {
            ++res.back().second;
        } else {
            res.push_back({p, 1});
        }
    }

    return res;
}

void main2() {
    cpp_int n("848539221987161640771731");

    auto f = factorize(n);

    for (size_t i = 0; i < f.size(); ++i) {
        if (i) std::cout << " * ";
        std::cout << f[i].first;
        if (f[i].second > 1) std::cout << "^" << f[i].second;
    }
    std::cout << '\n';
    
}