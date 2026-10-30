#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

using boost::multiprecision::cpp_int;

static std::mt19937_64 rng(
    (uint64_t)std::chrono::steady_clock::now().time_since_epoch().count());

cpp_int abs_cpp(cpp_int x)
{
    return x < 0 ? -x : x;
}

cpp_int gcd_cpp(cpp_int a, cpp_int b)
{
    a = abs_cpp(a);
    b = abs_cpp(b);
    while (b != 0)
    {
        cpp_int t = a % b;
        a = b;
        b = t;
    }
    return a;
}

cpp_int mod_pow(cpp_int a, cpp_int d, const cpp_int &mod)
{
    cpp_int r = 1;
    a %= mod;
    while (d > 0)
    {
        if ((d & 1) != 0)
            r = (r * a) % mod;
        a = (a * a) % mod;
        d >>= 1;
    }
    return r;
}

// Probabilistic Miller-Rabin for cpp_int
bool is_probable_prime(const cpp_int &n, int rounds = 20)
{
    if (n < 2)
        return false;

    static const int small_primes[] = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};

    for (int p : small_primes)
    {
        if (n == p)
            return true;
        if (n % p == 0)
            return false;
    }

    cpp_int d = n - 1;
    unsigned s = 0;
    while ((d & 1) == 0)
    {
        d >>= 1;
        ++s;
    }

    auto witness = [&](cpp_int a) -> bool
    {
        if (a % n == 0)
            return true;
        cpp_int x = mod_pow(a, d, n);
        if (x == 1 || x == n - 1)
            return true;

        for (unsigned r = 1; r < s; ++r)
        {
            x = (x * x) % n;
            if (x == n - 1)
                return true;
        }
        return false;
    };

    // Fixed small bases first
    for (int p : small_primes)
    {
        if (!witness(p))
            return false;
    }

    // Extra random rounds for larger numbers
    if (n > 40)
    {
        std::uniform_int_distribution<uint64_t> dist;
        for (int i = 0; i < rounds; ++i)
        {
            cpp_int a = 2 + cpp_int(dist(rng)) % (n - 3);
            if (!witness(a))
                return false;
        }
    }

    return true;
}

cpp_int pollard_rho(const cpp_int &n)
{
    if ((n & 1) == 0)
        return 2;
    if (n % 3 == 0)
        return 3;

    std::uniform_int_distribution<uint64_t> dist;

    while (true)
    {
        cpp_int c = 1 + cpp_int(dist(rng)) % (n - 1);
        cpp_int x = 2 + cpp_int(dist(rng)) % (n - 3);
        cpp_int y = x;
        cpp_int d = 1;

        auto f = [&](const cpp_int &v)
        {
            return (v * v + c) % n;
        };

        while (d == 1)
        {
            x = f(x);
            y = f(f(y));
            cpp_int diff = x > y ? x - y : y - x;
            d = gcd_cpp(diff, n);
        }

        if (d != n)
            return d;
    }
}

void factor_rec(const cpp_int &n, std::vector<cpp_int> &fac)
{
    if (n == 1)
        return;

    if (is_probable_prime(n))
    {
        fac.push_back(n);
        return;
    }

    cpp_int d = pollard_rho(n);
    factor_rec(d, fac);
    factor_rec(n / d, fac);
}

// Returns factors as (prime, exponent).
// For negative numbers, -1 is included first with exponent 1.
std::vector<std::pair<cpp_int, int>> factorize(cpp_int n)
{
    std::vector<std::pair<cpp_int, int>> res;

    if (n == 0)
    {
        throw std::invalid_argument("0 has no prime factorization");
    }

    if (n < 0)
    {
        res.push_back({-1, 1});
        n = -n;
    }

    if (n == 1)
    {
        return res;
    }

    std::vector<cpp_int> fac;
    factor_rec(n, fac);
    std::sort(fac.begin(), fac.end());

    for (const auto &p : fac)
    {
        if (!res.empty() && res.back().first == p)
        {
            ++res.back().second;
        }
        else
        {
            res.push_back({p, 1});
        }
    }

    return res;
}

void printFactorize(cpp_int n)
{

    auto f = factorize(n);

    for (size_t i = 0; i < f.size(); ++i)
    {
        if (i)
            std::cout << " . ";
        std::cout << f[i].first;
        if (f[i].second > 1)
            std::cout << "^" << f[i].second;
    }
    std::cout << '\n';
}

int main10()
{
// 68406/39919 + 398034*sqrt(5)/199595

    printFactorize(cpp_int("68406"));

    printFactorize(cpp_int("39919"));

    printFactorize(cpp_int("398034"));

    printFactorize(cpp_int("199595"));

    return 0;
}