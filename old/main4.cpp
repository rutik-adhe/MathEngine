// symbolic_complex_exact.cpp
// Exact symbolic algebra with rational and complex simplification.
// Features:
// - Exact Rational arithmetic
// - Exact complex literals a + b*i
// - Addition / multiplication / division / integer powers
// - Basic Euler-form handling via exp/log
// - sqrt / nth-root handling with exact symbolic output
// - No floating-point approximation in simplification
//
// Notes:
// - This is a symbolic engine, not a CAS that proves every algebraic identity.
// - When a form cannot be simplified exactly, it stays symbolic.
//
// Build:
//   g++ -std=c++17 -O2 -Wall -Wextra symbolic_complex_exact.cpp -o symbolic_complex_exact

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using std::make_shared;
using std::shared_ptr;
using std::string;
using std::vector;

// ============================================================
// Exact Rational
// ============================================================

class Rational
{
public:
    long long n; // numerator
    long long d; // denominator, always positive

    Rational(long long num = 0, long long den = 1) : n(num), d(den) { normalize(); }

    void normalize()
    {
        if (d == 0)
            throw std::runtime_error("Rational denominator cannot be zero");
        if (d < 0)
        {
            n = -n;
            d = -d;
        }
        long long g = std::gcd(std::llabs(n), std::llabs(d));
        if (g == 0)
        {
            d = 1;
            return;
        }
        n /= g;
        d /= g;
    }

    bool isZero() const { return n == 0; }
    bool isOne() const { return n == 1 && d == 1; }
    // isMi
    bool isMinusOne() const { return n == -1 && d == 1; }
    bool isInteger() const { return d == 1; }
    bool isPositive() const { return n > 0; }
    bool isNegative() const { return n < 0; }

    string str() const
    {
        if (d == 1)
            return std::to_string(n);
        return std::to_string(n) + "/" + std::to_string(d);
    }

    friend bool operator==(const Rational &a, const Rational &b)
    {
        return a.n == b.n && a.d == b.d;
    }
    friend bool operator!=(const Rational &a, const Rational &b) { return !(a == b); }

    friend Rational operator+(const Rational &a, const Rational &b)
    {
        return Rational(a.n * b.d + b.n * a.d, a.d * b.d);
    }
    friend Rational operator-(const Rational &a, const Rational &b)
    {
        return Rational(a.n * b.d - b.n * a.d, a.d * b.d);
    }
    friend Rational operator*(const Rational &a, const Rational &b)
    {
        return Rational(a.n * b.n, a.d * b.d);
    }
    friend Rational operator/(const Rational &a, const Rational &b)
    {
        if (b.n == 0)
            throw std::runtime_error("Division by zero rational");
        return Rational(a.n * b.d, a.d * b.n);
    }

    Rational operator-() const { return Rational(-n, d); }
};

static Rational powRationalInt(const Rational &base, long long exp)
{
    if (exp == 0)
    {
        if (base.isZero())
            throw std::runtime_error("0 cannot be raised to a 0 power");
        else
            return Rational(1);
    }
    if (base.isZero())
    {
        if (exp < 0)
            throw std::runtime_error("0 cannot be raised to a negative power");
        return Rational(0);
    }

    Rational result(1);
    Rational b = base;
    long long e = std::llabs(exp);
    while (e > 0)
    {
        if (e & 1LL)
            result = result * b;
        b = b * b;
        e >>= 1LL;
    }

    if (exp < 0)
    {
        if (result.isZero())
            throw std::runtime_error("division by zero in negative power");
        return Rational(result.d, result.n);
    }
    return result;
}

static bool exactIntegerNthRoot(long long x, int n, long long &root)
{
    if (n <= 0)
        return false;
    if (x == 0)
    {
        root = 0;
        return true;
    }

    bool neg = x < 0;
    long long ax = std::llabs(x);
    if (neg && (n % 2 == 0))
        return false;

    long long lo = 0, hi = ax;
    while (lo <= hi)
    {
        long long mid = lo + (hi - lo) / 2;
        __int128 p = 1;
        for (int i = 0; i < n; ++i)
        {
            p *= mid;
            if (p > ax)
                break;
        }
        if (p == ax)
        {
            root = neg ? -mid : mid;
            return true;
        }
        if (p < ax)
            lo = mid + 1;
        else
            hi = mid - 1;
    }
    return false;
}

static std::optional<Rational> exactRationalNthRoot(const Rational &x, int n)
{
    if (n <= 0)
        return std::nullopt;
    if (x.isZero())
        return Rational(0);

    bool neg = x.isNegative();
    if (neg && (n % 2 == 0))
        return std::nullopt;

    long long rn = 0, rd = 0;
    if (!exactIntegerNthRoot(std::llabs(x.n), n, rn))
        return std::nullopt;
    if (!exactIntegerNthRoot(x.d, n, rd))
        return std::nullopt;
    if (neg)
        rn = -rn;
    return Rational(rn, rd);
}

// ============================================================
// Complex exact (Gaussian rational)
// ============================================================

struct ComplexExact
{
    Rational re;
    Rational im;

    ComplexExact(Rational r = Rational(0), Rational i = Rational(0)) : re(r), im(i) {}

    bool isZero() const { return re.isZero() && im.isZero(); }
    bool isReal() const { return im.isZero(); }
    bool isPureImag() const { return re.isZero(); }

    string str() const
    {
        if (im.isZero())
            return re.str();
        if (re.isZero())
        {
            if (im.isOne())
                return "i";
            if (im.isMinusOne())
                return "-i";
            return im.str() + "*i";
        }

        string sign = im.isNegative() ? "-" : "+";
        Rational absIm = im.isNegative() ? -im : im;
        string imag = absIm.isOne() ? "i" : absIm.str() + "*i";
        return re.str() + sign + imag;
    }

    friend ComplexExact operator+(const ComplexExact &a, const ComplexExact &b)
    {
        return {a.re + b.re, a.im + b.im};
    }
    friend ComplexExact operator-(const ComplexExact &a, const ComplexExact &b)
    {
        return {a.re - b.re, a.im - b.im};
    }
    friend ComplexExact operator*(const ComplexExact &a, const ComplexExact &b)
    {
        Rational re = a.re * b.re - a.im * b.im;
        Rational im = a.re * b.im + a.im * b.re;
        return {re, im};
    }
    friend ComplexExact operator/(const ComplexExact &a, const ComplexExact &b)
    {
        Rational den = b.re * b.re + b.im * b.im;
        if (den.isZero())
            throw std::runtime_error("division by zero complex");
        ComplexExact num = a * ComplexExact(b.re, -b.im);
        return {num.re / den, num.im / den};
    }
};

static ComplexExact powComplexInt(ComplexExact base, long long exp)
{
    if (exp == 0)
        if (base.isZero())
            throw std::runtime_error("(0 complex) cannot be raised to a (0 complex) power");
        else
            return ComplexExact(1, 0);
    if (base.isZero())
    {
        if (exp < 0)
            throw std::runtime_error("(0 complex) cannot be raised to a negative power");
        return ComplexExact(0, 0);
    }
    bool neg = exp < 0;
    long long e = std::llabs(exp);
    ComplexExact result(1, 0);
    while (e > 0)
    {
        if (e & 1LL)
            result = result * base;
        base = base * base;
        e >>= 1LL;
    }
    if (neg)
        return ComplexExact(1, 0) / result;
    return result;
}

static std::string rationalAbsStr(const Rational &r)
{
    return (r.isNegative() ? (-r).str() : r.str());
}

// ============================================================
// Expression hierarchy
// ============================================================

class Expr
{
public:
    enum class Kind
    {
        Integer,
        Rational,
        Symbol,
        Constant,
        Complex,
        Add,
        Multiply,
        Power,
        Function
    };

    explicit Expr(Kind k) : kind(k) {}
    virtual ~Expr() = default;

    Kind getKind() const { return kind; }

    virtual shared_ptr<Expr> clone() const = 0;
    virtual string toString() const = 0;
    virtual void printPreorder(int indent = 0) const = 0;

private:
    Kind kind;
};

using ExprPtr = shared_ptr<Expr>;

class IntegerNode : public Expr
{
public:
    long long value;
    explicit IntegerNode(long long v) : Expr(Kind::Integer), value(v) {}
    shared_ptr<Expr> clone() const override { return make_shared<IntegerNode>(value); }
    string toString() const override { return std::to_string(value); }
    void printPreorder(int indent = 0) const override
    {
        std::cout << string(indent, ' ') << "Integer(" << value << ")\n";
    }
};

class RationalNode : public Expr
{
public:
    Rational value;
    explicit RationalNode(const Rational &r) : Expr(Kind::Rational), value(r) {}
    shared_ptr<Expr> clone() const override { return make_shared<RationalNode>(value); }
    string toString() const override { return value.str(); }
    void printPreorder(int indent = 0) const override
    {
        std::cout << string(indent, ' ') << "Rational(" << value.str() << ")\n";
    }
};

class SymbolNode : public Expr
{
public:
    string name;
    explicit SymbolNode(string n) : Expr(Kind::Symbol), name(std::move(n)) {}
    shared_ptr<Expr> clone() const override { return make_shared<SymbolNode>(name); }
    string toString() const override { return name; }
    void printPreorder(int indent = 0) const override
    {
        std::cout << string(indent, ' ') << "Symbol(" << name << ")\n";
    }
};

class ConstantNode : public Expr
{
public:
    string name;
    explicit ConstantNode(string n) : Expr(Kind::Constant), name(std::move(n)) {}
    shared_ptr<Expr> clone() const override { return make_shared<ConstantNode>(name); }
    string toString() const override { return name; }
    void printPreorder(int indent = 0) const override
    {
        std::cout << string(indent, ' ') << "Constant(" << name << ")\n";
    }
};

class ComplexNode : public Expr
{
public:
    ComplexExact value;
    explicit ComplexNode(const ComplexExact &v) : Expr(Kind::Complex), value(v) {}
    shared_ptr<Expr> clone() const override { return make_shared<ComplexNode>(value); }
    string toString() const override { return value.str(); }
    void printPreorder(int indent = 0) const override
    {
        std::cout << string(indent, ' ') << "Complex(" << value.str() << ")\n";
    }
};

class AddNode : public Expr
{
public:
    vector<ExprPtr> children;
    explicit AddNode(vector<ExprPtr> c) : Expr(Kind::Add), children(std::move(c)) {}
    shared_ptr<Expr> clone() const override
    {
        vector<ExprPtr> out;
        out.reserve(children.size());
        for (const auto &c : children)
            out.push_back(c->clone());
        return make_shared<AddNode>(out);
    }
    string toString() const override
    {
        std::ostringstream oss;
        oss << "Add(";
        for (size_t i = 0; i < children.size(); ++i)
        {
            if (i)
                oss << ",";
            oss << children[i]->toString();
        }
        oss << ")";
        return oss.str();
    }
    void printPreorder(int indent = 0) const override
    {
        std::cout << string(indent, ' ') << "Add\n";
        for (const auto &c : children)
            c->printPreorder(indent + 2);
    }
};

class MultiplyNode : public Expr
{
public:
    vector<ExprPtr> children;
    explicit MultiplyNode(vector<ExprPtr> c) : Expr(Kind::Multiply), children(std::move(c)) {}
    shared_ptr<Expr> clone() const override
    {
        vector<ExprPtr> out;
        out.reserve(children.size());
        for (const auto &c : children)
            out.push_back(c->clone());
        return make_shared<MultiplyNode>(out);
    }
    string toString() const override
    {
        std::ostringstream oss;
        oss << "Mul(";
        for (size_t i = 0; i < children.size(); ++i)
        {
            if (i)
                oss << ",";
            oss << children[i]->toString();
        }
        oss << ")";
        return oss.str();
    }
    void printPreorder(int indent = 0) const override
    {
        std::cout << string(indent, ' ') << "Multiply\n";
        for (const auto &c : children)
            c->printPreorder(indent + 2);
    }
};

class PowerNode : public Expr
{
public:
    ExprPtr base;
    ExprPtr exponent;
    PowerNode(ExprPtr b, ExprPtr e) : Expr(Kind::Power), base(std::move(b)), exponent(std::move(e)) {}
    shared_ptr<Expr> clone() const override { return make_shared<PowerNode>(base->clone(), exponent->clone()); }
    string toString() const override { return "Pow(" + base->toString() + "," + exponent->toString() + ")"; }
    void printPreorder(int indent = 0) const override
    {
        std::cout << string(indent, ' ') << "Power\n";
        base->printPreorder(indent + 2);
        exponent->printPreorder(indent + 2);
    }
};

class FunctionNode : public Expr
{
public:
    string name;
    vector<ExprPtr> args;
    FunctionNode(string n, vector<ExprPtr> a) : Expr(Kind::Function), name(std::move(n)), args(std::move(a)) {}
    shared_ptr<Expr> clone() const override
    {
        vector<ExprPtr> out;
        out.reserve(args.size());
        for (const auto &a : args)
            out.push_back(a->clone());
        return make_shared<FunctionNode>(name, out);
    }
    string toString() const override
    {
        std::ostringstream oss;
        oss << name << "(";
        for (size_t i = 0; i < args.size(); ++i)
        {
            if (i)
                oss << ",";
            oss << args[i]->toString();
        }
        oss << ")";
        return oss.str();
    }
    void printPreorder(int indent = 0) const override
    {
        std::cout << string(indent, ' ') << "Function(" << name << ")\n";
        for (const auto &a : args)
            a->printPreorder(indent + 2);
    }
};

// ============================================================
// Factory
// ============================================================

class ExprFactory
{
public:
    static ExprPtr integer(long long v) { return make_shared<IntegerNode>(v); }
    static ExprPtr rational(long long n, long long d)
    {
        Rational r(n, d);
        if (r.d == 1)
            return integer(r.n);
        return make_shared<RationalNode>(r);
    }
    static ExprPtr rational(const Rational &r)
    {
        if (r.d == 1)
            return integer(r.n);
        return make_shared<RationalNode>(r);
    }
    static ExprPtr complex(const ComplexExact &v)
    {
        if (v.im.isZero())
            return rational(v.re);
        return make_shared<ComplexNode>(v);
    }
    static ExprPtr symbol(const string &name) { return make_shared<SymbolNode>(name); }
    static ExprPtr constant(const string &name) { return make_shared<ConstantNode>(name); }
    static ExprPtr add(const vector<ExprPtr> &children)
    {
        if (children.size() == 1)
            return children[0];
        return make_shared<AddNode>(children);
    }
    static ExprPtr multiply(const vector<ExprPtr> &children)
    {
        if (children.size() == 1)
            return children[0];
        return make_shared<MultiplyNode>(children);
    }
    static ExprPtr power(ExprPtr base, ExprPtr exponent)
    {
        return make_shared<PowerNode>(std::move(base), std::move(exponent));
    }
    static ExprPtr function(const string &name, const vector<ExprPtr> &args)
    {
        return make_shared<FunctionNode>(name, args);
    }
};

// ============================================================
// Simplifier
// ============================================================

extern vector<std::pair<long long, int>> factorize(long long x);

class Simplifier
{
public:
    static ExprPtr simplify(const ExprPtr &e)
    {
        if (!e)
            return nullptr;
        switch (e->getKind())
        {
        case Expr::Kind::Integer:
        case Expr::Kind::Rational:
        case Expr::Kind::Symbol:
        case Expr::Kind::Constant:
        case Expr::Kind::Complex:
            return e;
        case Expr::Kind::Add:
            return simplifyAdd(std::static_pointer_cast<AddNode>(e));
        case Expr::Kind::Multiply:
            return simplifyMultiply(std::static_pointer_cast<MultiplyNode>(e));
        case Expr::Kind::Power:
            return simplifyPower(std::static_pointer_cast<PowerNode>(e));
        case Expr::Kind::Function:
            return simplifyFunction(std::static_pointer_cast<FunctionNode>(e));
        }
        return e;
    }

private:
    static ExprPtr evalulateArctan(const ExprPtr &e)
    {
        if (!e)
        {
            return e;
        }
    }

    static bool isNumericScalar(const ExprPtr &e, Rational &out)
    {
        if (!e)
            return false;
        if (e->getKind() == Expr::Kind::Integer)
        {
            out = Rational(std::static_pointer_cast<IntegerNode>(e)->value);
            return true;
        }
        if (e->getKind() == Expr::Kind::Rational)
        {
            out = std::static_pointer_cast<RationalNode>(e)->value;
            return true;
        }
        return false;
    }

    static bool isExactComplexLiteral(const ExprPtr &e, ComplexExact &out)
    {
        if (!e)
            return false;
        if (e->getKind() == Expr::Kind::Complex)
        {
            out = std::static_pointer_cast<ComplexNode>(e)->value;
            return true;
        }
        Rational r;
        if (isNumericScalar(e, r))
        {
            out = ComplexExact(r, 0);
            return true;
        }
        if (e->getKind() == Expr::Kind::Constant)
        {
            auto c = std::static_pointer_cast<ConstantNode>(e);
            if (c->name == "i")
            {
                out = ComplexExact(0, 1);
                return true;
            }
            if (c->name == "pi")
            {
                // symbolic constant; not numeric
                return false;
            }
            if (c->name == "e")
            {
                return false;
            }
        }
        return false;
    }

    static ExprPtr toExactExpr(const ComplexExact &z)
    {
        return ExprFactory::complex(z);
    }

    static string keyOf(const ExprPtr &e)
    {
        return e ? e->toString() : "null";
    }

    static void sortCanonical(vector<ExprPtr> &items)
    {
        std::sort(items.begin(), items.end(),
                  [](const ExprPtr &a, const ExprPtr &b)
                  {
                      return a->toString() < b->toString();
                  });
    }

    static bool isOne(const ExprPtr &e)
    {
        Rational r;
        return isNumericScalar(e, r) && r.isOne();
    }

    static bool isZero(const ExprPtr &e)
    {
        Rational r;
        return isNumericScalar(e, r) && r.isZero();
    }

    static std::pair<Rational, ExprPtr> splitCoefficient(const ExprPtr &term)
    {
        Rational r;
        if (isNumericScalar(term, r))
            return {r, nullptr};

        // c * x * y ... => coefficient c and the non-numeric product
        if (term->getKind() == Expr::Kind::Multiply)
        {
            auto m = std::static_pointer_cast<MultiplyNode>(term);
            Rational coeff(1);
            vector<ExprPtr> rest;
            for (const auto &f : m->children)
            {
                Rational fr;
                if (isNumericScalar(f, fr))
                {
                    coeff = coeff * fr;
                }
                else
                {
                    rest.push_back(f);
                }
            }
            if (rest.empty())
                return {coeff, nullptr};
            if (rest.size() == 1)
                return {coeff, rest[0]};
            return {coeff, ExprFactory::multiply(rest)};
        }

        return {Rational(1), term};
    }

    static std::optional<ComplexExact> evalExactComplex(const ExprPtr &e)
    {
        if (!e)
            return std::nullopt;

        ComplexExact z;
        if (isExactComplexLiteral(e, z))
            return z;

        switch (e->getKind())
        {
        case Expr::Kind::Integer:
            // already covered case in isExactComplexLiteral()
        case Expr::Kind::Rational:
            // already covered case in isExactComplexLiteral()
        case Expr::Kind::Complex:
            // already covered case in isExactComplexLiteral()
        case Expr::Kind::Constant:
            // already covered case in isExactComplexLiteral()
        case Expr::Kind::Symbol:
            // return null - Symbols are NON CONVERTABLE to numbers
            return std::nullopt;

        case Expr::Kind::Add:
        {
            auto a = std::static_pointer_cast<AddNode>(e);
            ComplexExact sum(0, 0);
            for (const auto &c : a->children)
            {
                auto v = evalExactComplex(c);
                if (!v)
                    return std::nullopt;
                sum = sum + *v;
            }
            return sum;
        }

        case Expr::Kind::Multiply:
        {
            auto m = std::static_pointer_cast<MultiplyNode>(e);
            ComplexExact prod(1, 0);
            for (const auto &c : m->children)
            {
                auto v = evalExactComplex(c);
                if (!v)
                    return std::nullopt;
                prod = prod * *v;
            }
            return prod;
        }

        case Expr::Kind::Power:
        {
            auto p = std::static_pointer_cast<PowerNode>(e);
            auto b = evalExactComplex(p->base);
            auto ex = simplify(p->exponent);
            if (!b)
                return std::nullopt;

            Rational er;
            if (isNumericScalar(ex, er) && er.isInteger())
            {
                return powComplexInt(*b, er.n);
            }
            return std::nullopt;
        }

        case Expr::Kind::Function:
            return std::nullopt;
        }
        return std::nullopt;
    }

    static ExprPtr simplifyFunction(const std::shared_ptr<FunctionNode> &f)
    {
        vector<ExprPtr> args;
        args.reserve(f->args.size());
        for (const auto &a : f->args)
            args.push_back(simplify(a));

        const string &n = f->name;

        // Basic function normalization
        if (n == "sqrt" && args.size() == 1)
        {
            return simplify(ExprFactory::power(args[0], ExprFactory::rational(1, 2)));
        }
        if (n == "cbrt" && args.size() == 1)
        {
            return simplify(ExprFactory::power(args[0], ExprFactory::rational(1, 3)));
        }
        if (n == "root" && args.size() == 2)
        {
            // root(x, n) = x^(1/n) when n is numeric
            Rational nrt;
            if (isNumericScalar(args[1], nrt) && nrt.isInteger() && nrt.n != 0)
            {
                return simplify(ExprFactory::power(args[0], ExprFactory::rational(1, nrt.n)));
            }
        }

        // exp(log(x)) -> x
        if (n == "exp" && args.size() == 1 && args[0]->getKind() == Expr::Kind::Function)
        {
            auto inner = std::static_pointer_cast<FunctionNode>(args[0]);
            if (inner->name == "log" && inner->args.size() == 1)
                return inner->args[0];
        }

        // log(exp(x)) -> x
        if (n == "log" && args.size() == 1 && args[0]->getKind() == Expr::Kind::Function)
        {
            auto inner = std::static_pointer_cast<FunctionNode>(args[0]);
            if (inner->name == "exp" && inner->args.size() == 1)
                return inner->args[0];
        }

        // exp(0) -> 1
        if (n == "exp" && args.size() == 1)
        {
            Rational r;
            if (isNumericScalar(args[0], r) && r.isZero())
                return ExprFactory::integer(1);
        }

        // log(1) -> 0
        if (n == "log" && args.size() == 1)
        {
            Rational r;
            if (isNumericScalar(args[0], r) && r.isOne())
                return ExprFactory::integer(0);
        }

        // Exact exp on pure imaginary multiples of pi.
        // exp(i*pi*k) collapses exactly for k in Z, and for k = ±1/2.
        if (n == "exp" && args.size() == 1)
        {
            if (args[0]->getKind() == Expr::Kind::Multiply)
            {
                auto m = std::static_pointer_cast<MultiplyNode>(args[0]);
                bool hasI = false, hasPi = false;
                Rational coeff(1);
                vector<ExprPtr> rest;
                for (const auto &c : m->children)
                {
                    if (c->getKind() == Expr::Kind::Constant &&
                        std::static_pointer_cast<ConstantNode>(c)->name == "i")
                    {
                        hasI = true;
                        continue;
                    }
                    if (c->getKind() == Expr::Kind::Constant &&
                        std::static_pointer_cast<ConstantNode>(c)->name == "pi")
                    {
                        hasPi = true;
                        continue;
                    }
                    Rational r;
                    if (isNumericScalar(c, r))
                        coeff = coeff * r;
                    else
                        rest.push_back(c);
                }

                if (hasI && hasPi && rest.empty())
                {
                    if (coeff.isInteger())
                    {
                        long long k = coeff.n;
                        if ((std::llabs(k) % 2LL) == 0LL)
                            return ExprFactory::integer(1);
                        return ExprFactory::integer(-1);
                    }
                    if (coeff == Rational(1, 2))
                        return ExprFactory::constant("i");
                    if (coeff == Rational(-1, 2))
                        return ExprFactory::multiply({ExprFactory::integer(-1), ExprFactory::constant("i")});
                    if (coeff == Rational(3, 2))
                        return ExprFactory::multiply({ExprFactory::integer(-1), ExprFactory::constant("i")});
                    if (coeff == Rational(-3, 2))
                        return ExprFactory::constant("i");
                }
            }
            else
            {
                Rational zero;
                if (isNumericScalar(args[0], zero) && zero.isZero())
                    return ExprFactory::integer(1);
            }
        }

        // Try exact evaluation if all args are exact complex literals and function is elementary numeric
        if (args.size() == 1)
        {
            auto v = evalExactComplex(args[0]);
            if (v)
            {
                if (n == "conj")
                    return toExactExpr(ComplexExact(v->re, -v->im));
                if (n == "re")
                    return toExactExpr(ComplexExact(v->re, 0));
                if (n == "im")
                    return toExactExpr(ComplexExact(v->im, 0));
                if (n == "neg")
                    return toExactExpr(ComplexExact(-v->re, -v->im));
                if (n == "abs")
                {
                    // exact only when purely real and non-negative perfect square or when imag=0
                    if (v->im.isZero())
                    {
                        Rational ar = v->re.isNegative() ? -v->re : v->re;
                        auto rt = exactRationalNthRoot(ar, 2);
                        if (rt)
                            return toExactExpr(ComplexExact(*rt, 0));
                    }
                }
                if (n == "sqrt")
                {
                    // exact principal sqrt of exact complex literal via symbolic radicals
                    // For real values, collapse when perfect squares.
                    if (v->im.isZero())
                    {
                        if (v->re.isZero())
                            return ExprFactory::integer(0);
                        if (v->re.isNegative())
                        {
                            Rational ar = -v->re;
                            auto rt = exactRationalNthRoot(ar, 2);
                            if (rt)
                                return ExprFactory::multiply({ExprFactory::constant("i"), ExprFactory::rational(*rt)});
                            return ExprFactory::multiply({ExprFactory::constant("i"), ExprFactory::function("sqrt", {ExprFactory::rational(ar)})});
                        }
                        auto rt = exactRationalNthRoot(v->re, 2);
                        if (rt)
                            return ExprFactory::rational(*rt);
                        return ExprFactory::function("sqrt", {ExprFactory::rational(v->re)});
                    }

                    // sqrt(a+bi) exact symbolic formula:
                    // u = sqrt((r+a)/2), v = sign(b)*sqrt((r-a)/2), r = sqrt(a^2+b^2)
                    // output u + i*v
                    Rational a = v->re;
                    Rational b = v->im;
                    ExprPtr aE = ExprFactory::rational(a);
                    ExprPtr bE = ExprFactory::rational(b);
                    ExprPtr a2 = ExprFactory::power(aE, ExprFactory::integer(2));
                    ExprPtr b2 = ExprFactory::power(bE, ExprFactory::integer(2));
                    ExprPtr r = ExprFactory::function("sqrt", {ExprFactory::add({a2, b2})});
                    ExprPtr half = ExprFactory::rational(1, 2);
                    ExprPtr uArg = ExprFactory::multiply({half, ExprFactory::add({r, aE})});
                    ExprPtr vArg = ExprFactory::multiply({half, ExprFactory::add({r, ExprFactory::multiply({ExprFactory::integer(-1), aE})})});
                    ExprPtr u = ExprFactory::function("sqrt", {uArg});
                    ExprPtr v = ExprFactory::function("sqrt", {vArg});
                    ExprPtr sign = b.isNegative() ? ExprFactory::integer(-1) : ExprFactory::integer(1);
                    ExprPtr imag = ExprFactory::multiply({sign, ExprFactory::constant("i"), v});
                    return simplify(ExprFactory::add({u, imag}));
                }
            }
        }

        // Leave as symbolic function after simplification of args
        return make_shared<FunctionNode>(n, args);
    }

    static bool extractRootLike(const ExprPtr &e, ExprPtr &base, Rational &exponentOut)
    {
        if (!e)
            return false;

        // Already canonical power: base^(1/n)
        if (e->getKind() == Expr::Kind::Power)
        {
            auto p = std::static_pointer_cast<PowerNode>(e);
            Rational exp;
            if (isNumericScalar(p->exponent, exp))
            {
                base = p->base;
                exponentOut = exp;
                return true;
            }
        }

        // sqrt(x) => x^(1/2)
        if (e->getKind() == Expr::Kind::Function)
        {
            auto f = std::static_pointer_cast<FunctionNode>(e);

            if (f->name == "sqrt" && f->args.size() == 1)
            {
                base = f->args[0];
                exponentOut = Rational(1, 2);
                return true;
            }

            if (f->name == "cbrt" && f->args.size() == 1)
            {
                base = f->args[0];
                exponentOut = Rational(1, 3);
                return true;
            }

            if (f->name == "root" && f->args.size() == 2)
            {
                Rational n;
                if (isNumericScalar(f->args[1], n) && n.isInteger() && n.n != 0)
                {
                    base = f->args[0];
                    exponentOut = Rational(1, n.n);
                    return true;
                }
            }
        }

        return false;
    }

    static ExprPtr simplifyPower(const std::shared_ptr<PowerNode> &p)
    {
        ExprPtr base = simplify(p->base);
        ExprPtr exp = simplify(p->exponent);

        Rational expNum;
        bool expIsNum = isNumericScalar(exp, expNum);

        // collapse nested powers if both exponents numeric
        if (base->getKind() == Expr::Kind::Power)
        {
            auto inner = std::static_pointer_cast<PowerNode>(base);
            Rational innerExp;
            if (isNumericScalar(simplify(inner->exponent), innerExp) && expIsNum)
            {
                return simplify(ExprFactory::power(inner->base, ExprFactory::rational(innerExp * expNum)));
            }
        }

        // x^1 = x
        if (expIsNum && expNum.isOne())
            return base;

        // x^0 = 1 (except 0^0 left symbolic)
        if (expIsNum && expNum.isZero())
        {
            if (!isZero(base))
                return ExprFactory::integer(1);
            else
                throw std::runtime_error("Cannot raise 0 to power 0");
        }

        // 1^x = 1
        Rational baseNum;
        if (isNumericScalar(base, baseNum))
        {
            if (baseNum.isOne())
                return ExprFactory::integer(1);
            if (baseNum.isZero() && expIsNum && expNum.n > 0)
                return ExprFactory::integer(0);
        }

        // i^n exact cycle
        if (base->getKind() == Expr::Kind::Constant &&
            std::static_pointer_cast<ConstantNode>(base)->name == "i" &&
            expIsNum && expNum.isInteger())
        {
            long long m = expNum.n % 4;
            if (m < 0)
                m += 4;
            switch (m)
            {
            case 0:
                return ExprFactory::integer(1);
            case 1:
                return ExprFactory::constant("i");
            case 2:
                return ExprFactory::integer(-1);
            case 3:
                return ExprFactory::multiply({ExprFactory::integer(-1), ExprFactory::constant("i")});
            }
        }

        // exact numeric/complex power
        if (expIsNum)
        {
            auto c = evalExactComplex(base);
            if (c)
            {
                if (expNum.isInteger())
                    return toExactExpr(powComplexInt(*c, expNum.n));
                if (expNum.n == 1 && expNum.d > 1)
                {
                    int n = static_cast<int>(expNum.d);
                    // exact nth root for exact real rational
                    if (c->im.isZero())
                    {
                        auto rt = exactRationalNthRoot(c->re, n);
                        if (rt)
                            return ExprFactory::rational(*rt);
                        if (c->re.isNegative() && (n % 2 == 0))
                        {
                            // principal even root of negative real: i * root(abs)
                            Rational pos = -c->re;
                            auto rt2 = exactRationalNthRoot(pos, n);
                            if (rt2)
                                return simplify(ExprFactory::multiply({ExprFactory::constant("i"), ExprFactory::rational(*rt2)}));
                            return ExprFactory::power(base, exp); // symbolic
                        }
                    }
                }
                else
                {
                    // general rational power (What if its negative)

                    if (c->im.isZero())
                    {
                        // only for base is rational real
                        long long coeff = 1;
                        bool factorTakenOut = false;

                        // what if base is negative
                        vector<std::pair<long long, int>> primeFactorsN = factorize(c->re.n);
                        vector<std::pair<long long, int>> primeFactorsD = factorize(c->re.d);

                        vector<std::pair<long long, ExprPtr>> newPrimeFactorsN;
                        vector<std::pair<long long, ExprPtr>> newPrimeFactorsD;

                        for (const auto &[prime, cnt] : primeFactorsN)
                        {
                            Rational temp = Rational(cnt) * expNum;
                            if (temp.isInteger())
                            {
                                coeff = coeff * temp.n;
                                factorTakenOut = true;
                            }
                            else if (temp.n > temp.d)
                            {
                                factorTakenOut = true;
                                long long factorOut = temp.n > temp.d;
                            }
                        }
                    }
                }
            }
        }

        // sqrt of perfect squares / exact negative square roots
        if (expIsNum && expNum == Rational(1, 2))
        {
            auto c = evalExactComplex(base);
            if (c)
            {
                if (c->isZero())
                    return ExprFactory::integer(0);
                if (c->im.isZero())
                {
                    if (!c->re.isNegative())
                    {
                        auto rt = exactRationalNthRoot(c->re, 2);
                        if (rt)
                            return ExprFactory::rational(*rt);
                    }
                    else
                    {
                        Rational pos = -c->re;
                        auto rt = exactRationalNthRoot(pos, 2);
                        if (rt)
                            return simplify(ExprFactory::multiply({ExprFactory::constant("i"), ExprFactory::rational(*rt)}));
                    }
                }
            }
        }

        // If exponent is 1/n and base exact rational perfect nth power.
        if (expIsNum && expNum.n == 1 && expNum.d > 1)
        {
            auto c = evalExactComplex(base);
            if (c && c->im.isZero())
            {
                auto rt = exactRationalNthRoot(c->re, static_cast<int>(expNum.d));
                if (rt)
                    return ExprFactory::rational(*rt);
            }
        }

        // (a^m)^n = a^(m*n) already handled above

        return ExprFactory::power(base, exp);
    }

    static ExprPtr simplifyMultiply(const std::shared_ptr<MultiplyNode> &m)
    {
        // remove unneccessary bracketing inside multiply operator ((a).(b).(c)) = (a.b.c)
        vector<ExprPtr> factors;
        for (const auto &c : m->children)
        {
            ExprPtr s = simplify(c);
            if (s->getKind() == Expr::Kind::Multiply)
            {
                auto mm = std::static_pointer_cast<MultiplyNode>(s);
                for (const auto &k : mm->children)
                    factors.push_back(k);
            }
            else
            {
                factors.push_back(s);
            }
        }

        // exact numeric complex fold first
        bool allExact = true;
        ComplexExact prod(1, 0);
        for (const auto &f : factors)
        {
            auto z = evalExactComplex(f);
            if (!z)
            {
                allExact = false;
                break;
            }
            prod = prod * *z;
        }
        if (allExact)
            return toExactExpr(prod);

        Rational coeff(1);
        int iCount = 0;

        // combine base^exponent terms, including root-like forms
        std::map<string, std::pair<ExprPtr, Rational>> powerMap;
        vector<ExprPtr> rest;

        for (const auto &f : factors)
        {
            Rational r;
            if (isNumericScalar(f, r))
            {
                coeff = coeff * r;
                continue;
            }

            if (f->getKind() == Expr::Kind::Constant &&
                std::static_pointer_cast<ConstantNode>(f)->name == "i")
            {
                iCount++;
                continue;
            }

            ExprPtr base;
            Rational exp;
            if (extractRootLike(f, base, exp))
            {
                string key = base->toString();
                auto it = powerMap.find(key);
                if (it == powerMap.end())
                {
                    powerMap[key] = {base, exp};
                }
                else
                {
                    it->second.second = it->second.second + exp;
                }
                continue;
            }

            rest.push_back(f);
        }

        if (coeff.isZero())
            return ExprFactory::integer(0);

        // normalize i^k
        iCount %= 4;
        if (iCount < 0)
            iCount += 4;
        switch (iCount)
        {
        case 0:
            break;
        case 1:
            rest.push_back(ExprFactory::constant("i"));
            break;
        case 2:
            coeff = coeff * Rational(-1);
            break;
        case 3:
            coeff = coeff * Rational(-1);
            rest.push_back(ExprFactory::constant("i"));
            break;
        }

        vector<ExprPtr> out;
        // powerMap.empty() implies there is no pair added
        if (!coeff.isOne() || (powerMap.empty() && rest.empty()))
        {
            out.push_back(ExprFactory::rational(coeff));
        }

        for (auto &kv : powerMap)
        {
            ExprPtr base = kv.second.first;
            Rational exp = kv.second.second;

            if (exp.isZero())
                continue;
            if (exp.isOne())
            {
                out.push_back(base);
            }
            else
            {
                out.push_back(simplify(ExprFactory::power(base, ExprFactory::rational(exp))));
            }
        }

        for (const auto &r : rest)
            out.push_back(r);

        sortCanonical(out);

        if (out.empty())
            return ExprFactory::integer(1);
        if (out.size() == 1)
            return out[0];
        return ExprFactory::multiply(out);
    }

    static ExprPtr simplifyAdd(const std::shared_ptr<AddNode> &a)
    {
        vector<ExprPtr> terms;
        for (const auto &c : a->children)
        {
            ExprPtr s = simplify(c);
            if (s->getKind() == Expr::Kind::Add)
            {
                auto aa = std::static_pointer_cast<AddNode>(s);
                for (const auto &k : aa->children)
                    terms.push_back(k);
            }
            else
            {
                terms.push_back(s);
            }
        }

        // Fold exact complex sum if all terms are exact
        vector<ExprPtr> NonExactTerms;

        bool allExact = true;
        ComplexExact sumTerms(0, 0);
        for (const auto &t : terms)
        {
            auto z = evalExactComplex(t);
            if (!z)
            {
                allExact = false;
                break;
                // NonExactTerms.push_back(t);
            }
            else
            {
                sumTerms = sumTerms + *z;
            }
        }
        if (allExact)
            return toExactExpr(sumTerms);

        // Combine like terms by exact symbolic key, with numeric coefficient extraction.
        // Rational numericSum(0);
        std::map<string, std::pair<Rational, ExprPtr>> groups;

        ComplexExact sum(0, 0);
        for (const auto &t : terms)
        {
            // Rational r;
            // if (isNumericScalar(t, r))
            // {
            //     numericSum = numericSum + r;
            //     continue;
            // }

            auto z = evalExactComplex(t);
            if (z)
            {
                sum = sum + *z;
                continue;
            }

            auto [coeff, rest] = splitCoefficient(t);
            if (!rest)
            {
                // numericSum = numericSum + coeff;
                sum = sum + *evalExactComplex(ExprFactory::rational(coeff));
                continue;
            }

            string key = keyOf(rest);
            auto it = groups.find(key);
            if (it == groups.end())
                groups[key] = {coeff, rest};
            else
                it->second.first = it->second.first + coeff;
        }

        vector<ExprPtr> out;
        if (!sum.isZero())
            if (sum.isReal())
            {
                out.push_back(ExprFactory::rational(sum.re));
            }
            else
            {
                out.push_back(ExprFactory::complex(sum));
            }
        // out.push_back(ExprFactory::rational(numericSum));

        for (auto &kv : groups)
        {
            Rational coeff = kv.second.first;
            ExprPtr rest = kv.second.second;
            if (coeff.isZero())
                continue;
            if (coeff.isOne())
            {
                out.push_back(rest);
            }
            else
            {
                out.push_back(simplify(ExprFactory::multiply({ExprFactory::rational(coeff), rest})));
            }
        }

        sortCanonical(out);
        if (out.empty())
            return ExprFactory::integer(0);
        if (out.size() == 1)
            return out[0];
        return ExprFactory::add(out);
    }
};

int main1()
{
    ExprPtr i = ExprFactory::constant("i");
    ExprPtr x = ExprFactory::symbol("x");
    ExprPtr y = ExprFactory::symbol("y");
    ExprPtr pi = ExprFactory::constant("pi");

    // Example 1: exact complex arithmetic
    ExprPtr complexPart = ExprFactory::power(
        ExprFactory::add({ExprFactory::integer(3),
                          ExprFactory::multiply({ExprFactory::integer(2), i})}),
        ExprFactory::integer(2));

    // Example 2: exact symbolic root
    ExprPtr sqrt9 = ExprFactory::function("sqrt", {ExprFactory::integer(9)});

    // Example 3: symbolic terms remain symbolic
    ExprPtr symbolic = ExprFactory::add({ExprFactory::multiply({ExprFactory::integer(6), x, y}),
                                         ExprFactory::multiply({ExprFactory::integer(2), x, y}),
                                         ExprFactory::add({x, x}),
                                         sqrt9});

    // Example 4: Euler form / exact exp rule
    ExprPtr euler = ExprFactory::function("exp", {ExprFactory::multiply({i, pi})});

    ExprPtr expr = ExprFactory::add({complexPart,
                                     symbolic,
                                     euler});

    std::cout << "Original tree (preorder):\n";
    expr->printPreorder();

    ExprPtr simp = Simplifier::simplify(expr);

    std::cout << "\nSimplified tree (preorder):\n";
    simp->printPreorder();

    std::cout << "\nCanonical string form:\n";
    std::cout << simp->toString() << "\n";

    return 0;
}

int main2()
{
    ExprPtr expr = ExprFactory::power(ExprFactory::add(
                                          {ExprFactory::integer(1),
                                           ExprFactory::power(ExprFactory::integer(2), ExprFactory::rational(1, 2)),
                                           ExprFactory::power(ExprFactory::integer(3), ExprFactory::rational(1, 2))}),
                                      ExprFactory::integer(3));

    std::cout << "Original tree (preorder):\n";
    expr->printPreorder();

    ExprPtr simp = Simplifier::simplify(expr);

    std::cout << "\nSimplified tree (preorder):\n";
    simp->printPreorder();

    return 0;
}

int main3()
{

    ExprPtr expr = ExprFactory::power(ExprFactory::integer(72), ExprFactory::rational(1, 2));

    std::cout << "Original tree (preorder):\n";
    expr->printPreorder();

    ExprPtr simp = Simplifier::simplify(expr);

    std::cout << "\nSimplified tree (preorder):\n";
    simp->printPreorder();

    return 0;
}

int main5()
{
    ExprPtr expr = ExprFactory::power(ExprFactory::complex({1, 2}), ExprFactory::complex({2, 3}));

    std::cout << "Original tree (preorder):\n";
    expr->printPreorder();

    ExprPtr simp = Simplifier::simplify(expr);

    std::cout << "\nSimplified tree (preorder):\n";
    simp->printPreorder();

    return 0;
}

int main6()
{
    ExprPtr expr = ExprFactory::power(ExprFactory::complex({1, 2}), ExprFactory::rational(1, 2));

    std::cout << "Original tree (preorder):\n";
    expr->printPreorder();

    ExprPtr simp = Simplifier::simplify(expr);

    std::cout << "\nSimplified tree (preorder):\n";
    simp->printPreorder();

    return 0;
}

int main()
{
    main3();
    return 0;
}