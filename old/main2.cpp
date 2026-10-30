#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <variant>
#include <vector>

using std::make_shared;
using std::shared_ptr;
using std::string;
using std::vector;

using C = std::complex<long double>;
static constexpr long double EPS = 1e-12L;

// ============================================================
// Helpers
// ============================================================

static bool nearZero(long double x)
{
    return std::fabsl(x) < EPS;
}

static bool nearEqual(long double a, long double b)
{
    return std::fabsl(a - b) < EPS;
}

static bool nearInt(long double x)
{
    return nearEqual(x, std::round(x));
}

static long long toIntExact(long double x)
{
    return static_cast<long long>(std::llround(x));
}

static string trimTrailingZeros(string s)
{
    auto pos = s.find('.');
    if (pos == string::npos)
        return s;

    while (!s.empty() && s.back() == '0')
        s.pop_back();
    if (!s.empty() && s.back() == '.')
        s.pop_back();
    if (s == "-0")
        s = "0";
    return s;
}

static string formatReal(long double x)
{
    if (nearInt(x))
        return std::to_string(toIntExact(x));
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(12) << x;
    return trimTrailingZeros(oss.str());
}

static string formatComplex(const C &z)
{
    long double re = z.real();
    long double im = z.imag();

    if (nearZero(im))
    {
        return formatReal(re);
    }

    if (nearZero(re))
    {
        if (nearEqual(im, 1.0L))
            return "i";
        if (nearEqual(im, -1.0L))
            return "-i";
        return formatReal(im) + "i";
    }

    string sign = (im >= 0 ? "+" : "");
    long double absIm = std::fabsl(im);

    if (nearEqual(absIm, 1.0L))
    {
        return formatReal(re) + sign + (im >= 0 ? "i" : "-i");
    }

    return formatReal(re) + sign + formatReal(im) + "i";
}

// ============================================================
// Exact Rational
// ============================================================

class Rational
{
public:
    long long n; // numerator
    long long d; // denominator, always positive

    Rational(long long num = 0, long long den = 1) : n(num), d(den)
    {
        normalize();
    }

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
    bool isMinusOne() const { return n == -1 && d == 1; }
    bool isInteger() const { return d == 1; }

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

    friend bool operator!=(const Rational &a, const Rational &b)
    {
        return !(a == b);
    }

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
};

static Rational powRationalInt(const Rational &base, long long exp)
{
    if (exp == 0)
        return Rational(1);
    if (base.isZero())
    {
        if (exp < 0)
            throw std::runtime_error("0 cannot be raised to negative power");
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
        return Rational(result.d, result.n);
    }
    return result;
}

// exact integer nth root for perfect powers
static bool exactNthRoot(long long x, int n, long long &root)
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
        else if (p < ax)
        {
            lo = mid + 1;
        }
        else
        {
            hi = mid - 1;
        }
    }

    return false;
}

static bool exactNthRootRational(const Rational &x, int n, Rational &root)
{
    if (n <= 0)
        return false;
    if (x.isZero())
    {
        root = Rational(0);
        return true;
    }

    long long rn = 0, rd = 0;
    bool neg = x.n < 0;

    if (neg && (n % 2 == 0))
        return false;

    if (!exactNthRoot(std::llabs(x.n), n, rn))
        return false;
    if (!exactNthRoot(x.d, n, rd))
        return false;

    if (neg)
        rn = -rn;
    root = Rational(rn, rd);
    return true;
}

// ============================================================
// Numeric value variant
// ============================================================

using Num = std::variant<Rational, C>;

static bool isRationalNum(const Num &n)
{
    return std::holds_alternative<Rational>(n);
}

static Rational asRational(const Num &n)
{
    return std::get<Rational>(n);
}

static C asComplex(const Num &n)
{
    if (std::holds_alternative<Rational>(n))
    {
        const auto &r = std::get<Rational>(n);
        return C(static_cast<long double>(r.n) / static_cast<long double>(r.d), 0.0L);
    }
    return std::get<C>(n);
}

static bool numIsZero(const Num &n)
{
    if (std::holds_alternative<Rational>(n))
        return std::get<Rational>(n).isZero();
    return nearZero(std::get<C>(n).real()) && nearZero(std::get<C>(n).imag());
}

static bool numIsOne(const Num &n)
{
    if (std::holds_alternative<Rational>(n))
        return std::get<Rational>(n).isOne();
    C z = std::get<C>(n);
    return nearEqual(z.real(), 1.0L) && nearZero(z.imag());
}

static bool numIsMinusOne(const Num &n)
{
    if (std::holds_alternative<Rational>(n))
        return std::get<Rational>(n).isMinusOne();
    C z = std::get<C>(n);
    return nearEqual(z.real(), -1.0L) && nearZero(z.imag());
}

static string numToString(const Num &n)
{
    if (std::holds_alternative<Rational>(n))
        return std::get<Rational>(n).str();
    return formatComplex(std::get<C>(n));
}

static Num numAdd(const Num &a, const Num &b)
{
    if (std::holds_alternative<Rational>(a) && std::holds_alternative<Rational>(b))
    {
        return std::get<Rational>(a) + std::get<Rational>(b);
    }
    return asComplex(a) + asComplex(b);
}

static Num numSub(const Num &a, const Num &b)
{
    if (std::holds_alternative<Rational>(a) && std::holds_alternative<Rational>(b))
    {
        return std::get<Rational>(a) - std::get<Rational>(b);
    }
    return asComplex(a) - asComplex(b);
}

static Num numMul(const Num &a, const Num &b)
{
    if (std::holds_alternative<Rational>(a) && std::holds_alternative<Rational>(b))
    {
        return std::get<Rational>(a) * std::get<Rational>(b);
    }
    return asComplex(a) * asComplex(b);
}

static Num numDiv(const Num &a, const Num &b)
{
    if (std::holds_alternative<Rational>(a) && std::holds_alternative<Rational>(b))
    {
        return std::get<Rational>(a) / std::get<Rational>(b);
    }
    return asComplex(a) / asComplex(b);
}

static Num numPow(const Num &base, const Num &exp)
{
    // Exact rational integer exponent on exact rational base
    if (std::holds_alternative<Rational>(base) && std::holds_alternative<Rational>(exp))
    {
        const Rational &b = std::get<Rational>(base);
        const Rational &e = std::get<Rational>(exp);

        if (e.isInteger())
        {
            return powRationalInt(b, e.n);
        }

        // exact root when exponent is +/- 1/n and base is perfect nth power
        if (std::llabs(e.n) == 1 && e.d > 1)
        {
            Rational root;
            if (exactNthRootRational(b, static_cast<int>(e.d), root))
            {
                if (e.n == 1)
                    return root;
                return Rational(1) / root;
            }
        }

        // otherwise principal complex power
        return std::pow(asComplex(base), asComplex(exp));
    }

    return std::pow(asComplex(base), asComplex(exp));
}

static Num numSqrt(const Num &x)
{
    return numPow(x, Rational(1, 2));
}

static std::optional<Num> makeNumericFromComplex(const C &z)
{
    if (nearZero(z.imag()))
    {
        long double re = z.real();
        if (nearInt(re))
        {
            return Rational(toIntExact(re));
        }
        return C(re, 0.0L);
    }
    return z;
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
        Complex,
        Symbol,
        Constant,
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

    shared_ptr<Expr> clone() const override
    {
        return make_shared<IntegerNode>(value);
    }

    string toString() const override
    {
        return std::to_string(value);
    }

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

    shared_ptr<Expr> clone() const override
    {
        return make_shared<RationalNode>(value);
    }

    string toString() const override
    {
        return value.str();
    }

    void printPreorder(int indent = 0) const override
    {
        std::cout << string(indent, ' ') << "Rational(" << value.str() << ")\n";
    }
};

class ComplexNode : public Expr
{
public:
    C value;

    explicit ComplexNode(const C &z) : Expr(Kind::Complex), value(z) {}

    shared_ptr<Expr> clone() const override
    {
        return make_shared<ComplexNode>(value);
    }

    string toString() const override
    {
        return formatComplex(value);
    }

    void printPreorder(int indent = 0) const override
    {
        std::cout << string(indent, ' ') << "Complex(" << formatComplex(value) << ")\n";
    }
};

class SymbolNode : public Expr
{
public:
    string name;

    explicit SymbolNode(string n) : Expr(Kind::Symbol), name(std::move(n)) {}

    shared_ptr<Expr> clone() const override
    {
        return make_shared<SymbolNode>(name);
    }

    string toString() const override
    {
        return name;
    }

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

    shared_ptr<Expr> clone() const override
    {
        return make_shared<ConstantNode>(name);
    }

    string toString() const override
    {
        return name;
    }

    void printPreorder(int indent = 0) const override
    {
        std::cout << string(indent, ' ') << "Constant(" << name << ")\n";
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
        vector<string> parts;
        for (const auto &c : children)
            parts.push_back(c->toString());
        std::ostringstream oss;
        oss << "Add(";
        for (size_t i = 0; i < parts.size(); ++i)
        {
            if (i)
                oss << ",";
            oss << parts[i];
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
        vector<string> parts;
        for (const auto &c : children)
            parts.push_back(c->toString());
        std::ostringstream oss;
        oss << "Mul(";
        for (size_t i = 0; i < parts.size(); ++i)
        {
            if (i)
                oss << ",";
            oss << parts[i];
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

    shared_ptr<Expr> clone() const override
    {
        return make_shared<PowerNode>(base->clone(), exponent->clone());
    }

    string toString() const override
    {
        return "Pow(" + base->toString() + "," + exponent->toString() + ")";
    }

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
        oss << "Fn(" << name << "(";
        for (size_t i = 0; i < args.size(); ++i)
        {
            if (i)
                oss << ",";
            oss << args[i]->toString();
        }
        oss << "))";
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
// Factory helpers
// ============================================================

class ExprFactory
{
public:
    static ExprPtr integer(long long v)
    {
        return make_shared<IntegerNode>(v);
    }

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

    static ExprPtr complex(const C &z)
    {
        if (nearZero(z.imag()))
        {
            long double re = z.real();
            if (nearInt(re))
                return integer(toIntExact(re));
            return make_shared<ComplexNode>(C(re, 0.0L));
        }
        return make_shared<ComplexNode>(z);
    }

    static ExprPtr symbol(const string &name)
    {
        return make_shared<SymbolNode>(name);
    }

    static ExprPtr constant(const string &name)
    {
        return make_shared<ConstantNode>(name);
    }

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
// Canonical string key / numeric conversion
// ============================================================

static string exprKey(const ExprPtr &e)
{
    return e ? e->toString() : "null";
}

static std::optional<Num> evalNumeric(const ExprPtr &e);

static ExprPtr numToExpr(const Num &n)
{
    if (std::holds_alternative<Rational>(n))
    {
        return ExprFactory::rational(std::get<Rational>(n));
    }
    return ExprFactory::complex(std::get<C>(n));
}

// ============================================================
// Numeric evaluation of supported expressions
// ============================================================

static std::optional<Num> evalFunctionNumeric(const std::shared_ptr<FunctionNode> &fn)
{
    const string &name = fn->name;

    if (name == "sqrt" && fn->args.size() == 1)
    {
        auto a = evalNumeric(fn->args[0]);
        if (!a)
            return std::nullopt;
        return numSqrt(*a);
    }

    if (name == "root" && fn->args.size() == 2)
    {
        auto degree = evalNumeric(fn->args[0]);
        auto radicand = evalNumeric(fn->args[1]);
        if (!degree || !radicand)
            return std::nullopt;

        // principal nth root: radicand^(1/degree)
        if (std::holds_alternative<Rational>(*degree) && std::get<Rational>(*degree).isInteger())
        {
            long long n = std::get<Rational>(*degree).n;
            if (n == 0)
                return std::nullopt;
            return numPow(*radicand, Rational(1, n));
        }

        C deg = asComplex(*degree);
        C rad = asComplex(*radicand);
        return std::pow(rad, C(1.0L, 0.0L) / deg);
    }

    auto unaryComplex = [&](auto fnc) -> std::optional<Num>
    {
        if (fn->args.size() != 1)
            return std::nullopt;
        auto a = evalNumeric(fn->args[0]);
        if (!a)
            return std::nullopt;
        return fnc(asComplex(*a));
    };

    if (name == "exp")
        return unaryComplex([](const C &z)
                            { return std::exp(z); });
    if (name == "log")
        return unaryComplex([](const C &z)
                            { return std::log(z); });
    if (name == "sin")
        return unaryComplex([](const C &z)
                            { return std::sin(z); });
    if (name == "cos")
        return unaryComplex([](const C &z)
                            { return std::cos(z); });
    if (name == "tan")
        return unaryComplex([](const C &z)
                            { return std::tan(z); });
    if (name == "asin")
        return unaryComplex([](const C &z)
                            { return std::asin(z); });
    if (name == "acos")
        return unaryComplex([](const C &z)
                            { return std::acos(z); });
    if (name == "atan")
        return unaryComplex([](const C &z)
                            { return std::atan(z); });
    if (name == "conj")
        return unaryComplex([](const C &z)
                            { return std::conj(z); });

    if (name == "abs" && fn->args.size() == 1)
    {
        auto a = evalNumeric(fn->args[0]);
        if (!a)
            return std::nullopt;
        long double v = std::abs(asComplex(*a));
        if (nearInt(v))
            return Rational(toIntExact(v));
        return C(v, 0.0L);
    }

    if (name == "arg" && fn->args.size() == 1)
    {
        auto a = evalNumeric(fn->args[0]);
        if (!a)
            return std::nullopt;
        long double v = std::arg(asComplex(*a));
        if (nearInt(v))
            return Rational(toIntExact(v));
        return C(v, 0.0L);
    }

    if (name == "real" && fn->args.size() == 1)
    {
        auto a = evalNumeric(fn->args[0]);
        if (!a)
            return std::nullopt;
        long double v = asComplex(*a).real();
        if (nearInt(v))
            return Rational(toIntExact(v));
        return C(v, 0.0L);
    }

    if (name == "imag" && fn->args.size() == 1)
    {
        auto a = evalNumeric(fn->args[0]);
        if (!a)
            return std::nullopt;
        long double v = asComplex(*a).imag();
        if (nearInt(v))
            return Rational(toIntExact(v));
        return C(v, 0.0L);
    }

    if (name == "cis" && fn->args.size() == 1)
    {
        auto a = evalNumeric(fn->args[0]);
        if (!a)
            return std::nullopt;
        C t = asComplex(*a);
        return std::polar(1.0L, t.real());
    }

    return std::nullopt;
}

static std::optional<Num> evalNumeric(const ExprPtr &e)
{
    if (!e)
        return std::nullopt;

    switch (e->getKind())
    {
    case Expr::Kind::Integer:
    {
        auto n = std::static_pointer_cast<IntegerNode>(e);
        return Rational(n->value);
    }

    case Expr::Kind::Rational:
    {
        auto n = std::static_pointer_cast<RationalNode>(e);
        return n->value;
    }

    case Expr::Kind::Complex:
    {
        auto n = std::static_pointer_cast<ComplexNode>(e);
        return n->value;
    }

    case Expr::Kind::Constant:
    {
        auto n = std::static_pointer_cast<ConstantNode>(e);
        if (n->name == "i")
            return C(0.0L, 1.0L);
        if (n->name == "pi")
            return C(std::acos(-1.0L), 0.0L);
        if (n->name == "e")
            return C(std::expl(1.0L), 0.0L);
        if (n->name == "tau")
            return C(2.0L * std::acos(-1.0L), 0.0L);
        return std::nullopt;
    }

    case Expr::Kind::Symbol:
        return std::nullopt;

    case Expr::Kind::Add:
    {
        auto a = std::static_pointer_cast<AddNode>(e);
        Num acc = Rational(0);
        for (const auto &c : a->children)
        {
            auto v = evalNumeric(c);
            if (!v)
                return std::nullopt;
            acc = numAdd(acc, *v);
        }
        return acc;
    }

    case Expr::Kind::Multiply:
    {
        auto m = std::static_pointer_cast<MultiplyNode>(e);
        Num acc = Rational(1);
        for (const auto &c : m->children)
        {
            auto v = evalNumeric(c);
            if (!v)
                return std::nullopt;
            acc = numMul(acc, *v);
        }
        return acc;
    }

    case Expr::Kind::Power:
    {
        auto p = std::static_pointer_cast<PowerNode>(e);
        auto b = evalNumeric(p->base);
        auto ex = evalNumeric(p->exponent);
        if (!b || !ex)
            return std::nullopt;
        return numPow(*b, *ex);
    }

    case Expr::Kind::Function:
    {
        auto fn = std::static_pointer_cast<FunctionNode>(e);
        return evalFunctionNumeric(fn);
    }
    }

    return std::nullopt;
}

// ============================================================
// Simplifier
// ============================================================

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
        case Expr::Kind::Complex:
        case Expr::Kind::Symbol:
        case Expr::Kind::Constant:
            return e;

        case Expr::Kind::Function:
            return simplifyFunction(std::static_pointer_cast<FunctionNode>(e));

        case Expr::Kind::Power:
            return simplifyPower(std::static_pointer_cast<PowerNode>(e));

        case Expr::Kind::Multiply:
            return simplifyMultiply(std::static_pointer_cast<MultiplyNode>(e));

        case Expr::Kind::Add:
            return simplifyAdd(std::static_pointer_cast<AddNode>(e));
        }

        return e;
    }

private:
    static void sortCanonical(vector<ExprPtr> &items)
    {
        std::sort(items.begin(), items.end(),
                  [](const ExprPtr &a, const ExprPtr &b)
                  {
                      return exprKey(a) < exprKey(b);
                  });
    }

    static std::pair<Num, ExprPtr> splitCoefficient(const ExprPtr &term)
    {
        auto nv = evalNumeric(term);
        if (nv)
        {
            return {*nv, nullptr};
        }

        if (term->getKind() == Expr::Kind::Multiply)
        {
            auto m = std::static_pointer_cast<MultiplyNode>(term);
            Num coeff = Rational(1);
            vector<ExprPtr> rest;

            for (const auto &f : m->children)
            {
                auto v = evalNumeric(f);
                if (v)
                {
                    coeff = numMul(coeff, *v);
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

    static ExprPtr simplifyFunction(const std::shared_ptr<FunctionNode> &f)
    {
        vector<ExprPtr> args;
        args.reserve(f->args.size());
        for (const auto &a : f->args)
            args.push_back(simplify(a));

        auto rebuilt = ExprFactory::function(f->name, args);

        if (auto nv = evalNumeric(rebuilt))
        {
            return numToExpr(*nv);
        }

        return rebuilt;
    }

    static ExprPtr simplifyPower(const std::shared_ptr<PowerNode> &p)
    {
        ExprPtr base = simplify(p->base);
        ExprPtr exp = simplify(p->exponent);

        auto full = evalNumeric(ExprFactory::power(base, exp));
        if (full)
        {
            return numToExpr(*full);
        }

        auto expNv = evalNumeric(exp);

        // x^1 = x
        if (expNv && numIsOne(*expNv))
            return base;

        // x^0 = 1  (except 0^0, left as-is)
        if (expNv && numIsZero(*expNv))
        {
            auto baseNv = evalNumeric(base);
            if (!(baseNv && numIsZero(*baseNv)))
            {
                return ExprFactory::integer(1);
            }
        }

        // 1^x = 1
        auto baseNv = evalNumeric(base);
        if (baseNv && numIsOne(*baseNv))
            return ExprFactory::integer(1);

        // 0^positive = 0
        if (baseNv && numIsZero(*baseNv) && expNv)
        {
            if (std::holds_alternative<Rational>(*expNv))
            {
                const Rational &r = std::get<Rational>(*expNv);
                if (r.n > 0)
                    return ExprFactory::integer(0);
            }
            else
            {
                C z = std::get<C>(*expNv);
                if (nearZero(z.imag()) && z.real() > 0)
                    return ExprFactory::integer(0);
            }
        }

        // (a^m)^n -> a^(m*n) when exponents are numeric
        if (base->getKind() == Expr::Kind::Power)
        {
            auto inner = std::static_pointer_cast<PowerNode>(base);
            auto innerExp = evalNumeric(inner->exponent);
            if (innerExp && expNv)
            {
                return simplify(ExprFactory::power(inner->base, numToExpr(numMul(*innerExp, *expNv))));
            }
        }

        return ExprFactory::power(base, exp);
    }

    static ExprPtr simplifyMultiply(const std::shared_ptr<MultiplyNode> &m)
    {
        // 1) simplify children and flatten nested multiplies
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

        // 2) collect coefficient and merge repeated symbolic bases
        Num coeff = Rational(1);

        // key -> (base, exponent)
        std::map<string, std::pair<ExprPtr, Rational>> powMap;
        vector<ExprPtr> otherFactors;

        for (const auto &f : factors)
        {
            auto v = evalNumeric(f);
            if (v)
            {
                coeff = numMul(coeff, *v);
                continue;
            }

            if (f->getKind() == Expr::Kind::Power)
            {
                auto p = std::static_pointer_cast<PowerNode>(f);
                auto expNv = evalNumeric(p->exponent);

                if (expNv && std::holds_alternative<Rational>(*expNv))
                {
                    string key = exprKey(p->base);
                    Rational ex = std::get<Rational>(*expNv);
                    auto it = powMap.find(key);
                    if (it == powMap.end())
                    {
                        powMap[key] = {p->base, ex};
                    }
                    else
                    {
                        it->second.second = it->second.second + ex;
                    }
                    continue;
                }
            }

            // treat factor as exponent 1
            string key = exprKey(f);
            auto it = powMap.find(key);
            if (it == powMap.end())
            {
                powMap[key] = {f, Rational(1)};
            }
            else
            {
                it->second.second = it->second.second + Rational(1);
            }
        }

        // zero annihilates everything
        if (numIsZero(coeff))
        {
            return ExprFactory::integer(0);
        }

        vector<ExprPtr> out;

        // rebuild repeated bases
        for (auto &kv : powMap)
        {
            ExprPtr base = kv.second.first;
            Rational ex = kv.second.second;

            if (ex.isZero())
                continue;
            if (ex.isOne())
            {
                out.push_back(base);
            }
            else
            {
                out.push_back(simplify(ExprFactory::power(base, ExprFactory::rational(ex))));
            }
        }

        // attach coefficient
        if (!numIsOne(coeff) || out.empty())
        {
            if (out.empty())
            {
                return numToExpr(coeff);
            }
            out.insert(out.begin(), numToExpr(coeff));
        }

        sortCanonical(out);

        if (out.size() == 1)
            return out[0];
        return ExprFactory::multiply(out);
    }

    static ExprPtr simplifyAdd(const std::shared_ptr<AddNode> &a)
    {
        // 1) simplify children and flatten nested adds
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

        // 2) combine numeric terms and like terms
        Num numericSum = Rational(0);
        bool haveNumeric = false;

        std::map<string, std::pair<Num, ExprPtr>> groups;

        for (const auto &t : terms)
        {
            auto nv = evalNumeric(t);
            if (nv)
            {
                if (!haveNumeric)
                {
                    numericSum = *nv;
                    haveNumeric = true;
                }
                else
                {
                    numericSum = numAdd(numericSum, *nv);
                }
                continue;
            }

            auto [coeff, rest] = splitCoefficient(t);

            if (!rest)
            {
                if (!haveNumeric)
                {
                    numericSum = coeff;
                    haveNumeric = true;
                }
                else
                {
                    numericSum = numAdd(numericSum, coeff);
                }
                continue;
            }

            string key = exprKey(rest);
            auto it = groups.find(key);
            if (it == groups.end())
            {
                groups[key] = {coeff, rest};
            }
            else
            {
                it->second.first = numAdd(it->second.first, coeff);
            }
        }

        // 3) rebuild
        vector<ExprPtr> out;

        if (haveNumeric && !numIsZero(numericSum))
        {
            out.push_back(numToExpr(numericSum));
        }

        for (auto &kv : groups)
        {
            Num coeff = kv.second.first;
            ExprPtr rest = kv.second.second;

            if (numIsZero(coeff))
                continue;

            if (numIsOne(coeff))
            {
                out.push_back(rest);
            }
            else
            {
                out.push_back(simplify(ExprFactory::multiply({numToExpr(coeff), rest})));
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

// ============================================================
// Example usage
// ============================================================

int main()
{
    // Example 1:
    // (3 + 2i)^2 + (x + x) + (6*x*y) + (2*x*y) + sqrt(9)
    ExprPtr i = ExprFactory::constant("i");
    ExprPtr x = ExprFactory::symbol("x");
    ExprPtr y = ExprFactory::symbol("y");
    ExprPtr pi = ExprFactory::constant("pi");
    ExprPtr e = ExprFactory::constant("e");

    ExprPtr complexPart = ExprFactory::power(
        ExprFactory::add({ExprFactory::integer(3),
                          ExprFactory::multiply({ExprFactory::integer(2), i})}),
        ExprFactory::integer(2));

    ExprPtr xPlusX = ExprFactory::add({x, x});
    ExprPtr sixxy = ExprFactory::multiply({ExprFactory::integer(6), x, y});
    ExprPtr twoxy = ExprFactory::multiply({ExprFactory::integer(2), x, y});
    ExprPtr sqrt9 = ExprFactory::function("sqrt", {ExprFactory::integer(9)});

    ExprPtr expr = ExprFactory::add({complexPart, xPlusX, sixxy, twoxy, sqrt9});

    std::cout << "Original tree (preorder):\n";
    expr->printPreorder();

    ExprPtr simp = Simplifier::simplify(expr);

    std::cout << "\nSimplified tree (preorder):\n";
    simp->printPreorder();

    std::cout << "\nCanonical string form:\n";
    std::cout << simp->toString() << "\n";

    // Example 2: Euler form
    ExprPtr euler = ExprFactory::power(
        e,
        ExprFactory::multiply({i, pi}));

    std::cout << "\nEuler form e^(i*pi):\n";
    auto eulerS = Simplifier::simplify(euler);
    eulerS->printPreorder();
    std::cout << eulerS->toString() << "\n";

    // Example 3: roots and complex power
    ExprPtr rootExpr = ExprFactory::function("sqrt", {ExprFactory::integer(-1)});
    // ExprPtr powExpr = ;

    std::cout << "\nSqrt(-1):\n";
    auto rootS = Simplifier::simplify(rootExpr);
    rootS->printPreorder();
    std::cout << rootS->toString() << "\n";

    std::cout << "\n(1+i)^5:\n";
    auto powS = Simplifier::simplify(ExprFactory::power(
        ExprFactory::add({ExprFactory::integer(1), i}),
        ExprFactory::integer(5)));
    powS->printPreorder();
    std::cout << powS->toString() << "\n";

    std::cout << "\nRoot 2 by 2:\n";
    auto powS2 = Simplifier::simplify(ExprFactory::multiply({ExprFactory::power(ExprFactory::integer(2), ExprFactory::rational(1, 2)), ExprFactory::integer(2)}));
    powS2->printPreorder();
    std::cout << powS2->toString() << "\n";

    return 0;
}