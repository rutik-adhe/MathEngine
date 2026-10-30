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
#include <stack>
#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/integer/common_factor.hpp>

using boost::multiprecision::cpp_int;
using std::make_shared;
using std::shared_ptr;
using std::string;
using std::vector;

// Exact Rational

class Rational
{
public:
    cpp_int n; // numerator
    cpp_int d; // denominator, always positive

    Rational(cpp_int num = 0, cpp_int den = 1) : n(num), d(den) { normalize(); }

    void normalize()
    {
        if (d == 0)
            throw std::runtime_error("Rational denominator cannot be zero");

        if (d < 0)
        {
            n = -n;
            d = -d;
        }
        cpp_int g = boost::integer::gcd(
            cpp_int(boost::multiprecision::abs(n)),
            cpp_int(boost::multiprecision::abs(d)));
        if (g == 0)
        {
            // handles the case n == 0
            d = 1;
            return;
        }
        n /= g;
        d /= g;
    }

    void inverse()
    {
        cpp_int temp = n;
        n = d;
        d = temp;
    }

    bool isZero() const { return n == 0; }
    bool isOne() const { return n == 1 && d == 1; }
    bool isMinusOne() const { return n == -1 && d == 1; }
    bool isInteger() const { return d == 1; }
    bool isPositive() const { return n > 0; }
    bool isNegative() const { return n < 0; }

    string str() const
    {
        if (d == 1)
            return n.str();
        return n.str() + "/" + d.str();
    }

    friend bool operator==(const Rational &a, const Rational &b)
    {
        return a.n == b.n && a.d == b.d;
    }
    friend bool operator<(const Rational &a, const Rational &b)
    {
        return a.n * b.d < b.n * a.d;
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

// Complex exact (Gaussian rational)

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

// Expression hierarchy

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

    explicit Expr(Kind k, bool s = false) : kind(k), simplified(s) {}
    virtual ~Expr() = default;

    Kind getKind() const { return kind; }

    bool isSimplified() const { return simplified; }

    virtual shared_ptr<Expr> clone() const = 0;
    virtual string toString() const = 0;
    virtual void printPreorder(int indent = 0) const = 0;

private:
    Kind kind;
    bool simplified;
};

using ExprPtr = shared_ptr<Expr>;

class IntegerNode : public Expr
{
public:
    cpp_int value;
    explicit IntegerNode(cpp_int v) : Expr(Kind::Integer, true), value(v) {}
    shared_ptr<Expr> clone() const override { return make_shared<IntegerNode>(value); }
    string toString() const override { return value.str(); }
    void printPreorder(int indent = 0) const override
    {
        std::cout << string(indent, ' ') << "Integer(" << value << ")\n";
    }
};

class RationalNode : public Expr
{
public:
    Rational value;
    explicit RationalNode(const Rational &r) : Expr(Kind::Rational, true), value(r) {}
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
    explicit SymbolNode(string n) : Expr(Kind::Symbol, true), name(std::move(n)) {}
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
    explicit ConstantNode(string n) : Expr(Kind::Constant, true), name(std::move(n)) {}
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
    explicit AddNode(vector<ExprPtr> c, bool s = false) : Expr(Kind::Add, s), children(std::move(c)) {}
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
    explicit MultiplyNode(vector<ExprPtr> c, bool s = false) : Expr(Kind::Multiply, s), children(std::move(c)) {}
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
    PowerNode(ExprPtr b, ExprPtr e, bool s = false) : Expr(Kind::Power, s), base(std::move(b)), exponent(std::move(e)) {}
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
    FunctionNode(string n, vector<ExprPtr> a) : Expr(Kind::Function, true), name(std::move(n)), args(std::move(a)) {}
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

extern vector<std::pair<cpp_int, int>> factorize(cpp_int x);

template <typename T>
T powerComputer(T base, const cpp_int &exp)
{
    // DANGER Assuming base and exp are valid numbers;
    T result{1};
    cpp_int e = boost::multiprecision::abs(exp);
    while (e > 0)
    {
        if (e & 1LL)
            result = result * base;
        base = base * base;
        e >>= 1LL;
    }
    return result;
}

class ExprFactory
{
public:
    static ExprPtr integer(cpp_int v) { return make_shared<IntegerNode>(v); }
    static ExprPtr rational(cpp_int n, cpp_int d)
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
    static ExprPtr add(const vector<ExprPtr> &children, bool simplified = false)
    {
        if (children.size() == 1)
            return children[0];
        return make_shared<AddNode>(children, simplified);
    }
    static ExprPtr multiply(const vector<ExprPtr> &children, bool simplified = false)
    {
        if (children.size() == 1)
            return children[0];
        return make_shared<MultiplyNode>(children, simplified);
    }
    static ExprPtr power(ExprPtr base, ExprPtr exponent, bool simplified = false)
    {
        return make_shared<PowerNode>(std::move(base), std::move(exponent), simplified);
    }
    static ExprPtr function(const string &name, const vector<ExprPtr> &args)
    {
        return make_shared<FunctionNode>(name, args);
    }
};

// ============================================================
// Function Pointers
// ============================================================

ExprPtr (*intg)(cpp_int) = ExprFactory::integer;
ExprPtr (*ratg)(cpp_int, cpp_int) = ExprFactory::rational;
ExprPtr (*powg)(ExprPtr, ExprPtr, bool) = ExprFactory::power;
ExprPtr (*addg)(const vector<ExprPtr> &, bool) = ExprFactory::add;
ExprPtr (*mulg)(const vector<ExprPtr> &, bool) = ExprFactory::multiply;

// ============================================================
// Commonly used Values
// ============================================================

ExprPtr INT_MINUS_THREE = ExprFactory::integer(-3);
ExprPtr INT_MINUS_TWO = ExprFactory::integer(-2);
ExprPtr INT_MINUS_ONE = ExprFactory::integer(-1);
ExprPtr INT_ZERO = ExprFactory::integer(0);
ExprPtr INT_ONE = ExprFactory::integer(1);
ExprPtr INT_TWO = ExprFactory::integer(2);
ExprPtr INT_THREE = ExprFactory::integer(3);
ExprPtr IOTA = ExprFactory::constant("i");
ExprPtr PIE = ExprFactory::constant("pi");
ExprPtr EULER = ExprFactory::constant("e");

ExprPtr FRACT_ONE_HALF = ExprFactory::rational(1, 2);

ExprPtr FRACT_MINUS_ONE_HALF = ExprFactory::rational(-1, 2);

ExprPtr FRACT_ROOT3_BYTWO = ExprFactory::multiply({ExprFactory::power(INT_THREE, FRACT_ONE_HALF, true),
                                                   FRACT_ONE_HALF},
                                                  true);

ExprPtr FRACT_MINUS_ROOT3_BYTWO = ExprFactory::multiply({INT_MINUS_ONE,
                                                         ExprFactory::power(INT_THREE, FRACT_ONE_HALF, true),
                                                         FRACT_ONE_HALF},
                                                        true);

ExprPtr FRACT_ONEBY_ROOT2 = ExprFactory::power(INT_TWO, FRACT_MINUS_ONE_HALF, true);

ExprPtr FRACT_MINUS_ONEBY_ROOT2 = ExprFactory::multiply({INT_MINUS_ONE,
                                                         ExprFactory::power(INT_TWO, FRACT_MINUS_ONE_HALF, true)},
                                                        true);

// ============================================================
// Operator Overloading
// ============================================================

ExprPtr operator*(ExprPtr a, ExprPtr b)
{
    return mulg({a, b}, false);
}

ExprPtr operator*(cpp_int a, ExprPtr b)
{
    return mulg({intg(a), b}, false);
}

ExprPtr operator*(ExprPtr a, cpp_int b)
{
    return mulg({a, intg(b)}, false);
}

ExprPtr operator^(ExprPtr a, ExprPtr b)
{
    return powg(a, b, false);
}

ExprPtr operator^(cpp_int a, ExprPtr b)
{
    return powg(intg(a), b, false);
}

ExprPtr operator^(ExprPtr a, cpp_int b)
{
    return powg(a, intg(b), false);
}

ExprPtr operator-(ExprPtr a, ExprPtr b)
{
    return addg({a, mulg({INT_MINUS_ONE, b}, false)}, false);
}

ExprPtr operator+(ExprPtr a, ExprPtr b)
{
    return addg({a, b}, false);
}

ExprPtr operator/(ExprPtr a, cpp_int b)
{
    return mulg({a, ratg(1, b)}, false);
}

// ============================================================
// Simplifier
// ============================================================

class MultinomialGenerator
{
    cpp_int n;
    vector<ExprPtr> vars;

public:
    MultinomialGenerator(cpp_int degree,
                         vector<ExprPtr> variables)
        : n(degree),
          vars(variables)
    {
    }

    ExprPtr generate()
    {
        vector<cpp_int> exp(vars.size());

        dfs(0, n, 1ULL, exp);

        if (mulNodes.empty())
        {
            return ExprFactory::integer(0);
        }
        else if (mulNodes.size() == 1)
        {
            return mulNodes[0];
        }
        else
        {
            return ExprFactory::add(mulNodes);
        }
    }

private:
    vector<ExprPtr> mulNodes;

    void dfs(size_t pos,
             cpp_int remaining,
             cpp_int coeff,
             vector<cpp_int> &exp)
    {

        // m is total number of variables
        cpp_int m = vars.size();

        if (pos == m - 1)
        {
            exp[pos] = remaining;
            mulNodes.push_back(printTerm(coeff, exp));
            return;
        }

        cpp_int currentCoeff = coeff;

        exp[pos] = remaining;

        dfs(pos + 1,
            0,
            currentCoeff,
            exp);

        for (cpp_int k = remaining - 1; k >= 0; --k)
        {
            currentCoeff =
                currentCoeff *
                (cpp_int)(k + 1) /
                (cpp_int)(remaining - k);

            exp[pos] = k;

            dfs(pos + 1,
                remaining - k,
                currentCoeff,
                exp);
        }
    }

    ExprPtr printTerm(cpp_int coeff,
                      const vector<cpp_int> &exp)
    {
        std::cout << coeff << '\n';

        vector<ExprPtr> nodes;

        nodes.push_back(ExprFactory::integer(coeff));

        for (size_t i = 0; i < exp.size(); ++i)
        {
            if (exp[i] == 0)
                continue;

            else if (exp[i] == 1)
            {
                nodes.push_back(vars[i]);
            }
            else
            {
                nodes.push_back(ExprFactory::power(vars[i], ExprFactory::integer(exp[i])));
            }
        }

        return ExprFactory::multiply(nodes);
    }
};

class Simplifier
{
public:
    static ExprPtr simplify(const ExprPtr &e)
    {
        if (!e)
            return nullptr;

        if (e->isSimplified())
            return e;

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

    static ExprPtr unitComplexNumberGen(Rational coeff)
    {
        // assuming the coeff = coefficient of PIE , that is theta = coeff * PIE
        // coeff belongs to [0, 2)

        // Handle negative angles, multiples and num > deno;

        if (coeff.isInteger())
        {
            coeff = Rational(coeff.n % 2);
        }

        cpp_int intPow = coeff.n / coeff.d;

        if (intPow)
        {
            Rational reducedCoeff(coeff.n - intPow * coeff.d, coeff.d);
            ExprPtr complexUnitTemp = unitComplexNumberGen(reducedCoeff);

            if (intPow % 2 == 0)
            {
                return complexUnitTemp;
            }
            else
            {
                return Simplifier::simplify(ExprFactory::multiply({INT_MINUS_ONE, complexUnitTemp}));
            }
        }
        else
        {
            // Convert to positive case
            if (coeff.isNegative())
            {
                return unitComplexNumberGen(Rational(2) + coeff);
            }

            // Egde cases
            if (coeff.isZero())
            {
                return INT_ONE;
            }
            else if (coeff.isOne())
            {
                return INT_MINUS_ONE;
            }

            if (Rational(1, 2) < coeff)
            {
                // Determine the Quadrant II CASE
                return Simplifier::simplify(ExprFactory::multiply({IOTA, unitComplexNumberGen(coeff - Rational(1, 2))}));
            }

            ExprPtr result = NULL;

            if (coeff == Rational(1, 6))
            {
                result = ExprFactory::add({FRACT_ROOT3_BYTWO, ExprFactory::multiply({FRACT_ONE_HALF, IOTA}, true)}, true);
            }
            else if (coeff == Rational(1, 4))
            {
                result = ExprFactory::add({FRACT_ONEBY_ROOT2, ExprFactory::multiply({FRACT_ONEBY_ROOT2, IOTA}, true)}, true);
            }
            else if (coeff == Rational(1, 3))
            {
                result = ExprFactory::add({FRACT_ONE_HALF, ExprFactory::multiply({FRACT_ROOT3_BYTWO, IOTA}, true)}, true);
            }
            else if (coeff == Rational(1, 2))
            {
                result = IOTA;
            }
            else
            {
                result = ExprFactory::power(EULER, ExprFactory::multiply({IOTA, PIE, ExprFactory::rational(coeff)}, true), true);
            }

            return result;
        }
    }

    static bool isZeroNode(ExprPtr ptr)
    {
        bool isZero = false;
        Expr::Kind coeffExprKind = ptr->getKind();
        if (coeffExprKind == Expr::Kind::Integer)
        {
            auto node = std::static_pointer_cast<IntegerNode>(ptr);
            if (node->value == 0)
                isZero = true;
        }
        return isZero;
    }

private:
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

    static ExprPtr simplifyFunction(const std::shared_ptr<FunctionNode> &f)
    {
        throw std::runtime_error("Invaild Call to Simplify Function");
    }

    static ExprPtr powIntInt(cpp_int base, const cpp_int &exp)
    {
        if (exp > 1)
        {
            if (base == 1)
            {
                return INT_ONE;
            }
            else if (base == 0)
            {
                return INT_ZERO;
            }
            else if (base == -1)
            {
                cpp_int e = boost::multiprecision::abs(exp);
                if (e % 2 == 0)
                {
                    return INT_ONE;
                }
                else
                    return INT_MINUS_ONE;
            }
            else
            {
                return ExprFactory::integer(powerComputer<cpp_int>(base, exp));
            }
        }
        else if (exp == 1)
        {
            return ExprFactory::integer(base);
        }
        else if (exp == 0)
        {
            if (base == 0)
            {
                throw std::runtime_error("Cannot Raise (Int) 0 to the power (Int) 0");
            }
            else
            {
                return INT_ONE;
            }
        }
        else
        {
            // Exponents (Negative Integers)
            if (base == 0)
                throw std::runtime_error("Cannot Divide by (Int) 0");

            // Special Case exp == -1
            if (exp == -1)
            {
                // We have to inverse only
                return ExprFactory::rational(1, base);
            }

            return ExprFactory::rational(1, powerComputer<cpp_int>(base, -exp));
        }
    }

    static ExprPtr powRationalInt(const Rational &base, const cpp_int &exp)
    {

        if (exp > 1)
        {
            if (base.isOne())
            {
                return INT_ONE;
            }
            else if (base.isZero())
            {
                return INT_ZERO;
            }
            else if (base.isMinusOne())
            {
                cpp_int e = boost::multiprecision::abs(exp);
                if (e % 2 == 0)
                {
                    return INT_ONE;
                }
                else
                    return INT_MINUS_ONE;
            }
            else
            {
                return ExprFactory::rational(powerComputer<Rational>(base, exp));
            }
        }
        else if (exp == 1)
        {
            return ExprFactory::rational(base);
        }
        else if (exp == 0)
        {
            if (base.isZero())
            {
                throw std::runtime_error("Cannot Raise (Int) 0 to the power (Int) 0");
            }
            else
            {
                return INT_ONE;
            }
        }
        else
        {
            // Exponents (Negative Integers)
            if (base.isZero())
                throw std::runtime_error("Cannot Divide by (Int) 0");

            // Special Case exp == -1
            if (exp == -1)
            {
                // We have to inverse only
                return ExprFactory::rational(base.d, base.n);
            }
            Rational rInverse = powerComputer<Rational>(base, -exp);
            return ExprFactory::rational(rInverse.d, rInverse.n);
        }
    }

    static ExprPtr powIotaInt(const cpp_int &expNum)
    {
        cpp_int m = expNum % 4;
        if (m < 0)
            m += 4;

        if (m == 0)
        {
            return INT_ONE;
        }
        else if (m == 1)
        {
            return IOTA;
        }
        else if (m == 2)
        {
            return INT_MINUS_ONE;
        }
        else if (m == 3)
        {
            return ExprFactory::multiply({INT_MINUS_ONE, IOTA}, true);
        }

        throw std::runtime_error("Don't call this undefined powIotaInt Function");
        return nullptr;
    }

    static vector<ExprPtr> extractorIntRational(const cpp_int &base, const Rational &exp, bool &splitted)
    {
        // std::cout << "Input Base : " << base << " Input Exp : " << exp.n << "/" << exp.d << std::endl;

        // ASSUMING base is int and exp is positive number
        if (exp.isNegative())
            throw std::runtime_error("Cannot Extract using Negative Rational Exponent");
        if (base <= 0)
            throw std::runtime_error("Cannot Extract from negative or Zero Base");

        if (base == 1)
            return {INT_ONE};

        if (base != 0 && exp.isZero())
            return {INT_ONE};

        vector<ExprPtr> newMultiplyNodes;

        if (exp.isInteger())
        {
            return {ExprFactory::integer(powerComputer<cpp_int>(base, exp.n))};
        }

        cpp_int intPow = exp.n / exp.d;

        // std::cout << "intPow : " << intPow << std::endl;

        if (intPow)
        {
            splitted = true;
            // std::cout << "Enter splitting Case " << std::endl;

            newMultiplyNodes.push_back(ExprFactory::integer(powerComputer<cpp_int>(base, intPow)));

            Rational newAbsExp = Rational(exp.n - exp.d * intPow, exp.d);
            vector<ExprPtr> nextMultiplyNodes = extractorIntRational(base, newAbsExp, splitted);

            newMultiplyNodes.insert(newMultiplyNodes.end(), nextMultiplyNodes.begin(), nextMultiplyNodes.end());
            // check factorization with ths newAbsExp and old base Integer
        }
        else
        {
            // std::cout << "Entered No splitting Case " << std::endl;
            // No Spliiting p < q

            vector<std::pair<cpp_int, int>> f = factorize(base);

            cpp_int coeff = 1;

            for (const auto &[prime, cnt] : f)
            {
                // std::cout << "Prime : " << prime << " Count : " << cnt << std::endl;

                Rational newAbsExp = Rational(cnt * exp.n, exp.d);

                if (newAbsExp.n >= newAbsExp.d)
                {
                    splitted = true;
                    vector<ExprPtr> nextMultiplyNodes = extractorIntRational(prime, newAbsExp, splitted);
                    newMultiplyNodes.insert(newMultiplyNodes.end(), nextMultiplyNodes.begin(), nextMultiplyNodes.end());
                }
                else
                {
                    newMultiplyNodes.push_back(ExprFactory::power(ExprFactory::integer(prime), ExprFactory::rational(exp), true));
                }
            }
            // this coeff will get multiplied with those prime whose count was not integral multiple do deno
            // the remaining prime generates a new same question having (+ve base) and (+ve exp) do same
        }

        return newMultiplyNodes;
    }

    static ExprPtr nthRootOfIota(cpp_int N)
    {
        if (N <= 0)
            throw std::runtime_error("Cannot take Negative/Zeroth Root of Iota");

        if (N == 1)
        {
            return IOTA;
        }
        else if (N == 2)
        {
            return ExprFactory::add({ExprFactory::power(ExprFactory::integer(2), ExprFactory::rational(-1, 2), true),
                                     ExprFactory::multiply({ExprFactory::power(ExprFactory::integer(2), ExprFactory::rational(-1, 2), true),
                                                            IOTA})});
        }
        else
        {
            throw std::runtime_error("No Definations : nthRootOfIota");
        }
    }

    static ExprPtr simplifyRational(const std::shared_ptr<RationalNode> &r)
    {
        if (r->value.isInteger())
        {
            return ExprFactory::integer(r->value.n);
        }
        return r;
    }

    static void IntPowSplitter(const ExprPtr num_childrens, cpp_int &numInt, std::map<string, std::pair<cpp_int, ExprPtr>> &numResetInt)
    {
        Expr::Kind num_Kind = num_childrens->getKind();

        if (num_Kind == Expr::Kind::Integer)
        {
            numInt = std::static_pointer_cast<IntegerNode>(num_childrens)->value;
        }
        else if (num_Kind == Expr::Kind::Multiply)
        {

            auto numMulNode = std::static_pointer_cast<MultiplyNode>(num_childrens);

            for (auto const childPtr : numMulNode->children)
            {
                Expr::Kind child_Kind = childPtr->getKind();

                if (child_Kind == Expr::Kind::Integer)
                {
                    numInt *= std::static_pointer_cast<IntegerNode>(childPtr)->value;
                }
                else if (child_Kind == Expr::Kind::Power)
                {
                    auto basePtr = std::static_pointer_cast<PowerNode>(childPtr)->base;
                    auto expoPtr = std::static_pointer_cast<PowerNode>(childPtr)->exponent;

                    auto baseNode = std::static_pointer_cast<IntegerNode>(basePtr);

                    string key = expoPtr->toString();
                    auto it = numResetInt.find(key);

                    if (it == numResetInt.end())
                    {
                        numResetInt[key] = {baseNode->value, expoPtr};
                    }
                    else
                    {
                        numResetInt[key].first *= baseNode->value;
                    }
                }
                else
                {
                    std::cout << std::endl
                              << "child_Kind : " << std::to_string(static_cast<int>(child_Kind)) << std::endl;
                    throw std::runtime_error("child Kind Unknown ");
                }
            }
        }
        else if (num_Kind == Expr::Kind::Power)
        {
            auto basePtr = std::static_pointer_cast<PowerNode>(num_childrens)->base;
            auto expoPtr = std::static_pointer_cast<PowerNode>(num_childrens)->exponent;

            auto baseNode = std::static_pointer_cast<IntegerNode>(basePtr);

            string key = expoPtr->toString();
            auto it = numResetInt.find(key);

            if (it == numResetInt.end())
            {
                numResetInt[key] = {baseNode->value, expoPtr};
            }
            else
            {
                numResetInt[key].first *= baseNode->value;
            }
        }
        else
        {
            std::cout << std::endl
                      << "numKind : " << std::to_string(static_cast<int>(num_Kind)) << std::endl;
            throw std::runtime_error("Num Kind Unknown");
        }
    }

    static ExprPtr simplifyPower(const std::shared_ptr<PowerNode> &p)
    {
        ExprPtr base = simplify(p->base);
        // Solves Base MultiplyNode = sqrt(-2)*sqrt(-3) = -sqrt(6)
        // Solves Base MultiplyNode = sqrt(-4)*sqrt(-3)sqrt(-6) = +6sqrt(2)i
        // Solves Base MultiplyNode = sqrt[(-2)*(-3)] = +sqrt(6)
        // Solves Base MultiplyNode = sqrt[(-4)*(-3)*(-6)] = +6sqrt(2)i

        // Try to take common out before powering upto group terms, reduce computation and then after powering up, insert the common terms inside again for further simplification consistency

        ExprPtr exp = simplify(p->exponent);

        Expr::Kind expKind = exp->getKind();
        Expr::Kind baseKind = base->getKind();

        if (expKind == Expr::Kind::Integer)
        {
            auto expNode = std::static_pointer_cast<IntegerNode>(exp);
            // first check for superier cases like is ZERO/ONE etc.
            // if exp is positive then only do repeativete computation

            // Universal cases
            if (expNode->value == 1)
                return base;

            if (expNode->value == 0 && (baseKind != Expr::Kind::Integer))
                return INT_ONE;

            // Different Cases for each base type
            if (baseKind == Expr::Kind::Integer)
            {
                auto baseNode = std::static_pointer_cast<IntegerNode>(base);

                return powIntInt(baseNode->value, expNode->value);
            }
            else if (baseKind == Expr::Kind::Rational)
            {
                auto baseNode = std::static_pointer_cast<RationalNode>(base);
                return powRationalInt(baseNode->value, expNode->value);
            }
            else if (baseKind == Expr::Kind::Constant)
            {
                auto baseNode = std::static_pointer_cast<ConstantNode>(base);
                if (baseNode->name == "i")
                {
                    return powIotaInt(expNode->value);
                }
                else
                {
                    return ExprFactory::power(base, exp, true);
                }
            }
            else if (baseKind == Expr::Kind::Symbol)
            {

                return ExprFactory::power(base, exp, true);
            }
            else if (baseKind == Expr::Kind::Complex)
            {
                auto baseNode = std::static_pointer_cast<ComplexNode>(base);
                throw std::runtime_error("No Implementation for base=complex in Power Function");
            }
            else if (baseKind == Expr::Kind::Add)
            {
                auto baseNode = std::static_pointer_cast<AddNode>(base);
                // Multinomial Theorem expansion

                cpp_int expValue = expNode->value;

                if (expValue > 1)
                {
                    // std::cout << "Doing Multinomial Expansion of Degree : " << expNode->value << " & Terms : " << baseNode->children.size() << std::endl;

                    MultinomialGenerator gen(expNode->value, baseNode->children);

                    return simplify(gen.generate());
                }
                else if (expValue == 1)
                {
                    return base;
                }
                else if (expValue == 0)
                {
                    return INT_ONE;
                }
                else if (expValue == -1)
                {
                    return ExprFactory::power(base, exp, true);
                }
                else
                {
                    // std::cout << "Doing Multinomial Expansion of Degree : " << expNode->value << " & Terms : " << baseNode->children.size() << std::endl;

                    MultinomialGenerator gen(-(expNode->value), baseNode->children);

                    return ExprFactory::power(simplify(gen.generate()), INT_MINUS_ONE);
                }
            }
            else if (baseKind == Expr::Kind::Multiply)
            {
                auto baseNode = std::static_pointer_cast<MultiplyNode>(base);

                vector<ExprPtr> raisedChildrens;

                for (auto const c : baseNode->children)
                {
                    raisedChildrens.push_back(simplify(ExprFactory::power(c, exp)));
                }

                return simplify(ExprFactory::multiply(raisedChildrens));
            }
            else if (baseKind == Expr::Kind::Power)
            {
                auto powerNode = std::static_pointer_cast<PowerNode>(base);

                // std::cout << "Reached this case = power^int" << std::endl;

                //  (2^a)^b = 2^(a.b)
                return simplify(ExprFactory::power(powerNode->base, simplify(ExprFactory::multiply({exp, powerNode->exponent}))));
            }
            else if (baseKind == Expr::Kind::Function)
            {
                // sin1^2, tan5^2, etc.
                return ExprFactory::power(base, exp, true);
            }
        }
        else if (expKind == Expr::Kind::Rational)
        {
            auto expNode = std::static_pointer_cast<RationalNode>(exp);

            if ((expNode->value).isInteger())
            {
                // edge cases of base 0 or 1 are handled inside the Exponent=Integer section of simplifyPower
                return simplify(ExprFactory::power(base, ExprFactory::integer(expNode->value.n)));
            }
            else
            {
                if (baseKind == Expr::Kind::Integer)
                {
                    auto baseNode = std::static_pointer_cast<IntegerNode>(base);

                    // std::cout << std::endl
                    //           << "(int)^(rat) Function Input Exponent : " << expNode->value.n << "/" << expNode->value.d << std::endl
                    //           << std::endl;

                    bool isNegExp = (expNode->value).isNegative();
                    Rational absExp = isNegExp ? -(expNode->value) : expNode->value;

                    // Handle many edge cases of base
                    if (baseNode->value > 1)
                    {
                        bool splitted = false;
                        vector<ExprPtr> extractorResultPtr = extractorIntRational(baseNode->value, absExp, splitted);

                        // Given Input (int)^(rational)
                        // Possible Outputs :- MultiplyNode or

                        if (extractorResultPtr.empty())
                        {
                            throw std::runtime_error("ExtractorIntRational No Results");
                        }

                        // std::cout << "Splitted " << splitted << std::endl;

                        if (splitted)
                        {
                            if (extractorResultPtr.size() == 1)
                            {
                                // std::cout << "One Element case observered" << std::endl;
                                if (isNegExp)
                                {
                                    return simplify(ExprFactory::power(extractorResultPtr[0], INT_MINUS_ONE));
                                }
                                else
                                {
                                    return extractorResultPtr[0];
                                }
                            }
                            else
                            {
                                if (isNegExp)
                                {
                                    return simplify(ExprFactory::power(ExprFactory::multiply(extractorResultPtr, true), INT_MINUS_ONE));
                                }
                                else
                                {
                                    return simplify(ExprFactory::multiply(extractorResultPtr));
                                }
                            }
                        }
                        else
                        {
                            if (isNegExp)
                            {
                                return ExprFactory::power(simplify(ExprFactory::power(baseNode, INT_MINUS_ONE)), ExprFactory::rational(absExp), true);
                            }
                            else
                            {
                                return ExprFactory::power(base, exp, true);
                            }
                        }
                    }
                    else if (baseNode->value == 1)
                    {
                        return INT_ONE;
                    }
                    else if (baseNode->value == 0)
                    {
                        if (expNode->value.isNegative())
                        {
                            throw std::runtime_error("Cannot Raise Negative (Rational) Power to 0 Base");
                        }
                        return INT_ZERO;
                    }
                    else if (baseNode->value == -1)
                    {
                        // NOTE (-1)^(-1) is same as (-1)
                        // Handle on bases of odd and even, and return the remaining part as it. Or check if exp = 1/2

                        if (absExp.n % 2 == 0)
                        {
                            return INT_ONE;
                        }
                        else
                        {
                            if (absExp.d == 2)
                            {
                                return IOTA;
                            }
                            else
                            {
                                // Nth Root of (-1) in terms of IOTA
                                return ExprFactory::power(IOTA, ExprFactory::rational(2, absExp.d), true);
                            }
                        }
                    }
                    else
                    {
                        // Negative Bases lessen than -1
                        cpp_int absBase = -(baseNode->value);

                        ExprPtr result1 = simplify(ExprFactory::power(ExprFactory::integer(absBase), exp));

                        if (absExp.d % 2 == 1)
                        {
                            // Odd -th Root of negative number
                            return simplify(ExprFactory::multiply({INT_MINUS_ONE, result1}));
                        }
                        else
                        {
                            // when even -th Root of negative number
                            ExprPtr iotaExpr = simplify(ExprFactory::power(IOTA, ExprFactory::rational(2, absExp.d)));
                            return simplify(ExprFactory::multiply({iotaExpr, result1}));
                        }
                    }
                }
                else if (baseKind == Expr::Kind::Rational)
                {
                    // exponent Rational Pure (Not Integer) and Base also Rational.
                    auto baseNode = std::static_pointer_cast<RationalNode>(base);
                    if ((baseNode->value).isInteger())
                    {
                        // edge cases of base 0 or 1 are handled inside the Exponent=Integer section of simplifyPower
                        return simplify(ExprFactory::power(ExprFactory::integer(baseNode->value.n), expNode));
                    }
                    else
                    {
                        // BOTH PURE RATIONALS
                        // When Exponent is negative Rational Inverse the numerator terms instead of denominator terms
                        // The numerator is itself a question of (int)^(Rational)

                        // std::cout << std::endl  << "REAACHED POWER END" << std::endl;

                        // HANDLE THE CASE when Base rational is also negative, its the same case as (-1)^(Rational)

                        bool isNegExp = (expNode->value).isNegative();
                        Rational absExp = isNegExp ? -(expNode->value) : expNode->value;

                        bool isNegBase = (baseNode->value).isNegative();
                        Rational absBase = isNegBase ? -(baseNode->value) : baseNode->value;

                        ExprPtr num_childrens = simplify(ExprFactory::power(ExprFactory::integer(absBase.n), ExprFactory::rational(absExp)));
                        // std::cout << "Numerator : " << absBase.n << " Rational Power : " << absExp.n << "/" << absExp.d << std::endl;

                        ExprPtr deno_childrens = simplify(ExprFactory::power(ExprFactory::integer(absBase.d), ExprFactory::rational(absExp)));
                        // std::cout << "Denominator : " << absBase.d << " Rational Power : " << absExp.n << "/" << absExp.d << std::endl;

                        // the output should be two rationals (rootExtractedOne)(non-extractable)

                        cpp_int numInt = 1;
                        cpp_int denoInt = 1;

                        std::map<string, std::pair<cpp_int, ExprPtr>> numResetInt;
                        std::map<string, std::pair<cpp_int, ExprPtr>> denResetInt;

                        IntPowSplitter(num_childrens, numInt, numResetInt);
                        IntPowSplitter(deno_childrens, denoInt, denResetInt);

                        Rational ExtractedPart(numInt, denoInt);

                        if (isNegExp)
                            ExtractedPart.inverse();

                        vector<ExprPtr> resultChildrens;

                        ExprPtr absExpNode = ExprFactory::rational(absExp);

                        resultChildrens.push_back(ExprFactory::rational(ExtractedPart));

                        for (const auto pair : numResetInt)
                        {

                            auto denoKey = denResetInt.find(pair.first);

                            if (denoKey == denResetInt.end())
                            {
                                if (isNegExp)
                                {
                                    resultChildrens.push_back(ExprFactory::power(ExprFactory::rational(1, pair.second.first), pair.second.second, true));
                                }
                                else
                                {
                                    resultChildrens.push_back(ExprFactory::power(ExprFactory::integer(pair.second.first), pair.second.second, true));
                                }
                            }
                            else
                            {
                                auto pair2_Value = denResetInt[pair.first];

                                if (isNegExp)
                                {
                                    resultChildrens.push_back(ExprFactory::power(ExprFactory::rational(pair2_Value.first, pair.second.first), pair.second.second, true));
                                }
                                else
                                {
                                    resultChildrens.push_back(ExprFactory::power(ExprFactory::rational(pair.second.first, pair2_Value.first), pair.second.second, true));
                                }
                            }
                        }

                        if (isNegBase)
                        {
                            ExprPtr minusOneExpr = simplify(ExprFactory::power(INT_MINUS_ONE, ExprFactory::rational(absExp)));
                            resultChildrens.push_back(minusOneExpr);
                        }

                        // std::cout << std::endl  << "REAACHED POWER END" << std::endl;

                        return simplify(ExprFactory::multiply(resultChildrens));
                    }
                }

                else if (baseKind == Expr::Kind::Constant)
                {
                    auto baseNode = std::static_pointer_cast<ConstantNode>(base);
                    if (baseNode->name == "i")
                    {
                        // deno-th root of IOTA^Num
                        return simplify(ExprFactory::power(powIotaInt(expNode->value.n), ExprFactory::integer(expNode->value.d)));
                    }
                    else
                    {
                        return ExprFactory::power(base, exp, true);
                    }
                }
                else if (baseKind == Expr::Kind::Symbol)
                {
                    // unchanged
                    return ExprFactory::power(base, exp, true);
                }
                else if (baseKind == Expr::Kind::Complex)
                {
                    // COMPLEX IS DATATYPE IS DELETE
                    // IT IS MERGED INTO FUNDAMENTAL OPERATIONS
                    throw std::runtime_error("Complex DATATYPE OCCURED");
                }
                else if (baseKind == Expr::Kind::Add)
                {
                    // extract the integral powers from rational power

                    // bool isNegExp = expNode->value.isNegative();
                    Rational expValue = expNode->value;

                    cpp_int intPow = expValue.n / expValue.d;
                    if (intPow)
                    {
                        Rational reducedExp(expValue.n - intPow * expValue.d, expValue.d);

                        return ExprFactory::multiply({simplify(ExprFactory::power(base, ExprFactory::integer(intPow))),
                                                      ExprFactory::power(base, ExprFactory::rational(reducedExp), true)},
                                                     true);
                    }
                    else
                    {
                        return ExprFactory::power(base, exp, true);
                    }
                }
                else if (baseKind == Expr::Kind::Multiply)
                {
                    auto baseNode = std::static_pointer_cast<MultiplyNode>(base);

                    vector<ExprPtr> raisedChildrens;

                    for (auto const c : baseNode->children)
                    {
                        raisedChildrens.push_back(simplify(ExprFactory::power(c, exp)));
                    }

                    return simplify(ExprFactory::multiply(raisedChildrens));
                }
                else if (baseKind == Expr::Kind::Power)
                {
                    auto baseNode = std::static_pointer_cast<PowerNode>(base);

                    //  (2^a)^b = 2^(a.b)
                    return simplify(ExprFactory::power(baseNode, simplify(ExprFactory::multiply({exp, baseNode->exponent}))));
                }
                else if (baseKind == Expr::Kind::Function)
                {
                    // sin1^2, tan5^2, etc.
                    return ExprFactory::power(base, exp, true);
                }
            }
        }

        else if (expKind == Expr::Kind::Complex)
        {
            // COMPLEX IS DATATYPE IS DELETE
            // IT IS MERGED INTO FUNDAMENTAL OPERATIONS
            throw std::runtime_error("Complex DATATYPE OCCURED");
        }
        else if (expKind == Expr::Kind::Constant)
        {
            // we need to define cos and sine functions and their basic values at known angles
            // we need to define ln(x) to support real to the power IOTA and its known values at x
            // NOTE 0^i doesn't exists

            auto expNode = std::static_pointer_cast<ConstantNode>(base);
            if (expNode->name == "i")
            {
                Expr::Kind baseKind = base->getKind();
                if (baseKind == Expr::Kind::Integer)
                {
                    auto baseNode = std::static_pointer_cast<IntegerNode>(base);

                    if (baseNode->value > 1)
                    {
                        // introduce ln(x), sin(x), cos(x) in order to solve this
                        return ExprFactory::power(base, exp, true);
                    }
                    else if (baseNode->value == 1)
                    {
                        return INT_ONE;
                    }
                    else if (baseNode->value == 0)
                    {
                        throw std::runtime_error("Cannot raise IOTA to base (int) 0-Zero");
                    }
                    else
                    {
                        // there are infinite many solutions for negative integer bases raised to IOTA Power
                        return ExprFactory::power(base, exp, true);
                    }
                }
            }
            else
            {
                // Handle the cases 0^pi, 1^pi, (-2)^pi = (2^2pi)*(i^2pi)
                return ExprFactory::power(base, exp, true);
            }
        }
        else if (expKind == Expr::Kind::Symbol)
        {
            // 2^x,  (root 2)^y , etc terms fall under this category which are unsolvable, they represent an equation.
            return ExprFactory::power(base, exp, true);
        }
        else if (expKind == Expr::Kind::Add)
        {
            // Distributive Property of Exponents
            auto expNode = std::static_pointer_cast<AddNode>(base);

            vector<ExprPtr> newNodes;

            for (auto const &e : expNode->children)
            {
                newNodes.push_back(simplify(ExprFactory::power(base, e)));
            }

            return simplify(ExprFactory::multiply(newNodes));
        }
        else if (expKind == Expr::Kind::Multiply)
        {
            // Int or Rational node from it can be solved reset remain as it
            // combinations Int/Rational * Symbol/Constant/Function/Power with Rational exponent

            // std::cout << "---- Entered the case Exp=MulKind ----" << std::endl;

            return ExprFactory::power(base, exp, true);
        }
        else if (expKind == Expr::Kind::Power)
        {
            // since we have already simplified exponent it means the
            // base - Symbol/Constant/Function
            // exponent - Symbol/Constant/Function/Power with Rational exponent and positive base
            // the base will not be negative because it will be simplified to other cases, IOTA root case, then it is expKind = Multiply

            // exp == 2^root(2), 2^i, i^root(2), e^pi, sin1^i, etc.
            // there are different answers depending upon POWER BASE AND EXP

            // When base is Integer positive and exponent is rational positive
            // Like 2^(1/3) == 2^(3^-1) is unchanged you cannot simplify furthur
            Expr::Kind baseKind = base->getKind();
            if (baseKind == Expr::Kind::Integer)
            {
                auto baseNode = std::static_pointer_cast<IntegerNode>(base);
                if (baseNode->value > 1)
                {
                    // in order to do this, it should have intelligence to select between the exp form or normal form so that final expression is most simplified
                    // a^(root(2)) = e^(root(2)*ln(a))
                    // return ExprFactory::power(EULER, ExprFactory::multiply({base, ExprFactory::function("ln", {exp})}));

                    return ExprFactory::power(base, exp, true);
                }
                else if (baseNode->value == 1)
                    return INT_ONE;

                else if (baseNode->value == 0)
                {
                    // Only if r > 0 positive real quantity
                    return INT_ZERO;
                }
                else
                    return ExprFactory::power(base, exp, true);
            }
            return ExprFactory::power(base, exp, true);
        }

        // std::cout << std::endl << "BaseKind : " << std::to_string(static_cast<int>(baseKind)) << "ExpKind : " << std::to_string(static_cast<int>(expKind)) << std::endl;

        throw std::runtime_error("POWER Case Not Found");
    }

    static ExprPtr simplifyMultiply(const std::shared_ptr<MultiplyNode> &m)
    {
        // collect the terms with same base in simplifyMultiply function treat each Bracket in multiply as base

        // remove unneccessary bracketing inside multiply operator ((a).(b).(c)) = (a.b.c)
        // cpp_int resultInt = 1;
        Rational resultRational(1);

        vector<ExprPtr> nonNumericFactorsPairs;

        // since we are only using factors vector only to solve for size(2) case, we can optimized this as well
        vector<ExprPtr> factors;
        for (const auto &c : m->children)
        {
            // This simplify on each children ensures flattening of all nested multiple nodes
            ExprPtr s = simplify(c);

            // std::cout << "mulFactor : " << std::to_string(static_cast<int>(s->getKind())) << std::endl;
            // 0-Integer , 1-Rational, 2-Symbol, 3-Constant, 4-Complex, 5-Add, 6-Multiply, 7-Power, 8-Function

            Expr::Kind kindS = s->getKind();
            if (kindS == Expr::Kind::Multiply)
            {
                auto mm = std::static_pointer_cast<MultiplyNode>(s);
                for (const auto &k : mm->children)
                {
                    Expr::Kind kindK = k->getKind();
                    if (kindK == Expr::Kind::Integer)
                    {
                        resultRational = resultRational * Rational(std::static_pointer_cast<IntegerNode>(k)->value);
                    }
                    else if (kindK == Expr::Kind::Rational)
                    {
                        resultRational = resultRational * std::static_pointer_cast<RationalNode>(k)->value;
                    }
                    else
                    {
                        // if symbol, constant, function, add, power
                        nonNumericFactorsPairs.push_back(k);
                    }
                    factors.push_back(k);
                }
            }
            else
            {
                if (kindS == Expr::Kind::Integer)
                {
                    resultRational = resultRational * Rational(std::static_pointer_cast<IntegerNode>(s)->value);
                }
                else if (kindS == Expr::Kind::Rational)
                {
                    resultRational = resultRational * std::static_pointer_cast<RationalNode>(s)->value;
                }
                else
                {
                    nonNumericFactorsPairs.push_back(s);
                }
                factors.push_back(s);
            }
        }

        if (factors.size() == 0)
        {
            // confusing answer 0 or 1
            return INT_ZERO;
        }
        else if (factors.size() == 1)
            return factors[0];
        // else if (factors.size() == 2)
        // {
        //     std::cout << std::endl
        //               << "call TwoTermMultiplier for them" << std::endl;
        //     return TwoTermMultiplier(factors[0], factors[1]);
        // }
        else
        {

            // std::cout << "running MainMultiplier for them" << std::endl;

            // Group terms with same bases and add their power
            // counts how much time each term as repeated other than numbers & rational
            std::map<string, std::pair<ExprPtr, cpp_int>> powerNodesMapInt;
            std::map<string, std::pair<ExprPtr, std::pair<Rational, bool>>> powerNodesMapRational;

            // We don't need to solve powerNodesMapOther, they are other unsolvable power nodes
            std::map<string, std::pair<ExprPtr, std::vector<ExprPtr>>> powerNodesMapOther;

            // possible nodoes = symbol, constants, addNode, powerNode, function

            // out of this 5 nodes constant(i), powerNode (unsolvable power terms like x^2, e^i, (1+root(2))^(1/2), i^(1/2), i^i   ), addNode (unsolvable addNode) can only generate numeric terms
            // symbol and function terms will remain as it is!
            // Reoccurance of other than powerNode produces pure integer power terms
            // we don't do integer exponent extraction for symbols, only for addNodes (powernodes having base addNode)

            for (const auto &factorNode : nonNumericFactorsPairs)
            {

                Expr::Kind factorNodeKind = factorNode->getKind();
                string factorKey = (factorNode)->toString();

                if (factorNodeKind == Expr::Kind::Power)
                {
                    auto powerNode = std::static_pointer_cast<PowerNode>(factorNode);
                    string key = (powerNode->base)->toString();

                    // Their can exists some powerNodes with Integer power and some with rational power and some with other stuff
                    // we need 3 groups within this as well. And our aim is to extract integer powers terms out of those bunch terms.

                    // Below is for BaseKind constant == i, AddNode
                    Expr::Kind baseKind = powerNode->base->getKind();
                    Expr::Kind expKind = powerNode->exponent->getKind();

                    if (baseKind == Expr::Kind::Integer || baseKind == Expr::Kind::Rational || baseKind == Expr::Kind::Add || (baseKind == Expr::Kind::Constant && std::static_pointer_cast<ConstantNode>(powerNode->base)->name == "i"))
                    {
                        if (expKind == Expr::Kind::Integer)
                        {
                            if (baseKind == Expr::Kind::Integer || baseKind == Expr::Kind::Rational)
                            {
                                throw std::runtime_error("Pre-simplify fails - MultiplyFunc");
                            }

                            cpp_int expValue = std::static_pointer_cast<IntegerNode>(powerNode->exponent)->value;

                            auto it = powerNodesMapInt.find(key);
                            if (it == powerNodesMapInt.end())
                            {
                                powerNodesMapInt[key] = {powerNode->base, expValue};
                            }
                            else
                            {
                                powerNodesMapInt[key].second += expValue;
                            }
                        }
                        else if (expKind == Expr::Kind::Rational)
                        {
                            auto it = powerNodesMapRational.find(key);
                            Rational expValue = std::static_pointer_cast<RationalNode>(powerNode->exponent)->value;

                            if (it == powerNodesMapRational.end())
                            {
                                powerNodesMapRational[key] = {powerNode->base, {expValue, false}};
                            }
                            else
                            {
                                powerNodesMapRational[key].second.first = powerNodesMapRational[key].second.first + expValue;
                                powerNodesMapRational[key].second.second = true;
                            }
                        }
                        else
                        {
                            auto it = powerNodesMapOther.find(key);
                            if (it == powerNodesMapOther.end())
                            {
                                powerNodesMapOther[key] = {powerNode->base, {powerNode->exponent}};
                            }
                            else
                            {
                                powerNodesMapOther[key].second.push_back(powerNode->exponent);
                            }
                        }
                    }
                    else if (baseKind == Expr::Kind::Function || baseKind == Expr::Kind::Symbol || baseKind == Expr::Kind::Constant)
                    {
                        // take care about the thing that don't sum exponents of e^(non IOTA) and e^(IOTA) terms

                        auto it = powerNodesMapOther.find(key);
                        if (it == powerNodesMapOther.end())
                        {
                            powerNodesMapOther[key] = {powerNode->base, {powerNode->exponent}};
                        }
                        else
                        {
                            powerNodesMapOther[key].second.push_back(powerNode->exponent);
                        }
                    }
                    else
                    {
                        throw std::runtime_error("Basekkind different observed");
                    }
                }
                else if (factorNodeKind == Expr::Kind::Add || (factorNodeKind == Expr::Kind::Constant && std::static_pointer_cast<ConstantNode>(factorNode)->name == "i"))
                {
                    // possible nodoes = symbol, constants, addNode, function
                    // this block captures repeation of nodes
                    auto it = powerNodesMapInt.find(factorKey);
                    if (it == powerNodesMapInt.end())
                    {
                        powerNodesMapInt[factorKey] = {factorNode, 1};
                    }
                    else
                    {
                        powerNodesMapInt[factorKey].second += 1;
                    }
                }
                else if (factorNodeKind == Expr::Kind::Symbol || factorNodeKind == Expr::Kind::Function || factorNodeKind == Expr::Kind::Constant)
                {
                    auto it = powerNodesMapOther.find(factorKey);
                    if (it == powerNodesMapOther.end())
                    {
                        powerNodesMapOther[factorKey] = {factorNode, {INT_ONE}};
                    }
                    else
                    {
                        powerNodesMapOther[factorKey].second.push_back({INT_ONE});
                    }
                }
                else
                {
                    factorNode->printPreorder();

                    throw std::runtime_error("factorkind different observed");
                }
            }

            // Now Evalute Irrational Maps first find normal terms merge them with Integer Maps, push the still remaining Rational Map to

            vector<ExprPtr> unsolvablePowerTerms;
            for (auto const &mapPair : powerNodesMapRational)
            {

                ExprPtr resultExpr = NULL;
                if (mapPair.second.second.second)
                {
                    resultExpr = simplify(ExprFactory::power(mapPair.second.first, ExprFactory::rational(mapPair.second.second.first)));
                }
                else
                {
                    resultExpr = ExprFactory::power(mapPair.second.first, ExprFactory::rational(mapPair.second.second.first), true);
                }

                Expr::Kind resultExprKind = resultExpr->getKind();

                if (resultExprKind == Expr::Kind::Multiply)
                {
                    auto mulNode = std::static_pointer_cast<MultiplyNode>(resultExpr);

                    for (const auto &child : mulNode->children)
                    {
                        Expr::Kind childExprKind = child->getKind();
                        if (childExprKind == Expr::Kind::Power)
                        {

                            // check if exponent is integer like  (add)^-1 or (add)^(-1/2)
                            // i^(-1/2) * i^(-1/2) will not occur in this case
                            // (add)^(-1/4) * (add)^(-1/4) * (add)^(-1/4) * (add)^(-1/4) * (add)^(-1/4) = (add)^(-1) * (add)^(-1/2)

                            auto powerNode = std::static_pointer_cast<PowerNode>(child);
                            ExprPtr expPtr = powerNode->exponent;
                            if (expPtr->getKind() == Expr::Kind::Integer)
                            {
                                cpp_int expValue = std::static_pointer_cast<IntegerNode>(expPtr)->value;

                                string key = powerNode->base->toString();
                                auto it = powerNodesMapInt.find(key);
                                if (it == powerNodesMapInt.end())
                                {
                                    powerNodesMapInt[key] = {powerNode->base, expValue};
                                }
                                else
                                {
                                    powerNodesMapInt[key].second += expValue;
                                }
                            }
                            else
                            {
                                unsolvablePowerTerms.push_back(child);
                            }
                        }
                        else if (childExprKind == Expr::Kind::Integer)
                        {
                            resultRational = resultRational * Rational(std::static_pointer_cast<IntegerNode>(child)->value);
                        }
                        else if (childExprKind == Expr::Kind::Rational)
                        {
                            resultRational = resultRational * std::static_pointer_cast<RationalNode>(child)->value;
                        }
                        else if (childExprKind == Expr::Kind::Add || childExprKind == Expr::Kind::Constant)
                        {

                            // push to integer maps
                            string key = child->toString();
                            auto it = powerNodesMapInt.find(key);
                            if (it == powerNodesMapInt.end())
                            {
                                powerNodesMapInt[key] = {child, 0};
                            }
                            else
                            {
                                powerNodesMapInt[key].second += 1;
                            }
                        }
                        else
                        {
                            throw std::runtime_error("Program hasn't overed the case : ");
                        }
                    }
                }

                else if (resultExprKind == Expr::Kind::Integer)
                {
                    resultRational = resultRational * Rational(std::static_pointer_cast<IntegerNode>(resultExpr)->value);
                }
                else if (resultExprKind == Expr::Kind::Rational)
                {
                    resultRational = resultRational * std::static_pointer_cast<RationalNode>(resultExpr)->value;
                }
                else
                {
                    unsolvablePowerTerms.push_back(resultExpr);
                }
            }

            vector<std::shared_ptr<AddNode>> addNodesExpOne;
            bool isExternalIota = false;

            for (auto const &mapPair : powerNodesMapInt)
            {
                cpp_int nodePower = mapPair.second.second;

                ExprPtr resultExpr = NULL;
                if (nodePower != 1)
                {
                    resultExpr = simplify(ExprFactory::power(mapPair.second.first, ExprFactory::integer(nodePower)));
                }
                else
                {
                    resultExpr = mapPair.second.first;
                }

                Expr::Kind resultExprKind = resultExpr->getKind();

                if (resultExprKind == Expr::Kind::Multiply)
                {
                    // ex : mul(i, x, i, x, y, i)
                    auto mulNode = std::static_pointer_cast<MultiplyNode>(resultExpr);

                    for (auto const childNode : mulNode->children)
                    {
                        Expr::Kind childNodeKind = childNode->getKind();

                        if (childNodeKind == Expr::Kind::Integer)
                        {
                            resultRational = resultRational * Rational(std::static_pointer_cast<IntegerNode>(childNode)->value);
                        }
                        else if (childNodeKind == Expr::Kind::Constant)
                        {
                            auto constNode = std::static_pointer_cast<ConstantNode>(childNode);

                            if (constNode->name == "i")
                            {
                                unsolvablePowerTerms.push_back(IOTA);
                            }
                            else
                            {
                                throw std::runtime_error("Program hasn't overed the case constant (NOT IOTA) : ");
                            }
                        }
                        else
                        {
                            throw std::runtime_error("Program hasn't overed the case MULTIPLY : ");
                        }
                    }
                }
                else if (resultExprKind == Expr::Kind::Power)
                {
                    throw std::runtime_error("Program hasn't overed the case (NOT POSSIBLE POWER) : ");
                    unsolvablePowerTerms.push_back(resultExpr);
                }
                else if (resultExprKind == Expr::Kind::Integer)
                {
                    resultRational = resultRational * Rational(std::static_pointer_cast<IntegerNode>(resultExpr)->value);
                }
                else if (resultExprKind == Expr::Kind::Rational)
                {
                    resultRational = resultRational * std::static_pointer_cast<RationalNode>(resultExpr)->value;
                }
                else if (resultExprKind == Expr::Kind::Add)
                {
                    addNodesExpOne.push_back(std::static_pointer_cast<AddNode>(resultExpr));
                }
                else if (resultExprKind == Expr::Kind::Constant)
                {
                    auto constNode = std::static_pointer_cast<ConstantNode>(resultExpr);

                    if (constNode->name == "i")
                    {
                        if (isExternalIota == true)
                            throw std::runtime_error("Reoccurance of IOTA");
                        isExternalIota = true;
                    }

                    // prevent the pushBack of IOTA
                    unsolvablePowerTerms.push_back(resultExpr);
                }
                else
                {
                    throw std::runtime_error("Program hasn't overed the case : ");
                }
            }

            ExprPtr numericResult = ExprFactory::rational(resultRational);

            // Single AddNode Generator
            ExprPtr addNodeFinal = NULL; // null is being as flag to check latter

            if (!addNodesExpOne.empty())
            {
                size_t totalAddNodes = addNodesExpOne.size();
                addNodeFinal = addNodesExpOne[0];
                addNodesExpOne.erase(addNodesExpOne.begin());

                for (auto const addNode : addNodesExpOne)
                {
                    vector<ExprPtr> resultsAry;

                    Expr::Kind addNodeFinalKind = addNodeFinal->getKind();

                    for (auto const child1 : addNode->children)
                    {
                        child1->printPreorder();

                        if (addNodeFinalKind == Expr::Kind::Add)
                        {
                            auto const addNodeFinal_Mul = std::static_pointer_cast<MultiplyNode>(addNodeFinal);
                            for (auto const child2 : addNodeFinal_Mul->children)
                            {
                                resultsAry.push_back(simplify(ExprFactory::multiply({child1, child2})));
                            }
                        }
                        else
                        {
                            resultsAry.push_back(simplify(ExprFactory::multiply({child1, addNodeFinal})));
                        }
                    }

                    addNodeFinal = Simplifier::simplify(ExprFactory::add(resultsAry));
                }
            }

            // check if external iota is present is so introduce it in addNodeFinal

            if (addNodeFinal && addNodeFinal->getKind() == Expr::Kind::Add && isExternalIota)
            {

                auto addNodeFinal_Node = std::static_pointer_cast<AddNode>(addNodeFinal);

                vector<ExprPtr> mulNodes;
                for (const auto node : addNodeFinal_Node->children)
                {
                    mulNodes.push_back(Simplifier::simplify(ExprFactory::multiply({node, IOTA})));
                }

                addNodeFinal = Simplifier::simplify(ExprFactory::add(mulNodes));

                auto it = std::find(unsolvablePowerTerms.begin(), unsolvablePowerTerms.end(), IOTA);

                if (it != unsolvablePowerTerms.end())
                {
                    unsolvablePowerTerms.erase(it);
                }
            }

            // AtLast check rational formed is integer, Check if groupresults are empty
            // multiply both Rational and Int and new generaeted integers from other groups , then find if the result of all is integer/rational, return accordingly.

            // result = (int)*(rational)*(unsolvable terms)*(mul(add Terms))*(symbols/function terms)
            vector<ExprPtr> results;

            // std::cout << "Multiplication Numeric : " << resultRational.n << "/" << resultRational.d << std::endl;
            // std::cout << "unsolvablePowerTerms : " << unsolvablePowerTerms.size() << std::endl;
            // std::cout << "addNodeFinal : " << static_cast<int>(addNodeFinal != NULL) << std::endl;
            // std::cout << "powerNodesMapOther : " << powerNodesMapOther.size() << std::endl;
            // std::cout << "powerNodesMapOther : " << powerNodesMapOther.size() << std::endl << std::endl;

            if (!unsolvablePowerTerms.empty())
            {
                results.insert(results.end(), unsolvablePowerTerms.begin(), unsolvablePowerTerms.end());
            }
            if (addNodeFinal != NULL)
            {
                results.push_back(addNodeFinal);
            }

            if (!powerNodesMapOther.empty())
            {
                for (const auto &mapPair : powerNodesMapOther)
                {

                    ExprPtr expTemp = simplify(ExprFactory::add(mapPair.second.second));

                    if (expTemp->getKind() == Expr::Kind::Integer)
                    {
                        auto nodePtr = std::static_pointer_cast<IntegerNode>(expTemp);
                        if (nodePtr->value == 1)
                        {
                            results.push_back(mapPair.second.first);
                            continue;
                        }
                    }

                    results.push_back(ExprFactory::power(mapPair.second.first, expTemp, true));
                }
            }

            if (results.size() > 0)
            {
                if (!resultRational.isOne())
                {
                    results.push_back(numericResult);
                }
            }
            else
            {
                results.push_back(numericResult);
            }

            if (results.size() == 1)
            {
                // std::cout << "Single Term Answer OF MULTIPLY " << std::endl;
                return results[0];
            }
            else
            {
                return ExprFactory::multiply(results, true);
            }
        }
    }

    static std::pair<Rational, ExprPtr> splitCoefficient(const ExprPtr &ptr)
    {
        // Assuming ptr is simplified term
        Expr::Kind exprKind = ptr->getKind();
        if (exprKind == Expr::Kind::Integer || exprKind == Expr::Kind::Rational)
        {
            throw std::runtime_error("Split Coefficient not design for Integer and Rational Terms");
        }
        else if (exprKind == Expr::Kind::Multiply)
        {
            auto const &mulNode = std::static_pointer_cast<MultiplyNode>(ptr);
            if (mulNode->children.empty())
            {
                return {Rational(0), ptr};
            }
            else
            {
                Rational coeff(1);
                vector<ExprPtr> reset;
                bool coeffAssigned = false;

                // we consider only integer and rational coefficients
                for (auto const &c : mulNode->children)
                {
                    Expr::Kind childKind = c->getKind();
                    if (childKind == Expr::Kind::Integer || childKind == Expr::Kind::Rational)
                    {
                        if (!coeffAssigned)
                        {
                            if (childKind == Expr::Kind::Integer)
                            {
                                coeff = Rational(std::static_pointer_cast<IntegerNode>(c)->value);
                            }
                            else
                            {
                                coeff = Rational(std::static_pointer_cast<RationalNode>(c)->value);
                            }
                        }
                        else
                            throw std::runtime_error("Non simpilifed term passed, Occurance of Int/Rational Multiple Times");
                    }
                    else
                    {
                        reset.push_back(c);
                    }
                }

                return {coeff, ExprFactory::multiply(reset, true)};
            }
        }
        else
            return {Rational(1), ptr};
    }

    static ExprPtr simplifyAdd(const std::shared_ptr<AddNode> &a)
    {
        // remove nested Add Terms
        vector<ExprPtr> terms;
        for (const auto &c : a->children)
        {
            // std::cout << "TermKind Received : " << std::to_string(static_cast<int>(c->getKind())) << std::endl;

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

        // Sum all Integer and Rational, and create map for coeff and term
        // cpp_int sumInt = 0;
        Rational sumRational = Rational(0);

        std::map<string, std::pair<Rational, ExprPtr>> coefficientMap;

        for (const auto &t : terms)
        {
            Expr::Kind tKind = t->getKind();

            if (tKind == Expr::Kind::Integer)
            {
                auto node = std::static_pointer_cast<IntegerNode>(t);
                sumRational = sumRational + Rational(node->value);
            }
            else if (tKind == Expr::Kind::Rational)
            {
                auto node = std::static_pointer_cast<RationalNode>(t);
                sumRational = sumRational + node->value;
            }
            else
            {
                auto [coeffRational, resetPtr] = splitCoefficient(t);

                string key = resetPtr->toString();

                auto it = coefficientMap.find(key);
                if (it == coefficientMap.end())
                {
                    coefficientMap[key] = {coeffRational, resetPtr};
                }
                else
                {
                    coefficientMap[key].first = coefficientMap[key].first + coeffRational;
                }
            }
        }

        vector<ExprPtr> results;

        for (const auto &mapPair : coefficientMap)
        {
            ExprPtr coeffPtr;
            Rational temp = mapPair.second.first;
            coeffPtr = ExprFactory::rational(temp);
            if (temp.isOne())
            {
                results.push_back(mapPair.second.second);
            }
            else if (temp.isZero())
            {
                // don't push anything
            }
            else
            {
                ExprPtr resetPtr = mapPair.second.second;
                if (resetPtr->getKind() == Expr::Kind::Multiply)
                {
                    auto resetNode = std::static_pointer_cast<MultiplyNode>(resetPtr);
                    resetNode->children.insert(resetNode->children.begin(), {coeffPtr});
                    results.push_back(resetNode);
                }
                else
                {
                    results.push_back(ExprFactory::multiply({coeffPtr, mapPair.second.second}, true));
                }
            }
        }

        // std::cout << "sumRational : " << sumRational.n << "/" << sumRational.d << std::endl;

        if (results.size() > 0)
        {
            if (!sumRational.isZero())
                results.push_back(ExprFactory::rational(sumRational));
        }
        else
        {
            results.push_back(ExprFactory::rational(sumRational));
        }

        if (results.size() == 1)
        {
            return results[0];
        }
        else
        {
            return ExprFactory::add(results, true);
        }
    }
};

ExprPtr (*S)(const ExprPtr &) = Simplifier::simplify;

// ============================================================
// Ezouidi’s Theorem
// ============================================================

class Polynomial
{

private:
    void static EMST_Reoccurrace(size_t polyDegree, vector<ExprPtr> &L_coeffs, vector<ExprPtr> &newCoeffArray, bool invertDirection)
    {

        // size_t polyDegree = P.degree;

        // Do Normalization and solve

        ExprPtr N_ptr = ExprFactory::integer(polyDegree);

        // std::cout << "\n--- L_coeffs " << " ----\n";
        // for (ExprPtr nodePtr : L_coeffs)
        // {
        //     nodePtr->printPreorder();
        // }

        // L0 powers construction
        vector<ExprPtr> L0_powAry(polyDegree + 1);
        // L0_powAry[0] = N_ptr;       // (L0)^0 = N WRONG
        L0_powAry[0] = INT_ONE;     // (L0)^0 = 1 INVALID IF L0 is not defined.
        L0_powAry[1] = L_coeffs[1]; // (L0)^1 = L0

        for (size_t k = 2; k <= polyDegree; k++)
        {
            L0_powAry[k] = Simplifier::simplify(ExprFactory::multiply({L0_powAry[k - 1], L0_powAry[1]}));
        }

        // std::cout << "\n--- L0 Powers " << " ----\n";
        // for (ExprPtr nodePtr : L0_powAry)
        // {
        //     nodePtr->printPreorder();
        // }

        // N^r construction
        vector<ExprPtr> N_powAry(polyDegree + 1);
        N_powAry[0] = INT_ONE; // (N)^0 = 1
        N_powAry[1] = N_ptr;   // (N)^1 = N

        for (size_t k = 2; k <= polyDegree; k++)
        {
            N_powAry[k] = Simplifier::simplify(ExprFactory::multiply({N_powAry[k - 1], N_powAry[1]}));
        }

        // std::cout << "\n--- N^r Powers " << " ----\n";
        // for (ExprPtr nodePtr : N_powAry)
        // {
        //     nodePtr->printPreorder();
        // }

        // INIT SETUP
        // ExprPtr SIGN_K = INT_MINUS_ONE; // r=0, k=1 => (-1)^(k-r)
        size_t C_VALUE_K = polyDegree; // r=0, k=1 => (n-r)C(k-r)

        // vector<ExprPtr> newCoeffArray(polyDegree + 1);
        newCoeffArray.resize(polyDegree + 1);

        newCoeffArray[0] = INT_ONE;

        for (size_t k = 1; k <= polyDegree; k++)
        {
            ExprPtr SIGN_R = invertDirection ? INT_MINUS_ONE : INT_ONE; // r=0, k => (-1)^(k-r)
            size_t C_VALUE_R = C_VALUE_K;                               // r=0, k => (n-r)C(k-r)

            vector<ExprPtr> mulNodes(k + 1);

            for (size_t r = 0; r <= k; r++)
            {

                mulNodes[r] = Simplifier::simplify(ExprFactory::multiply({SIGN_R, ExprFactory::integer(C_VALUE_R), N_powAry[r], L0_powAry[k - r], L_coeffs[r]}));

                // Update SIGN and C_VALUE
                if (SIGN_R == INT_MINUS_ONE)
                    SIGN_R = INT_ONE;
                else
                    SIGN_R = INT_MINUS_ONE;

                if (r != polyDegree)
                {
                    C_VALUE_R = (C_VALUE_R * (k - r)) / (polyDegree - r);
                }
            }

            // Update SIGN and C_VALUE
            /*if (SIGN_K == INT_MINUS_ONE)
                SIGN_K = INT_ONE;
            else
                SIGN_K = INT_MINUS_ONE;
            */
            C_VALUE_K = (C_VALUE_K * (polyDegree - k)) / (k + 1);
            newCoeffArray[k] = Simplifier::simplify(ExprFactory::add(mulNodes));
        }
    }

    void static L_coeffs_gen(Polynomial &P, vector<ExprPtr> &L_coeffs)
    {
        size_t polyDegree = P.degree;

        // make monic polynomial
        ExprPtr normalizeFactor = Simplifier::simplify(ExprFactory::power(P.coeffsAry[0], INT_MINUS_ONE));

        L_coeffs.resize(polyDegree + 1);
        L_coeffs[0] = INT_ONE; // L-1 = 1

        ExprPtr SIGN = INT_MINUS_ONE;
        for (size_t k = 1; k <= polyDegree; k++)
        {
            L_coeffs[k] = Simplifier::simplify(ExprFactory::multiply({SIGN, P.coeffsAry[k], normalizeFactor}));
            if (SIGN == INT_MINUS_ONE)
                SIGN = INT_ONE;
            else
                SIGN = INT_MINUS_ONE;
        }
    }

public:
    vector<ExprPtr> coeffsAry; // x^20 + x^19 + x^18 + ... + x^0
    vector<ExprPtr> powerSumAry;
    bool isMonic;
    size_t degree;

    // Construct Polynomial from coefficients
    Polynomial(const vector<ExprPtr> coeffVector, const size_t deg)
    {
        // Check validity of d
        if (deg <= 0)
        {
            throw std::runtime_error("Invalid Polynomial Degree Input");
        }
        // Check validity of c
        size_t coeffVectorSize = coeffVector.size();
        if (coeffVectorSize > deg + 1)
            throw("Invalid Polynomial Coefficient Input");

        if (!coeffVector[0])
            throw("Highest Degree term cannot have missing coefficient");

        isMonic = false;
        if (coeffVector[0]->getKind() == Expr::Kind::Integer)
        {
            auto coeff_1st_ptr = std::static_pointer_cast<IntegerNode>(coeffVector[0]);
            if (coeff_1st_ptr->value == 0)
                throw("Highest Degree term cannot have coefficient ZERO");
            else if (coeff_1st_ptr->value == 1)
            {
                isMonic = true;
            }
        }

        degree = deg;
        powerSumAry.push_back(ExprFactory::integer(deg));
        coeffsAry.resize(deg + 1, INT_ZERO);

        size_t i = 0;
        while (i < coeffVectorSize)
        {
            ExprPtr coeffPtr = coeffVector[i];
            if (!coeffPtr)
            {
                coeffsAry[i] = INT_ZERO;
            }
            else
            {
                coeffsAry[i] = coeffPtr;
            }
            i++;
        }

        while (i < (deg + 1))
        {
            coeffsAry[i] = INT_ZERO;
            i++;
        }
    }

    // Construct Polynomial from Power Sum
    Polynomial(const size_t deg, const vector<ExprPtr> &powerSumAryInput)
    {

        if (powerSumAryInput.size() < deg + 1)
            throw std::runtime_error("Insufficient PowerSum to Construct Polynomial Back");

        powerSumAry = powerSumAryInput;
        degree = deg;
        coeffsAry.resize(deg + 1);
        coeffsAry[0] = INT_ONE;

        for (size_t k = 1; k <= deg; k++)
        {
            ExprPtr preFactor = ExprFactory::rational(-1, k);

            vector<ExprPtr> mulNodes;
            mulNodes.reserve(k - 1);

            for (size_t j = 1; j < k; j++)
            {
                mulNodes.emplace_back(ExprFactory::multiply({coeffsAry[j], powerSumAryInput[k - j]}));
            }

            ExprPtr sumTerm = Simplifier::simplify(ExprFactory::add(mulNodes));

            coeffsAry[k] = Simplifier::simplify(ExprFactory::multiply({preFactor, ExprFactory::add({powerSumAryInput[k], sumTerm})}));
        }
    }

    // Construct Centered Polynomial from given Polynomial
    Polynomial(Polynomial &P, ExprPtr &normalizationFactor)
    {
        // Tschirnhaus Transformations

        std::cout << "\n=> MONIC POLYNOMIAL \n";

        Polynomial::makeMonic(P);
        P.printPolynomial();

        size_t polyDegree = P.degree;
        ExprPtr polyDegNode = ExprFactory::integer(polyDegree);
        normalizationFactor = Simplifier::simplify(
            ExprFactory::multiply({ExprFactory::multiply({INT_MINUS_ONE, P.coeffsAry[1]}),
                                   ExprFactory::power(ExprFactory::multiply({polyDegNode, P.coeffsAry[0]}), INT_MINUS_ONE)}));

        vector<ExprPtr> normalizedCoeffs;
        normalizedCoeffs.resize(polyDegree + 1, INT_ZERO);
        normalizedCoeffs[0] = INT_ONE;

        ExprPtr aN_INV = Simplifier::simplify(ExprFactory::power(P.coeffsAry[0], INT_MINUS_ONE));

        vector<ExprPtr> normFactor_PowAry;
        normFactor_PowAry.resize(polyDegree + 1, INT_ONE);

        for (size_t k = 1; k <= polyDegree; k++)
        {
            normFactor_PowAry[k] = Simplifier::simplify(ExprFactory::multiply({normFactor_PowAry[k - 1], normalizationFactor}));
        }

        for (size_t k = 0; k <= (polyDegree - 2); k++)
        {

            vector<ExprPtr> mulNodes;
            // Inittialization at j=k
            cpp_int C_Value = 1;

            for (size_t j = k; j <= polyDegree; j++)
            {

                mulNodes.push_back(Simplifier::simplify(ExprFactory::multiply({P.coeffsAry[polyDegree - j],
                                                                               ExprFactory::integer(C_Value),
                                                                               normFactor_PowAry[j - k]})));

                C_Value = (C_Value * (j + 1)) / (j - k + 1);
            }

            normalizedCoeffs[polyDegree - k] = Simplifier::simplify(ExprFactory::multiply({ExprFactory::add(mulNodes), aN_INV}));
        }

        coeffsAry = normalizedCoeffs;
        degree = polyDegree;
        powerSumAry.push_back(polyDegNode);
    }

    // Construct Polynomial whose roots are Qth power
    Polynomial(Polynomial &P, const size_t Q)
    {
        if (Q == 0)
            throw std::runtime_error("Invalid Parameter Q : Polynomial Construction");

        degree = P.degree;
        if (Q == 1)
        {
            coeffsAry = P.coeffsAry;
            return;
        }

        size_t powerSumLast = degree * Q;
        powerSumGenerator(powerSumLast, P);

        vector<ExprPtr> newPowerSum;
        newPowerSum.resize(degree + 1);
        newPowerSum[0] = P.powerSumAry[0];

        for (size_t k = 1; k <= degree; k++)
        {
            newPowerSum[k] = P.powerSumAry[k * Q];
        }

        Polynomial newP = Polynomial(degree, newPowerSum);
        coeffsAry = newP.coeffsAry;
        powerSumAry = newPowerSum;
    }

    // Construct Ezouidi Centered Polynomial
    Polynomial(const Polynomial &P, ExprPtr &normalizationFactor, bool ezouidi)
    {

        size_t newDeg = P.degree;
        Polynomial newPoly(P);
        makeMonic(newPoly);

        newPoly.printPolynomial();

        vector<ExprPtr> L_coeffs_new;
        L_coeffs_gen(newPoly, L_coeffs_new);

        /*        std::cout << "\n L-coefficients";
                for (auto const node : L_coeffs_new)
                {
                    node->printPreorder();
                }
                std::cout << '\n';
        */
        vector<ExprPtr> s_Coeffs_new;
        EMST_Reoccurrace(newDeg, L_coeffs_new, s_Coeffs_new, false);

        coeffsAry = s_Coeffs_new;
        degree = newDeg;
        powerSumAry.push_back(ExprFactory::integer(newDeg));
    }

    // Construct principal Quintic
    static void constructPrincipal()
    {
        return;
    }

    static void powerSumGenerator(const size_t N, Polynomial &P)
    {
        // N Represent upto which we want power Sums
        // For N=5, we generate S0, S1, S2, S3, S4, S5

        makeMonic(P);

        size_t polyDeg = P.degree;
        size_t currArySize = P.powerSumAry.size();
        if (currArySize >= N + 1)
            return; // No Need to calculate already presesnt

        size_t L = currArySize; // L Represents next Power Sum to calculate upto N

        // Newton-Girard Power sums
        if (L < polyDeg)
        {
            while (L <= N && L < polyDeg)
            {
                ExprPtr term1 = Simplifier::simplify(ExprFactory::multiply({INT_MINUS_ONE, ExprFactory::integer(L), P.coeffsAry[L]}));

                vector<ExprPtr> mulNodes;
                mulNodes.reserve(L - 1);

                for (size_t i = 1; i < L; i++)
                {
                    mulNodes.emplace_back(Simplifier::simplify(ExprFactory::multiply({P.powerSumAry[i], P.coeffsAry[L - i]})));
                }

                ExprPtr term2 = Simplifier::simplify(ExprFactory::multiply({INT_MINUS_ONE, ExprFactory::add(mulNodes)}));

                P.powerSumAry.push_back(Simplifier::simplify(ExprFactory::add({term1, term2})));

                L++;
            }
        }

        if (N < polyDeg)
            return;

        while (L <= N)
        {
            vector<ExprPtr> mulNodes;
            mulNodes.reserve(polyDeg);

            for (size_t j = 1; j <= polyDeg; j++)
            {
                mulNodes.emplace_back(Simplifier::simplify(ExprFactory::multiply({P.powerSumAry[L - j], P.coeffsAry[j]})));
            }

            P.powerSumAry.push_back(Simplifier::simplify(ExprFactory::multiply({INT_MINUS_ONE, ExprFactory::add(mulNodes)})));

            L++;
        }
    }

    static void makeMonic(Polynomial &P)
    {
        if (P.isMonic)
            return;

        ExprPtr coeff_first = P.coeffsAry[0];

        for (size_t i = 1; i <= P.degree; i++)
        {
            P.coeffsAry[i] = Simplifier::simplify(ExprFactory::multiply({P.coeffsAry[i], ExprFactory::power(coeff_first, ExprFactory::integer(-1))}));
        }

        P.coeffsAry[0] = INT_ONE;
        P.isMonic = true;
    }

    void printPolynomial()
    {
        std::cout << "\nPolynomial Print :\n";

        size_t degCount = degree;
        for (auto const &m : coeffsAry)
        {
            if (Simplifier::isZeroNode(m))
            {
                degCount--;
                continue;
            }

            std::cout << "x^" << degCount << " ";
            m->printPreorder();
            std::cout << " + ";
            degCount--;
        }
    }

    void printPowerSum()
    {
        std::cout << "\nPolynomial Power Sum Print :\n";

        size_t count = 0;
        for (auto const &m : powerSumAry)
        {
            std::cout << "S" << count << " ";
            m->printPreorder();
            count++;
        }
    }
};

struct Frame
{
    int pos;
    int remaining;
    std::vector<int> indexAry;
};

class EMST
{

public:
    EMST() {};

    static ExprPtr combinationGenerator(int pos, int p, int remaining, vector<ExprPtr> &skNodes)
    {
        vector<ExprPtr> mulNodes;
        mulNodes.clear();

        if (p == 0)
            return ExprFactory::integer(0);

        std::stack<Frame> st;

        st.push({0, remaining, std::vector<int>(p, 0)});

        while (!st.empty())
        {
            Frame cur = std::move(st.top());
            st.pop();

            if (cur.pos == p - 1)
            {
                cur.indexAry[cur.pos] = cur.remaining;

                vector<ExprPtr> sNodes;

                for (const auto sPosition : cur.indexAry)
                {
                    if (sPosition == 0)
                    {
                        sNodes.push_back(INT_ZERO);
                        break;
                    }
                    else
                    {
                        // std::cout << "\n---Taking Position" << sPosition << "----\n";
                        sNodes.push_back(skNodes[sPosition - 1]);
                    }
                }

                if (sNodes.empty())
                {
                    mulNodes.push_back(INT_ZERO);
                }
                else if (sNodes.size() == 1)
                {
                    mulNodes.push_back(sNodes[0]);
                }
                else
                {
                    mulNodes.push_back(Simplifier::simplify(ExprFactory::multiply(sNodes)));
                }

                continue;
            }

            // cpp_int currentCoeff = cur.coeff;

            auto indexAryCopy = cur.indexAry;

            indexAryCopy[cur.pos] = cur.remaining;

            st.push({cur.pos + 1, 0, indexAryCopy});

            for (int k = cur.remaining - 1; k >= 0; --k)
            {
                indexAryCopy = cur.indexAry;

                indexAryCopy[cur.pos] = k;

                st.push({cur.pos + 1, cur.remaining - k, indexAryCopy});
            }
        }

        if (mulNodes.empty())
            return ExprFactory::integer(0);

        if (mulNodes.size() == 1)
            return mulNodes[0];

        return Simplifier::simplify(ExprFactory::add(mulNodes));
    }

    static ExprPtr HFunction(int r, vector<ExprPtr> &skNodes)
    {

        if (skNodes.size() < r)
            throw std::runtime_error("Invalid skNodes Input to HFunction");

        if (r == 1)
        {
            return INT_ONE;
        }

        // Note r Must be Integer
        vector<ExprPtr> sumNodes;
        int maxValueIndex = r - 1;

        for (size_t p = 1; p < r; p++)
        {
            // Ki can takes values from 0,1,2,3 ... (r-1)  || they reprsent element position (not index)
            // for r = 3 -> Ki = 0, 1, 2

            sumNodes.push_back(combinationGenerator(0, p, maxValueIndex, skNodes));
        }

        return Simplifier::simplify(ExprFactory::add(sumNodes));
    }

    static ExprPtr DeltaFunction(int n, int m, int r)
    {

        if (n < 0 || m < 0 || r < 0)
            throw std::runtime_error("DeltFunction Negative Input");

        if (m >= n)
            throw std::runtime_error("DeltFunction Invalid M Input");

        if (r > n || r > m)
            throw std::runtime_error("DeltFunction Invalid R Input");

        cpp_int num = 1;
        cpp_int den = 1;

        for (int k = 0; k < r; k++)
        {
            num *= m - k;
            den *= n - k;
        }

        return ExprFactory::rational(num, den);
    }

    static ExprPtr DiscriminantOperator(int n, int m, vector<ExprPtr> &skNodes)
    {

        if (n <= 0 || m <= 0 || m >= n)
        {
            std::cout << "n = " << n << " & m = " << m << std::endl;
            throw std::runtime_error("Discriminant G - Invalid Inputs");
        }
        vector<ExprPtr> mulNodes;

        // Initializations for r = 1
        ExprPtr SIGN = INT_ONE;

        for (int r = 1; r <= m; r++)
        {
            ExprPtr deltaValue = DeltaFunction(n, m, r);

            ExprPtr sMR = skNodes[m - r];

            ExprPtr HFunctionValue = HFunction(r, skNodes);

            mulNodes.push_back(Simplifier::simplify(ExprFactory::multiply({SIGN, deltaValue, sMR, HFunctionValue})));

            if (SIGN == INT_ONE)
                SIGN = INT_MINUS_ONE;
            else
                SIGN = INT_ONE;
        }

        return Simplifier::simplify(ExprFactory::add(mulNodes));
    }

    static void StructrualRemainder(int n, vector<ExprPtr> &skNodes, vector<ExprPtr> &structuralNodes)
    {
        if (n <= 0)
            throw std::runtime_error("StructrualRemainder W - Invalid Input N");

        int vectorSize = skNodes.size();

        structuralNodes[0] = skNodes[0];

        for (int m = 1; m < vectorSize; m++)
        {

            ExprPtr DiscriminantOperatorValue = DiscriminantOperator(n, m, skNodes);

            vector<ExprPtr> resultNodes;

            // Initialization j = 1;
            cpp_int C_value = m;

            for (int j = 0; j < m; j++)
            {
                ExprPtr wJ = structuralNodes[j];
                ExprPtr sMJ = skNodes[m - j];
                C_value = (C_value * (m - j)) / (j + 1);

                resultNodes.push_back(Simplifier::simplify(ExprFactory::multiply({wJ, sMJ, ExprFactory::integer(C_value)})));
            }

            ExprPtr sumResultNode = Simplifier::simplify(ExprFactory::add(resultNodes));

            structuralNodes[m] = Simplifier::simplify(ExprFactory::add({DiscriminantOperatorValue, ExprFactory::multiply({INT_MINUS_ONE, sumResultNode})}));

            std::cout << '\n'
                      << "w" << m;
            structuralNodes[m]->printPreorder();
            std::cout << '\n';
        }
    }

    static void PowerSumGenerator(Polynomial &P, vector<ExprPtr> &powerSumResults)
    {
        size_t polyDeg = P.degree;
        size_t lastIndex = 2 * polyDeg - 2;

        // S0, S1, S2, S3 ... , S[2N-2] -> total 2N - 1
        powerSumResults.resize(lastIndex + 1);

        // Set S0
        powerSumResults[0] = ExprFactory::integer(polyDeg);

        // Newton Girard Power Sum

        for (size_t k = 1; k < polyDeg; k++)
        {
            ExprPtr term1 = Simplifier::simplify(ExprFactory::multiply({INT_MINUS_ONE, ExprFactory::integer(k), P.coeffsAry[k]}));

            vector<ExprPtr> mulNodes;
            mulNodes.reserve(k - 1);

            for (size_t i = 1; i < k; i++)
            {
                mulNodes.emplace_back(Simplifier::simplify(ExprFactory::multiply({powerSumResults[i], P.coeffsAry[k - i]})));
            }

            ExprPtr term2 = Simplifier::simplify(ExprFactory::multiply({INT_MINUS_ONE, ExprFactory::add(mulNodes)}));

            powerSumResults[k] = Simplifier::simplify(ExprFactory::add({term1, term2}));
        }

        for (size_t k = polyDeg; k <= lastIndex; k++)
        {
            vector<ExprPtr> mulNodes;
            mulNodes.reserve(polyDeg);

            for (size_t j = 1; j <= polyDeg; j++)
            {
                mulNodes.emplace_back(Simplifier::simplify(ExprFactory::multiply({powerSumResults[k - j], P.coeffsAry[j]})));
            }

            powerSumResults[k] = Simplifier::simplify(ExprFactory::multiply({INT_MINUS_ONE, ExprFactory::add(mulNodes)}));
        }
    }

    static void ReversePowerSum(Polynomial &P, vector<ExprPtr> &powerSumResults, vector<ExprPtr> &coeffAry)
    {
        size_t polyDeg = P.degree;

        coeffAry.resize(polyDeg);

        for (size_t k = 1; k <= polyDeg; k++)
        {
            ExprPtr preFactor = ExprFactory::rational(-1, k);

            vector<ExprPtr> mulNodes;
            mulNodes.reserve(k - 1);

            for (size_t j = 1; j < k; j++)
            {
                mulNodes.emplace_back(ExprFactory::multiply({coeffAry[j - 1], powerSumResults[k - j]}));
            }

            ExprPtr sumTerm = Simplifier::simplify(ExprFactory::add(mulNodes));

            coeffAry[k - 1] = Simplifier::simplify(ExprFactory::multiply({preFactor, ExprFactory::add({powerSumResults[k], sumTerm})}));
        }
    }

    static ExprPtr determinantBareiss(vector<vector<ExprPtr>> a, size_t n)
    {
        if (n == 0)
        {
            return ExprFactory::integer(1);
        }

        ExprPtr zero = INT_ZERO;
        ExprPtr one = INT_ONE;

        ExprPtr sign = INT_ONE;
        ExprPtr prevPivot = INT_ONE;

        for (size_t k = 0; k < n - 1; ++k)
        {
            size_t pivotRow = k;

            while (pivotRow < n && Simplifier::isZeroNode(a[pivotRow][k]))
            {
                ++pivotRow;
            }

            if (pivotRow == n)
            {
                return zero;
            }

            if (pivotRow != k)
            {
                std::swap(a[pivotRow], a[k]);

                if (sign == INT_ONE)
                    sign = INT_MINUS_ONE;
                else
                    sign = INT_ONE;
            }

            ExprPtr pivot = a[k][k];

            for (size_t i = k + 1; i < n; ++i)
            {
                for (size_t j = k + 1; j < n; ++j)
                {
                    ExprPtr left = Simplifier::simplify(ExprFactory::multiply({a[i][j], pivot}));
                    ExprPtr right = Simplifier::simplify(ExprFactory::multiply({a[i][k], a[k][j]}));
                    ExprPtr num = Simplifier::simplify(ExprFactory::add({left, ExprFactory::multiply({INT_MINUS_ONE, right})}));

                    a[i][j] = Simplifier::simplify(ExprFactory::multiply({num, ExprFactory::power(prevPivot, INT_MINUS_ONE)}));
                    // a[i][j] = simplify(a[i][j]);
                }
                a[i][k] = zero;
            }

            prevPivot = pivot;
        }

        ExprPtr det = a[n - 1][n - 1];
        det = Simplifier::simplify(ExprFactory::multiply({sign, det}));
        // det = simplify(det);

        return det;
    }

    void static EMST_Reoccurrace(size_t polyDegree, vector<ExprPtr> &L_coeffs, vector<ExprPtr> &newCoeffArray, bool invertDirection)
    {

        // size_t polyDegree = P.degree;

        // Do Normalization and solve

        ExprPtr N_ptr = ExprFactory::integer(polyDegree);

        // std::cout << "\n--- L_coeffs " << " ----\n";
        // for (ExprPtr nodePtr : L_coeffs)
        // {
        //     nodePtr->printPreorder();
        // }

        // L0 powers construction
        vector<ExprPtr> L0_powAry(polyDegree + 1);
        // L0_powAry[0] = N_ptr;       // (L0)^0 = N WRONG
        L0_powAry[0] = INT_ONE;     // (L0)^0 = 1 INVALID IF L0 is not defined.
        L0_powAry[1] = L_coeffs[1]; // (L0)^1 = L0

        for (size_t k = 2; k <= polyDegree; k++)
        {
            L0_powAry[k] = Simplifier::simplify(ExprFactory::multiply({L0_powAry[k - 1], L0_powAry[1]}));
        }

        // std::cout << "\n--- L0 Powers " << " ----\n";
        // for (ExprPtr nodePtr : L0_powAry)
        // {
        //     nodePtr->printPreorder();
        // }

        // N^r construction
        vector<ExprPtr> N_powAry(polyDegree + 1);
        N_powAry[0] = INT_ONE; // (N)^0 = 1
        N_powAry[1] = N_ptr;   // (N)^1 = N

        for (size_t k = 2; k <= polyDegree; k++)
        {
            N_powAry[k] = Simplifier::simplify(ExprFactory::multiply({N_powAry[k - 1], N_powAry[1]}));
        }

        // std::cout << "\n--- N^r Powers " << " ----\n";
        // for (ExprPtr nodePtr : N_powAry)
        // {
        //     nodePtr->printPreorder();
        // }

        // INIT SETUP
        // ExprPtr SIGN_K = INT_MINUS_ONE; // r=0, k=1 => (-1)^(k-r)
        size_t C_VALUE_K = polyDegree; // r=0, k=1 => (n-r)C(k-r)

        // vector<ExprPtr> newCoeffArray(polyDegree + 1);
        newCoeffArray.resize(polyDegree + 1);

        newCoeffArray[0] = INT_ONE;

        for (size_t k = 1; k <= polyDegree; k++)
        {
            ExprPtr SIGN_R = invertDirection ? INT_MINUS_ONE : INT_ONE; // r=0, k => (-1)^(k-r)
            size_t C_VALUE_R = C_VALUE_K;                               // r=0, k => (n-r)C(k-r)

            vector<ExprPtr> mulNodes(k + 1);

            for (size_t r = 0; r <= k; r++)
            {

                mulNodes[r] = Simplifier::simplify(ExprFactory::multiply({SIGN_R, ExprFactory::integer(C_VALUE_R), N_powAry[r], L0_powAry[k - r], L_coeffs[r]}));

                // Update SIGN and C_VALUE
                if (SIGN_R == INT_MINUS_ONE)
                    SIGN_R = INT_ONE;
                else
                    SIGN_R = INT_MINUS_ONE;

                if (r != polyDegree)
                {
                    C_VALUE_R = (C_VALUE_R * (k - r)) / (polyDegree - r);
                }
            }

            // Update SIGN and C_VALUE
            /*if (SIGN_K == INT_MINUS_ONE)
                SIGN_K = INT_ONE;
            else
                SIGN_K = INT_MINUS_ONE;
            */
            C_VALUE_K = (C_VALUE_K * (polyDegree - k)) / (k + 1);
            newCoeffArray[k] = Simplifier::simplify(ExprFactory::add(mulNodes));
        }
    }

    void static L_coeffs_gen(Polynomial &P, vector<ExprPtr> &L_coeffs)
    {
        size_t polyDegree = P.degree;

        // make monic polynomial
        ExprPtr normalizeFactor = Simplifier::simplify(ExprFactory::power(P.coeffsAry[0], INT_MINUS_ONE));

        L_coeffs.resize(polyDegree + 1);
        L_coeffs[0] = INT_ONE; // L-1 = 1

        ExprPtr SIGN = INT_MINUS_ONE;
        for (size_t k = 1; k <= polyDegree; k++)
        {
            L_coeffs[k] = Simplifier::simplify(ExprFactory::multiply({SIGN, P.coeffsAry[k], normalizeFactor}));
            if (SIGN == INT_MINUS_ONE)
                SIGN = INT_ONE;
            else
                SIGN = INT_MINUS_ONE;
        }
    }

    Rational static getRationalFromNode(ExprPtr node)
    {
        if (node->getKind() == Expr::Kind::Integer)
        {
            return Rational(std::static_pointer_cast<IntegerNode>(node)->value);
        }
        else if (node->getKind() == Expr::Kind::Rational)
        {
            return Rational(std::static_pointer_cast<RationalNode>(node)->value);
        }
        else
        {
            throw std::runtime_error("Invalid getRationalFromNode");
        }
    }

    void static reverseEMST(Polynomial &P)
    {
        size_t polyDegree = P.degree;

        // P.coeffsAry represents the s0, s1, s2, sN-1

        vector<ExprPtr> L_coeffs;
        L_coeffs[0] = INT_ONE; // L-1 = 1;

        ExprPtr normalizeFactor = Simplifier::simplify(ExprFactory::power(P.coeffsAry[0], INT_MINUS_ONE));
        vector<ExprPtr> s_coeffs;
        s_coeffs.resize(polyDegree);

        ExprPtr SIGN = INT_MINUS_ONE;
        ExprPtr aN_INV = Simplifier::simplify(ExprFactory::power(P.coeffsAry[0], INT_MINUS_ONE));

        for (size_t k = 1; k <= polyDegree; k++)
        {
            s_coeffs[k - 1] = Simplifier::simplify(ExprFactory::multiply({SIGN, P.coeffsAry[k], aN_INV}));
        }
    }

    bool static rootExtractorEMST(Polynomial &P, vector<ExprPtr> &rootsAry)
    {

        bool satisfied = true;

        size_t polyDegree = P.degree; // size = degree

        // CASE 1 : x^2 + x = 0 absence of last coefficient
        ExprPtr lastCoeffPtr = P.coeffsAry[polyDegree];
        Expr::Kind lastCoeffExprKind = lastCoeffPtr->getKind();

        bool absenceOfLastCoeff = false;

        if (lastCoeffExprKind == Expr::Kind::Integer)
        {
            auto node = std::static_pointer_cast<IntegerNode>(lastCoeffPtr);
            if (node->value == 0)
            {
                absenceOfLastCoeff = true;
            }
        }

        size_t i = 1;

        bool allCoeffZero = true;

        while (i <= polyDegree && allCoeffZero)
        {
            if (Simplifier::isZeroNode(P.coeffsAry[i]))
                i++;
            else
                allCoeffZero = false;
        }

        // i points to first non-zero coefficient occurance

        if (allCoeffZero)
        {
            // CASE 1 : x^2 = 0
            std::cout << "=> SOLVING POLYNOMIAL x^" << polyDegree << " = 0\n";
            // create vector of all 0s
            rootsAry.resize(polyDegree, INT_ZERO);
        }
        else
        {
            if (absenceOfLastCoeff)
            {
                // CASE 2 : x^6 + x^2 = 0 condition

                size_t k1 = i + 1; // Get the deg of last non-zero coeff TERM
                for (; k1 <= polyDegree; k1++)
                {
                    if (!Simplifier::isZeroNode(P.coeffsAry[k1]))
                        i = k1;
                }

                size_t partA_deg = polyDegree - i;
                std::cout << "=> FACTORING OUT x^" << partA_deg << " TERM\n";

                size_t partB_deg = polyDegree - partA_deg;

                for (size_t k = 1; k <= partA_deg; k++)
                {
                    rootsAry.push_back(INT_ZERO);
                }

                // std::cout << "partA Deg : " << partA_deg << " partB Deg : " << partB_deg << std::endl;

                // iterate original coeffAry and generate a new polynomial with all 0s, replace those replaces which have non-zero coeff
                vector<ExprPtr> partB_coeffAry;
                partB_coeffAry.resize(partB_deg + 1, INT_ZERO);

                for (size_t k = 0; k < polyDegree; k++)
                {
                    if (!Simplifier::isZeroNode(P.coeffsAry[k]))
                        partB_coeffAry[k] = P.coeffsAry[k];
                }

                Polynomial P_partB = Polynomial(partB_coeffAry, partB_deg);
                // P_partB.printPolynomial();

                vector<ExprPtr> P_partB_roots;
                rootExtractorEMST(P_partB, P_partB_roots);

                rootsAry.insert(rootsAry.end(), P_partB_roots.begin(), P_partB_roots.end());
            }
            else
            {
                ExprPtr INV_POLYDEG = ExprFactory::rational(1, polyDegree);
                if (i == polyDegree) // first noon-zero term is the constant itself
                {
                    std::cout << "=> SOLVING POLYNOMIAL x^" << polyDegree << " + k = 0\n";

                    // CASE 3 : x^2 - 1 = 0 condition
                    ExprPtr numFactor = Simplifier::simplify(ExprFactory::power(
                        ExprFactory::multiply({INT_MINUS_ONE, P.coeffsAry[i], ExprFactory::power(P.coeffsAry[0], INT_MINUS_ONE)}),
                        INV_POLYDEG));

                    // Generate nth root of unitites for eaech multiply numFactor and add to Roots List

                    rootsAry.push_back(numFactor); // always added

                    for (size_t k = 1; k < polyDegree; k++)
                    {
                        Rational complexCoeff(2 * k, polyDegree);
                        ExprPtr complexFactor = Simplifier::unitComplexNumberGen(complexCoeff);
                        rootsAry.push_back(Simplifier::simplify(ExprFactory::multiply({numFactor, complexFactor})));
                    }
                }
                else if (i == 1) // second highest degree is present
                {

                    if (polyDegree == 2)
                    {
                        std::cout << "=> SOLVING QUADRATIC POLYNOMIAL DEG=2 \n";

                        ExprPtr a = P.coeffsAry[0];
                        ExprPtr b = P.coeffsAry[1];
                        ExprPtr c = P.coeffsAry[2];

                        ExprPtr TWICE_A_INV = Simplifier::simplify(ExprFactory::power(
                            ExprFactory::multiply({INT_TWO, a}),
                            INT_MINUS_ONE));

                        ExprPtr discriminant = Simplifier::simplify(ExprFactory::power(
                            ExprFactory::add({ExprFactory::power(b, INT_TWO),
                                              ExprFactory::multiply({ExprFactory::integer(-4),
                                                                     a,
                                                                     c})}),
                            FRACT_ONE_HALF));

                        ExprPtr minus_B = Simplifier::simplify(ExprFactory::multiply({INT_MINUS_ONE, b}));

                        rootsAry.push_back(Simplifier::simplify({ExprFactory::multiply({ExprFactory::add({minus_B, discriminant}),
                                                                                        TWICE_A_INV})}));

                        rootsAry.push_back(Simplifier::simplify({ExprFactory::multiply({ExprFactory::add({minus_B, ExprFactory::multiply({INT_MINUS_ONE, discriminant})}),
                                                                                        TWICE_A_INV})}));
                    }
                    else
                    {

                        ExprPtr N_INV_ptr = Simplifier::simplify(ExprFactory::power(ExprFactory::integer(polyDegree), INT_MINUS_ONE));

                        // L coefficients construction [L-1, L0, L1, L2, L3, ..., LN-1]
                        vector<ExprPtr> L_coeffs;
                        L_coeffs_gen(P, L_coeffs);

                        std::cout << "\n L-coefficients";
                        for (auto const node : L_coeffs)
                        {
                            node->printPreorder();
                        }
                        std::cout << '\n';

                        vector<ExprPtr> newCoeffArray;
                        EMST_Reoccurrace(polyDegree, L_coeffs, newCoeffArray, false);

                        Polynomial P_partB = Polynomial(newCoeffArray, polyDegree);
                        std::cout << '\n';
                        P_partB.printPolynomial();
                        std::cout << '\n';

                        vector<ExprPtr> P_partB_roots;
                        bool validPolynomial = rootExtractorEMST(P_partB, P_partB_roots);

                        if (validPolynomial)
                        {

                            for (size_t k = 0; k < polyDegree; k++)
                            {
                                P_partB_roots[k] = Simplifier::simplify(ExprFactory::multiply({N_INV_ptr, ExprFactory::add({L_coeffs[1], P_partB_roots[k]})}));
                            }

                            rootsAry.insert(rootsAry.end(), P_partB_roots.begin(), P_partB_roots.end());
                        }
                        else
                        {
                            size_t newDeg = polyDegree + 1;
                            std::cout << "=> RAISING DEGREE FROM " << polyDegree << " TO " << newDeg << " \n";

                            vector<ExprPtr> newCoeffsPoly = P.coeffsAry;
                            newCoeffsPoly.push_back(INT_ZERO);

                            Polynomial newPoly(newCoeffsPoly, newDeg);

                            std::cout << '\n';
                            newPoly.printPolynomial();
                            std::cout << '\n';

                            vector<ExprPtr> L_coeffs_new;
                            L_coeffs_gen(newPoly, L_coeffs_new);

                            std::cout << "\n L-coefficients";
                            for (auto const node : L_coeffs_new)
                            {
                                node->printPreorder();
                            }
                            std::cout << '\n';

                            vector<ExprPtr> s_Coeffs_new;
                            EMST_Reoccurrace(newDeg, L_coeffs_new, s_Coeffs_new, false);

                            Polynomial P_partC = Polynomial(s_Coeffs_new, newDeg);

                            std::cout << '\n';
                            P_partC.printPolynomial();
                            std::cout << '\n';

                            vector<ExprPtr> P_partC_roots;
                            bool validPolynomial_C = rootExtractorEMST(P_partC, P_partC_roots);

                            if (validPolynomial_C)
                            {
                                ExprPtr rootSumPtr = L_coeffs[1];
                                Expr::Kind rootSumKind = rootSumPtr->getKind();

                                // Equality Builder?
                                if (rootSumKind == Expr::Kind::Integer || rootSumKind == Expr::Kind::Rational)
                                {

                                    Rational rootSum = getRationalFromNode(rootSumPtr);

                                    ExprPtr newPolyRootSum = Simplifier::simplify(ExprFactory::add(P_partC_roots));

                                    bool rootsFound = false;
                                    ExprPtr rootExclude = NULL;

                                    for (const ExprPtr rootNode : P_partC_roots)
                                    {
                                        // rootNode->printPreorder();

                                        ExprPtr rootMinusNode = Simplifier::simplify(ExprFactory::multiply({INT_MINUS_ONE, rootNode}));
                                        ExprPtr testSum = Simplifier::simplify(ExprFactory::add({newPolyRootSum, rootMinusNode}));
                                        Expr::Kind testSumKind = testSum->getKind();

                                        if (testSumKind == Expr::Kind::Integer || testSumKind == Expr::Kind::Rational)
                                        {
                                            if (rootSum == getRationalFromNode(testSum))
                                            {
                                                rootsFound = true;
                                                rootExclude = rootNode;
                                                break;
                                            }
                                        }
                                    }

                                    if (rootsFound)
                                    {
                                        auto it = std::find(P_partC_roots.begin(), P_partC_roots.end(), rootExclude);

                                        if (it != P_partC_roots.end())
                                        {
                                            P_partC_roots.erase(it);

                                            // generate actual roots now
                                            for (const ExprPtr rootNode : P_partC_roots)
                                            {
                                                ExprPtr actualRoot = Simplifier::simplify({ExprFactory::multiply({ExprFactory::add({rootSumPtr, rootNode}),
                                                                                                                  ExprFactory::rational(1, newDeg)})});
                                                rootsAry.push_back(actualRoot);
                                            }
                                        }
                                        else
                                        {
                                            throw std::runtime_error("Cannot Detect rootExclude");
                                        }
                                    }
                                    else
                                    {
                                        rootSumPtr->printPreorder();
                                        throw std::runtime_error("Cannot Find Valid rootSum");
                                    }
                                }
                                else
                                {
                                    std::cout << "=> FAILED ROOT_SUM_KIND NOT HANDLED\n";
                                    satisfied = false;
                                }
                            }
                            else
                            {
                                std::cout << "=> FAILED RAISING DEGREE FROM " << polyDegree << " TO " << newDeg << " \n";
                                satisfied = false;
                            }
                        }
                    }
                }
                else
                {
                    // first non-zero coefficient occurs between the constant constant and ZERO-(N-1 Deg Term)

                    // GCD Degree filter
                    size_t gcdAll = polyDegree - i;
                    size_t totalNonTerms = 1;
                    vector<ExprPtr> NonZeroPtrAry;

                    NonZeroPtrAry.push_back(P.coeffsAry[0]);
                    NonZeroPtrAry.push_back(P.coeffsAry[i]);

                    for (size_t k = i + 1; k < polyDegree; k++)
                    {
                        ExprPtr termPtr = P.coeffsAry[k];
                        if (!Simplifier::isZeroNode(termPtr))
                        {
                            gcdAll = std::gcd(gcdAll, polyDegree - k);
                            totalNonTerms++;
                            NonZeroPtrAry.push_back(termPtr);
                        }
                    }

                    NonZeroPtrAry.push_back(P.coeffsAry[polyDegree]);

                    if (gcdAll == 0)
                        throw std::runtime_error("GCD ZERO OCCURRED");

                    if (gcdAll > 1)
                    {
                        // create new polynomial with that and solve

                        size_t subDeg = polyDegree / gcdAll; // Degree pf new polynomial
                        std::cout << "=> POLYNOMIAL SUBSTITUTION z = x^" << gcdAll << " \n";

                        vector<ExprPtr> P_partB_coeffs = NonZeroPtrAry;

                        // std::cout << P_partB_coeffs.size() << '\n';

                        Polynomial P_partB = Polynomial(P_partB_coeffs, subDeg);
                        vector<ExprPtr> P_partB_roots;
                        rootExtractorEMST(P_partB, P_partB_roots);

                        for (const auto rootPtr : P_partB_roots)
                        {
                            ExprPtr rootTmp = Simplifier::simplify(ExprFactory::power(rootPtr, ExprFactory::rational(1, subDeg)));
                            rootsAry.push_back(rootTmp);

                            for (size_t k = 1; k < subDeg; k++)
                            {
                                Rational complexCoeff(2 * k, subDeg);
                                ExprPtr complexFactor = Simplifier::unitComplexNumberGen(complexCoeff);
                                rootsAry.push_back(Simplifier::simplify(ExprFactory::multiply({rootTmp, complexFactor})));
                            }
                        }
                    }
                    else
                    {
                        satisfied = false;
                    }
                }
            }
        }

        return satisfied;
    }

    static void solveChebyshev(const Polynomial &P)
    {

        long long polyDeg = P.degree;
        ExprPtr degNode = ExprFactory::integer(polyDeg);
        ExprPtr INV_degNode = Simplifier::simplify(ExprFactory::power(degNode, INT_MINUS_ONE));

        ExprPtr aN = P.coeffsAry[0];
        ExprPtr aN_1 = P.coeffsAry[1];
        ExprPtr aN_2 = P.coeffsAry[2];
        ExprPtr a0 = P.coeffsAry[polyDeg];

        if (!Simplifier::isZeroNode(aN_1))
            throw std::runtime_error("solveChebyshev NOT-SIMPLIFIED FORM");

        ExprPtr amplitude = Simplifier::simplify(
            ExprFactory::multiply({INT_TWO,
                                   ExprFactory::power(
                                       ExprFactory::multiply(
                                           {INT_MINUS_ONE,
                                            aN_2, INV_degNode,
                                            ExprFactory::power(aN, INT_MINUS_ONE)}),
                                       FRACT_ONE_HALF)}));

        std::cout << "\nAmplitude :- \n";

        amplitude->printPreorder();

        ExprPtr scalingFactor = Simplifier::simplify(
            ExprFactory::multiply(
                {aN,
                 ExprFactory::power(amplitude, degNode),
                 ExprFactory::power(INT_TWO, ExprFactory::integer(1 - polyDeg))}));

        ExprPtr INV_scalingFactor = ExprFactory::power(scalingFactor, INT_MINUS_ONE);
        std::cout << "\nScale :- \n";

        scalingFactor->printPreorder();

        ExprPtr angleTerm = Simplifier::simplify(ExprFactory::multiply({INT_MINUS_ONE, a0, INV_scalingFactor}));
        if (polyDeg % 2 == 0)
        {
            ExprPtr negFactor = Simplifier::simplify(ExprFactory::power(INT_MINUS_ONE, ExprFactory::integer(polyDeg / 2)));
            angleTerm = Simplifier::simplify(ExprFactory::add({negFactor, angleTerm}));
        }

        std::cout << "\nFor arcos(theta), THETA :- \n";
        angleTerm->printPreorder();
    }

    bool static rootExtractor(Polynomial &P, vector<ExprPtr> &rootsAry, bool simpleCheck, bool debug)
    {
        bool satisfied = true;

        size_t polyDegree = P.degree; // size = degree

        // CASE 1 : x^2 + x = 0 absence of last coefficient
        ExprPtr lastCoeffPtr = P.coeffsAry[polyDegree];
        Expr::Kind lastCoeffExprKind = lastCoeffPtr->getKind();

        bool absenceOfLastCoeff = false;

        if (lastCoeffExprKind == Expr::Kind::Integer)
        {
            auto node = std::static_pointer_cast<IntegerNode>(lastCoeffPtr);
            if (node->value == 0)
            {
                absenceOfLastCoeff = true;
            }
        }

        size_t i = 1;

        bool allCoeffZero = true;

        while (i <= polyDegree && allCoeffZero)
        {
            if (Simplifier::isZeroNode(P.coeffsAry[i]))
                i++;
            else
                allCoeffZero = false;
        }

        // i points to first non-zero coefficient occurance

        if (allCoeffZero)
        {
            // CASE 1 : x^2 = 0

            if (debug)
            {
                std::cout << "\n=> SOLVING POLYNOMIAL x^" << polyDegree << " = 0\n\n";

                P.printPolynomial();
                std::cout << '\n';
            }

            // create vector of all 0s
            rootsAry.resize(polyDegree, INT_ZERO);
        }
        else
        {
            if (absenceOfLastCoeff)
            {
                // CASE 2 : x^6 + x^2 = 0 condition

                size_t k1 = i + 1; // Get the deg of last non-zero coeff TERM
                for (; k1 <= polyDegree; k1++)
                {
                    if (!Simplifier::isZeroNode(P.coeffsAry[k1]))
                        i = k1;
                }

                size_t partA_deg = polyDegree - i;

                if (debug)
                {
                    std::cout << "\n=> FACTORING OUT x^" << partA_deg << " TERM\n";
                }

                size_t partB_deg = polyDegree - partA_deg;

                for (size_t k = 1; k <= partA_deg; k++)
                {
                    rootsAry.push_back(INT_ZERO);
                }

                // iterate original coeffAry and generate a new polynomial with all 0s, replace those replaces which have non-zero coeff
                vector<ExprPtr> partB_coeffAry;
                partB_coeffAry.resize(partB_deg + 1, INT_ZERO);

                for (size_t k = 0; k < polyDegree; k++)
                {
                    if (!Simplifier::isZeroNode(P.coeffsAry[k]))
                        partB_coeffAry[k] = P.coeffsAry[k];
                }

                Polynomial P_partB = Polynomial(partB_coeffAry, partB_deg);
                // P_partB.printPolynomial();

                vector<ExprPtr> P_partB_roots;
                satisfied = rootExtractor(P_partB, P_partB_roots, simpleCheck, debug);

                rootsAry.insert(rootsAry.end(), P_partB_roots.begin(), P_partB_roots.end());
            }
            else
            {
                ExprPtr INV_POLYDEG = ExprFactory::rational(1, polyDegree);
                if (i == polyDegree) // first noon-zero term is the constant itself
                {

                    if (debug)
                    {
                        std::cout << "\n=> SOLVING POLYNOMIAL x^" << polyDegree << " + k = 0\n\n";

                        P.printPolynomial();
                        std::cout << '\n';
                    }

                    // CASE 3 : x^2 - 1 = 0 condition
                    ExprPtr numFactor = Simplifier::simplify(ExprFactory::power(
                        ExprFactory::multiply({INT_MINUS_ONE, P.coeffsAry[i], ExprFactory::power(P.coeffsAry[0], INT_MINUS_ONE)}),
                        INV_POLYDEG));

                    // Generate nth root of unitites for eaech multiply numFactor and add to Roots List

                    rootsAry.push_back(numFactor); // always added

                    for (size_t k = 1; k < polyDegree; k++)
                    {
                        Rational complexCoeff(2 * k, polyDegree);
                        ExprPtr complexFactor = Simplifier::unitComplexNumberGen(complexCoeff);

                        rootsAry.push_back(Simplifier::simplify(ExprFactory::multiply({numFactor, complexFactor})));
                    }
                }
                else if (i == 1) // second highest degree is present
                {

                    if (polyDegree == 2)
                    {
                        if (debug)
                        {
                            std::cout << "\n=> SOLVING QUADRATIC POLYNOMIAL DEG=2 \n\n";

                            P.printPolynomial();
                            std::cout << '\n';
                        }

                        ExprPtr a = P.coeffsAry[0];
                        ExprPtr b = P.coeffsAry[1];
                        ExprPtr c = P.coeffsAry[2];

                        ExprPtr TWICE_A_INV = Simplifier::simplify(ExprFactory::power(
                            ExprFactory::multiply({INT_TWO, a}),
                            INT_MINUS_ONE));

                        ExprPtr discriminant = Simplifier::simplify(ExprFactory::power(
                            ExprFactory::add({ExprFactory::power(b, INT_TWO),
                                              ExprFactory::multiply({ExprFactory::integer(-4),
                                                                     a,
                                                                     c})}),
                            FRACT_ONE_HALF));

                        ExprPtr minus_B = Simplifier::simplify(ExprFactory::multiply({INT_MINUS_ONE, b}));

                        rootsAry.push_back(Simplifier::simplify({ExprFactory::multiply({ExprFactory::add({minus_B, discriminant}),
                                                                                        TWICE_A_INV})}));

                        rootsAry.push_back(Simplifier::simplify({ExprFactory::multiply({ExprFactory::add({minus_B, ExprFactory::multiply({INT_MINUS_ONE, discriminant})}),
                                                                                        TWICE_A_INV})}));
                    }
                    else
                    {

                        if (debug)
                        {
                            std::cout << "\n=> CENTERING THE POLYNOMIAL DEG=" << polyDegree << " \n\n";

                            P.printPolynomial();
                            std::cout << '\n';
                        }

                        ExprPtr normalizationFactor;

                        Polynomial P_partB = Polynomial(P, normalizationFactor);

                        if (debug)
                        {
                            std::cout << "\n=> Normalization Factor" << polyDegree << " \n\n";
                            normalizationFactor->printPreorder();

                            std::cout << "\n=> CENTERED POLYNOMIAL \n\n";

                            P_partB.printPolynomial();
                            std::cout << '\n';
                        }

                        vector<ExprPtr> P_partB_roots;
                        satisfied = rootExtractor(P_partB, P_partB_roots, simpleCheck, debug);

                        for (size_t k = 0; k < polyDegree; k++)
                        {
                            P_partB_roots[k] = Simplifier::simplify(ExprFactory::add({P_partB_roots[k], normalizationFactor}));
                        }

                        rootsAry.insert(rootsAry.end(), P_partB_roots.begin(), P_partB_roots.end());

                        /*

                        if (polyDegree % 2 == 0)
                        {

                         }
                        else
                        {

                            size_t newDeg = polyDegree + 1;

                            if (debug)
                            {
                                std::cout << "=> RAISING DEGREE FROM " << polyDegree << " TO " << newDeg << " \n";
                            }

                            Polynomial::makeMonic(P);
                            vector<ExprPtr> newCoeffsPoly = P.coeffsAry;
                            newCoeffsPoly.push_back(INT_ZERO);

                            Polynomial newPoly(newCoeffsPoly, newDeg);

                            if (debug)
                            {

                                std::cout << '\n';
                                newPoly.printPolynomial();
                                std::cout << '\n';
                            }

                            vector<ExprPtr> L_coeffs_new;
                            L_coeffs_gen(newPoly, L_coeffs_new);

                            vector<ExprPtr> s_Coeffs_new;
                            EMST_Reoccurrace(newDeg, L_coeffs_new, s_Coeffs_new, false);

                            Polynomial P_partC = Polynomial(s_Coeffs_new, newDeg);

                            if (debug)
                            {
                                std::cout << "=> EMST RECURRANCE DEG=" << polyDegree << " \n";
                                std::cout << '\n';
                                P_partC.printPolynomial();
                                std::cout << '\n';
                            }

                            vector<ExprPtr> P_partC_roots;
                            satisfied = rootExtractor(P_partC, P_partC_roots, simpleCheck, debug);
                        }
                    */

                        // !------------ APPROACH 2 EZOUIDI NORMALIZATION -------!
                        /*
                                                if (debug)
                                                {
                                                    std::cout << "\n=> APPLYING EMST THEOREM TO DEG=" << polyDegree << " \n\n";

                                                    P.printPolynomial();
                                                    std::cout << '\n';
                                                }

                                                // Do ezoudi Normalization
                                                ExprPtr ezouidiL0;
                                                Polynomial P_EC = Polynomial(P, ezouidiL0, true);
                                                // P_EC.printPolynomial();

                                                if (debug)
                                                {
                                                    std::cout << "\n=> SIMPLIFIED EMST POLYNOMIAL \n\n";

                                                    P_EC.printPolynomial();
                                                    std::cout << '\n';
                                                }

                                                vector<ExprPtr> P_partB_roots;
                                                satisfied = rootExtractor(P_EC, P_partB_roots, simpleCheck, debug);
                                                */
                    }
                }
                else
                {
                    // first non-zero coefficient occurs between the constant constant and ZERO-(N-1 Deg Term)

                    // GCD Degree filter
                    size_t gcdAll = polyDegree - i;
                    size_t totalNonTerms = 1;
                    vector<ExprPtr> NonZeroPtrAry;

                    NonZeroPtrAry.push_back(P.coeffsAry[0]);
                    NonZeroPtrAry.push_back(P.coeffsAry[i]);

                    for (size_t k = i + 1; k < polyDegree; k++)
                    {
                        ExprPtr termPtr = P.coeffsAry[k];
                        if (!Simplifier::isZeroNode(termPtr))
                        {
                            gcdAll = std::gcd(gcdAll, polyDegree - k);
                            totalNonTerms++;
                            NonZeroPtrAry.push_back(termPtr);
                        }
                    }

                    NonZeroPtrAry.push_back(P.coeffsAry[polyDegree]);

                    if (gcdAll == 0)
                        throw std::runtime_error("GCD ZERO OCCURRED");

                    if (gcdAll > 1)
                    {
                        // create new polynomial with that and solve

                        size_t subDeg = polyDegree / gcdAll;
                        if (debug)
                            std::cout << "\n=> POLYNOMIAL SUBSTITUTION z = x^" << gcdAll << " \n\n";

                        vector<ExprPtr> P_partB_coeffs = NonZeroPtrAry;
                        Polynomial P_partB = Polynomial(P_partB_coeffs, subDeg);

                        if (debug)
                        {
                            P_partB.printPolynomial();
                            std::cout << '\n';
                        }

                        vector<ExprPtr> P_partB_roots;
                        satisfied = rootExtractor(P_partB, P_partB_roots, simpleCheck, debug);

                        for (const auto rootPtr : P_partB_roots)
                        {
                            ExprPtr rootTmp = Simplifier::simplify(ExprFactory::power(rootPtr, ExprFactory::rational(1, subDeg)));
                            rootsAry.push_back(rootTmp);

                            for (size_t k = 1; k < subDeg; k++)
                            {
                                Rational complexCoeff(2 * k, subDeg);
                                ExprPtr complexFactor = Simplifier::unitComplexNumberGen(complexCoeff);
                                rootsAry.push_back(Simplifier::simplify(ExprFactory::multiply({rootTmp, complexFactor})));
                            }
                        }
                    }
                    else
                    {

                        /*
                        if (simpleCheck)
                            return false;
                        // pure non-substitutional polynomial do methods

                        // check possiblities for trinomials

                        // check for quadratic form (even degree polynomials only)

                        // if not apply decoupled decomposition, also print for what degree you enconuntered this thing

                        if (debug)
                        {
                            std::cout << "\n=> IRREDUCIBLE POLYNOMIAL \n\n";

                            P.printPolynomial();
                            std::cout << "\n=> MAKING MONIC \n\n";
                        }

                        // ensure divding every term by aN term
                        Polynomial::makeMonic(P);

                        if (debug)
                        {
                            P.printPolynomial();

                            std::cout << "\n=> RUNNING POWER SUM CHECKER \n\n";
                        }

                        bool foundGoodOne = false;
                        size_t Q = 2;

                        while (!foundGoodOne && Q < 50)
                        {
                            std::cout << "Q=" << Q << " \n";
                            Polynomial P_partB = Polynomial(P, Q);

                            vector<ExprPtr> P_partB_roots;
                            foundGoodOne = rootExtractor(P_partB, P_partB_roots, true, false);

                            if (foundGoodOne)
                            {
                                std::cout << "\nSOLVABLE POWER POLYNOMIAL Q=" << Q << " \n";
                                P_partB.printPolynomial();

                                for (const auto m : P_partB_roots)
                                {
                                    std::cout << "\nRoot \n";
                                    m->printPreorder();
                                }
                            }

                            Q++;
                        }
                        */

                        // find weighted invarients check cases

                        /*
                        vector<ExprPtr> powerSumAry;
                        PowerSumGenerator(P, powerSumAry);

                        std::cout << "\n--- Power Sums DEG=" << polyDegree << " ----\n";
                        for (ExprPtr nodePtr : powerSumAry)
                        {
                            nodePtr->printPreorder();
                        }

                        // directly solve for cubic/quadratic systems

                        if (polyDegree == 3)
                        {
                        }
                        else if (polyDegree == 4)
                        {
                        }

                        // Build Hankel system and confirm for distinct roots

                        std::cout << "\n--- Power Sums DEG=" << polyDegree << " ----\n";

                        size_t d = polyDegree;
                        vector<vector<ExprPtr>> HankelMatrixFinal;

                        while (d >= 2 && HankelMatrixFinal.empty())
                        {
                            vector<vector<ExprPtr>> HankelMatrix(d, vector<ExprPtr>(d));

                            for (size_t r = 0; r < d; r++)
                            {
                                for (size_t c = 0; c < d; c++)
                                {
                                    HankelMatrix[r][c] = powerSumAry[r + c];
                                }
                            }

                            std::cout << "\nDetermiant of H(" << d << ")\n";
                            ExprPtr detValue = determinantBareiss(HankelMatrix, d);
                            detValue->printPreorder();
                            if (!Simplifier::isZeroNode(detValue))
                                HankelMatrixFinal = HankelMatrix;
                            else
                                d--;
                        }

                        std::cout << "\n--- Hankel Matrix Calculated ----\n";

                        for (size_t r = 0; r < d; r++)
                        {
                            for (size_t c = 0; c < d; c++)
                            {
                                HankelMatrixFinal[r][c]->printPreorder();
                            }
                            std::cout << "\n---\n";
                        }

                        if (HankelMatrixFinal.empty())
                            throw std::runtime_error("NO DET WAS NON-ZERO IN HANKEL SCAN");

                        if (d == polyDegree)
                        {
                            std::cout << "=> ALL ROOTS ARE DISTINCT - HANKEL SCAN \n";
                        }
                        else
                        {
                            std::cout << "=> REPEATED ROOTS - HANKEL SCAN \n";
                        }
                        */

                        /* bool isRotationalSym = true;
                         bool isAllZero = false;
                         size_t k = 0; // first non-zero coefficient

                         while (isRotationalSym && k < polyDegree - 1)
                         {
                             if (!Simplifier::isZeroNode(weightInvAry[k]))
                                 isRotationalSym = false;
                             else
                                 k++;
                         }

                         if (Simplifier::isZeroNode(weightInvAry[polyDegree - 1]))
                         {
                             isRotationalSym = false;
                             isAllZero = true;
                         }

                         if (isRotationalSym)
                         {
                             std::cout << "=> ROTATIONAL SYMMETRY IN ROOTS DEG=" << polyDegree << " \n";

                             vector<ExprPtr> partB_coeffs(polyDegree + 1, INT_ZERO);
                             partB_coeffs[0] = INT_ONE;
                             partB_coeffs[polyDegree] = Simplifier::simplify(ExprFactory::multiply({INT_MINUS_ONE, weightInvAry[polyDegree - 1]}));

                             Polynomial P_partB(partB_coeffs, polyDegree);
                             vector<ExprPtr> partB_roots;
                             rootExtractor(P_partB, partB_roots);

                             rootsAry.insert(rootsAry.end(), partB_roots.begin(), partB_roots.end());
                         }
                         else
                         {
                             vector<ExprPtr> nonZeroWeights;
                             nonZeroWeights.push_back(weightInvAry[k]);

                             // check partial symmetry, measure distances between nodes, if consistent than true else chaotic

                             bool isPartialSym = true;

                             size_t gap = 0;
                             size_t prevIndex = k;

                             k = k + 1;
                             while (Simplifier::isZeroNode(weightInvAry[k]))
                                 k++;
                             gap = k - prevIndex;
                             prevIndex = k;
                             nonZeroWeights.push_back(weightInvAry[k]);

                             for (k = k + 1; k < polyDegree && isPartialSym; k++)
                             {
                                 if (!Simplifier::isZeroNode(weightInvAry[k]))
                                 {
                                     size_t newGap = k - prevIndex;
                                     if (newGap != gap)
                                         isPartialSym = false;
                                     else
                                         prevIndex = k;

                                     nonZeroWeights.push_back(weightInvAry[k]);
                                 }
                             }

                             if (gap == 1)
                             {
                                 std::cout << "=> CHAOTIC/REPEATED ROOTS DEG=" << polyDegree << " \n";

                                 // do check Hankel matrix Repeated roots else its chaotic system
                             }
                             else if (gap > 1 && isPartialSym)
                             {
                                 // partial symmetry observed
                                 std::cout << "=> PARTIAL SYMMETRY IN ROOTS DEG=" << polyDegree << " \n";

                                 // check if only two nonzeroInvarients?
                             }
                             else
                             {
                                 std::cout << "=> CHAOTIC ROOTS DEG=" << polyDegree << " \n";

                                 // chaotic system
                             }
                         }
    */

                        if (debug)
                        {
                            std::cout << "\n=> CHEBYSHEV POLYNOMIAL SOLUTION \n\n";
                        }

                        solveChebyshev(P);

                        satisfied = true;

                        rootsAry.resize(polyDegree, INT_ZERO);

                        // check for TRINOMIAL CASES and quadratic, cubic, quintic in any x^k ;

                        // we express z = u + v;

                        // vector<ExprPtr> strcturalNodes;

                        // std::vector<ExprPtr> skNodes(P.coeffsAry.begin() + 1, P.coeffsAry.end());
                        // std::reverse(skNodes.begin(), skNodes.end());

                        // StructrualRemainder(polyDegree, skNodes, strcturalNodes);

                        // Quadratic of phi
                        // ExprPtr A =
                    }
                }
            }
        }

        return satisfied;
    }
};

// ============================================================
// MAIN PROGRAM
// ============================================================

// Simplifier
int simplifierTest()
{

    ExprPtr expr0 = ExprFactory::power(ExprFactory::integer(144), ExprFactory::rational(1, 2));

    ExprPtr expr1 = ExprFactory::add({ExprFactory::integer(20), ExprFactory::integer(-2), ExprFactory::rational(1, 2), ExprFactory::multiply({ExprFactory::integer(2), IOTA}), IOTA});

    ExprPtr expr2 = ExprFactory::multiply(
        {ExprFactory::add({INT_ONE, ExprFactory::multiply({ExprFactory::integer(2), IOTA})}),
         ExprFactory::add({INT_ONE, ExprFactory::multiply({ExprFactory::integer(3), IOTA})})});

    ExprPtr expr3 = ExprFactory::multiply(
        {ExprFactory::add({INT_ONE, ExprFactory::power(ExprFactory::integer(2), ExprFactory::rational(1, 2))}),
         ExprFactory::add({INT_ONE, ExprFactory::power(ExprFactory::integer(2), ExprFactory::rational(1, 2))})});

    ExprPtr expr4 = ExprFactory::power(INT_ZERO, INT_ZERO);

    ExprPtr expr5 = ExprFactory::power(ExprFactory::integer(12), ExprFactory::integer(23)); // negative answer

    ExprPtr expr6 = ExprFactory::power(ExprFactory::integer(12), ExprFactory::integer(2));

    ExprPtr expr7 = ExprFactory::power(ExprFactory::rational(12, 24), ExprFactory::integer(-2));

    ExprPtr expr8 = ExprFactory::power(ExprFactory::rational(12, 24), ExprFactory::rational(2, 4));

    ExprPtr expr9 = ExprFactory::power(IOTA, ExprFactory::integer(4));

    ExprPtr expr10 = ExprFactory::multiply({ExprFactory::rational(12, 24), IOTA, IOTA});

    ExprPtr expr11 = ExprFactory::power(ExprFactory::integer(144), ExprFactory::rational(5, 4));

    ExprPtr expr12 = ExprFactory::power(ExprFactory::integer(-2), ExprFactory::rational(-1, 2));

    ExprPtr expr13 = ExprFactory::power(ExprFactory::add({INT_ONE, ExprFactory::power(ExprFactory::integer(2), ExprFactory::rational(1, 2))}),

                                        ExprFactory::rational(5, 2));

    ExprPtr expr14 = ExprFactory::power(ExprFactory::add(
                                            {INT_ONE,
                                             ExprFactory::power(ExprFactory::integer(2), ExprFactory::rational(1, 2)),
                                             ExprFactory::power(ExprFactory::integer(3), ExprFactory::rational(1, 2))}),
                                        ExprFactory::integer(3));

    // wrong - (complex)^(complex)
    ExprPtr expr15 = ExprFactory::power(ExprFactory::add({INT_ONE, ExprFactory::multiply({ExprFactory::integer(2), IOTA})}),

                                        ExprFactory::add({INT_ONE, ExprFactory::multiply({ExprFactory::integer(2), IOTA})}));

    ExprPtr expr00 = ExprFactory::rational(13440, 256);

    ExprPtr expr16 = (ExprFactory::multiply({INT_MINUS_ONE, INT_ONE}));

    ExprPtr expr17 = (ExprFactory::power(
        ExprFactory::multiply({INT_MINUS_ONE, INT_ONE, ExprFactory::power(INT_TWO, INT_MINUS_ONE)}),
        ExprFactory::rational(1, 2)));

    ExprPtr expr18 = (ExprFactory::multiply({INT_MINUS_ONE, INT_ONE, ExprFactory::power(INT_TWO, INT_MINUS_ONE)}));

    ExprPtr expr19 = ExprFactory::power(INT_TWO, INT_MINUS_ONE);

    ExprPtr expr20 = ExprFactory::multiply({ExprFactory::power(INT_THREE, FRACT_ONE_HALF, true),
                                            ExprFactory::power(INT_TWO, INT_MINUS_ONE)});

    ExprPtr ROOT2 = ExprFactory::power(INT_TWO, FRACT_ONE_HALF);
    ExprPtr addTerm1 = ExprFactory::add({INT_ONE, ROOT2});
    ExprPtr MINUS_ONE_FOURTH = ExprFactory::rational(-1, 4);
    ExprPtr addTermRat = ExprFactory::power(addTerm1, MINUS_ONE_FOURTH);

    ExprPtr expr21 = ExprFactory::multiply({addTerm1, addTermRat, addTermRat, addTermRat, addTermRat, addTermRat});

    int R = 6;

    vector<ExprPtr> skNodes;
    skNodes.resize(R);

    for (int k = 0; k < R; k++)
    {
        skNodes[k] = ExprFactory::symbol("s" + k);
    }

    ExprPtr expr22 = EMST::HFunction(R, skNodes);

    ExprPtr expr23 = ExprFactory::multiply({IOTA, ExprFactory::multiply({IOTA, FRACT_ROOT3_BYTWO})});
    ExprPtr expr24 = ExprFactory::multiply({IOTA, ExprFactory::add({FRACT_ONE_HALF, ExprFactory::multiply({IOTA, FRACT_ROOT3_BYTWO})})});

    ExprPtr x = ExprFactory::symbol("x");
    ExprPtr y = ExprFactory::symbol("y");
    ExprPtr expr25 = ExprFactory::multiply({IOTA, x, y, IOTA, x, y, IOTA});

    ExprPtr expr26 = ExprFactory::multiply({ExprFactory::add({FRACT_ONE_HALF, ExprFactory::multiply({IOTA, FRACT_ROOT3_BYTWO})}), IOTA, ExprFactory::multiply({IOTA, IOTA, IOTA})});

    ExprPtr expr27 = ExprFactory::multiply({ExprFactory::add({INT_ONE, IOTA}), ExprFactory::add({FRACT_ONE_HALF, IOTA})});

    ExprPtr expr28 = ExprFactory::power(INT_MINUS_ONE, FRACT_ONE_HALF);

    ExprPtr ROOT6 = ExprFactory::power(ExprFactory::integer(6), FRACT_ONE_HALF);

    ExprPtr expr = ExprFactory::multiply({IOTA, IOTA, ROOT6, ROOT6, ExprFactory::rational(405, 2), ExprFactory::rational(135, 2)});

    // EXAMPLE 1 - works

    ExprPtr complexPart = ExprFactory::power(
        ExprFactory::add({ExprFactory::integer(3),
                          ExprFactory::multiply({ExprFactory::integer(2), IOTA})}),
        ExprFactory::integer(2));

    ExprPtr xPlusX = ExprFactory::add({x, x});
    ExprPtr sixxy = ExprFactory::multiply({ExprFactory::integer(6), x, y});
    ExprPtr twoxy = ExprFactory::multiply({ExprFactory::integer(2), x, y});
    ExprPtr sqrt9 = ExprFactory::power(ExprFactory::integer(9), ExprFactory::rational(1, 2));

    ExprPtr exprEx1 = ExprFactory::add({complexPart, xPlusX, sixxy, twoxy, sqrt9});

    // Example 2 - Not simplified form

    ExprPtr symbolic = ExprFactory::add({ExprFactory::multiply({ExprFactory::integer(6), x, y}),
                                         ExprFactory::multiply({ExprFactory::integer(2), x, y}),
                                         ExprFactory::add({x, x}),
                                         sqrt9});

    ExprPtr euler = ExprFactory::power(ExprFactory::constant("e"), ExprFactory::multiply({IOTA, PIE}));

    ExprPtr exprEx2 = ExprFactory::add({complexPart,
                                        symbolic,
                                        euler});

    std::cout << "Original tree (preorder):\n";
    expr->printPreorder();

    ExprPtr simp = Simplifier::simplify(expr);

    std::cout << "\nSimplified tree (preorder):\n";
    simp->printPreorder();

    return 0;
}

// Root Validator
void RootValidator()
{
    size_t D = 20;

    // they are exact coefficients of polynomial
    vector<ExprPtr> coeffs = {INT_ONE,
                              ExprFactory::integer(-14),
                              ExprFactory::integer(84),
                              ExprFactory::integer(-623),
                              ExprFactory::integer(3304),
                              ExprFactory::integer(-8904),
                              ExprFactory::integer(11424),
                              ExprFactory::integer(-5616)};

    ExprPtr root = Simplifier::simplify(ExprFactory::add({ExprFactory::integer(2), ExprFactory::power(ExprFactory::integer(343), ExprFactory::rational(1, 3))}));

    vector<ExprPtr> mulNodes;

    ExprPtr rootPrevPowerTerm = INT_ONE;

    for (int i = D; i >= 0; i--)
    {
        ExprPtr currRootPowerTerm = Simplifier::simplify(ExprFactory::multiply({rootPrevPowerTerm, root}));

        mulNodes.push_back(Simplifier::simplify(ExprFactory::multiply({currRootPowerTerm, coeffs[i]})));

        rootPrevPowerTerm = currRootPowerTerm;
    }

    ExprPtr result = Simplifier::simplify(ExprFactory::add(mulNodes));

    std::cout << "\nRoot Substitution Result :\n";
    result->printPreorder();
}

Polynomial generateHermite(size_t degree)
{
    std::vector<cpp_int> H0 = {1};
    std::vector<cpp_int> H1 = {2};

    if (degree == 0)
    {
        std::vector<ExprPtr> coeffs;
        coeffs.push_back(ExprFactory::integer(1));

        return Polynomial(coeffs, 0);
    }

    if (degree == 1)
    {
        std::vector<ExprPtr> coeffs;
        coeffs.push_back(ExprFactory::integer(2));
        coeffs.push_back(ExprFactory::integer(0));

        return Polynomial(coeffs, 1);
    }

    std::vector<cpp_int> Hnm2 = H0;
    std::vector<cpp_int> Hnm1 = H1;
    std::vector<cpp_int> Hn;

    for (size_t n = 2; n <= degree; n++)
    {
        Hn.assign(n + 1, 0);

        // 2x * H(n-1)
        for (size_t i = 0; i < Hnm1.size(); i++)
        {
            Hn[i] += 2 * Hnm1[i];
        }

        // -2(n-1) * H(n-2)
        cpp_int factor = -2 * cpp_int(n - 1);

        for (size_t i = 0; i < Hnm2.size(); i++)
        {
            Hn[i + 2] += factor * Hnm2[i];
        }

        Hnm2 = Hnm1;
        Hnm1 = Hn;
    }

    std::vector<ExprPtr> coeffs;

    for (const cpp_int &c : Hnm1)
    {
        coeffs.push_back(ExprFactory::integer(c));
    }

    return Polynomial(coeffs, degree);
}

// PolyTestes

void PolyTester()
{

    // size_t degree = 11;

    // Polynomial HermitePoly = generateHermite(degree);

    // HermitePoly.printPolynomial();

    // Polynomial P5 = Polynomial({intg(1024),
    //                             intg(-23040),
    //                             intg(161280),
    //                             intg(-403200),
    //                             intg(302400),
    //                             intg(-30240)},
    //                            5);

    ExprPtr root6 = powg(intg(6), FRACT_ONE_HALF, true);

    Polynomial P5 = Polynomial({INT_ONE,
                                INT_ZERO,
                                INT_ZERO,
                                intg(21780),
                                addg({ratg(2190375, 8), mulg({ratg(-289575, 2), root6, IOTA}, true)}, true),
                                addg({ratg(-9497169, 4), mulg({ratg(-7798923, 4), root6, IOTA}, true)}, true)},
                               5);

    Polynomial P5_2 = Polynomial({INT_ONE,
                                  intg(55), intg(-1150), intg(-65500), intg(405000), intg(9000000)

                                 },
                                 5);

    // P5.printPolynomial();

    // Polynomial::makeMonic(P5);

    // P5.printPolynomial();

    std::cout
        << "\n ----LOGS --- :\n\n";

    try
    {

        // ExprPtr ezouidiL0;
        // Polynomial P12_ezoudi = Polynomial(P5, ezouidiL0, true);
        // P12_ezoudi.printPolynomial();

        vector<ExprPtr> rootsAry;

        EMST::rootExtractor(P5_2, rootsAry, false, true);
        /*
                vector<ExprPtr> rootsAry;

                EMST::rootExtractor(HermitePoly, rootsAry, false, true);
        */
        // P7_ez.printPolynomial();

        /*
        std::cout << "\n ----Roots of Polynomial --- :\n\n";

        for (size_t k = 0; k < 2; k++)
        {
            std::cout << "Root " << k + 1 << " : \n";
            rootsAry[k]->printPreorder();
            std::cout << "\n";
        }
*/

        //
        // Polynomial P_Q2 = Polynomial(P4, 8);
        //   P_Q2.printPolynomial();
        //
        //   Polynomial P_Q2_EC = Polynomial(P_Q2, ezouidiL0, true);
        //   P_Q2_EC.printPolynomial();

        /*  // TEST 2
          ExprPtr ezouidiL0;
          Polynomial P_EC = Polynomial(P3, ezouidiL0, true);
          P_EC.printPolynomial();
        */

        /* vector<ExprPtr> rootsAry;

             EMST::rootExtractor(HermitePoly, rootsAry, false, true);

                   std::cout << "\n ----Roots of Polynomial --- :\n\n";

                  for (size_t k = 0; k < degree; k++)
                   {
                       std::cout << "Root " << k + 1 << " : \n";
                       rootsAry[k]->printPreorder();
                       std::cout << "\n";
                   }
        */
    }
    catch (std::runtime_error e)
    {
        std::cout << e.what();
    }
}

// Polynomial Root Extractor
int RootExtractorDriver()
{

    Polynomial P3_EMST({ExprFactory::integer(1),
                        ExprFactory::integer(-11),
                        ExprFactory::integer(96),
                        ExprFactory::rational(-2893, 8)},
                       3);

    Polynomial P4({ExprFactory::integer(1),
                   ExprFactory::integer(-12),
                   ExprFactory::integer(54),
                   ExprFactory::integer(-108),
                   ExprFactory::integer(81)},
                  4);

    Polynomial P4_Hankel({ExprFactory::integer(1),
                          ExprFactory::integer(-1),
                          ExprFactory::integer(-7),
                          ExprFactory::integer(13),
                          ExprFactory::integer(-6)},
                         4);

    Polynomial P6_H({ExprFactory::integer(64),
                     ExprFactory::integer(0),
                     ExprFactory::integer(-480),
                     ExprFactory::integer(0),
                     ExprFactory::integer(720),
                     ExprFactory::integer(0),
                     ExprFactory::integer(-120)},
                    6);

    Polynomial P7_H({ExprFactory::integer(128),
                     ExprFactory::integer(-1344),
                     ExprFactory::integer(3360),
                     ExprFactory::integer(0),
                     ExprFactory::integer(720),
                     ExprFactory::integer(0),
                     ExprFactory::integer(-120)},
                    7);

    Polynomial P7({INT_ONE,
                   ExprFactory::integer(14),
                   ExprFactory::integer(84),
                   ExprFactory::integer(623),
                   ExprFactory::integer(3304),
                   ExprFactory::integer(8904),
                   ExprFactory::integer(11424),
                   ExprFactory::integer(5616)},
                  7);

    Polynomial P8_EMST({INT_ONE,
                        ExprFactory::integer(-8),
                        ExprFactory::integer(28),
                        ExprFactory::integer(-56),
                        ExprFactory::integer(70),
                        ExprFactory::integer(-56),
                        ExprFactory::integer(28),
                        ExprFactory::integer(-8), INT_ONE},
                       8);

    Polynomial P8_EMST2({INT_ONE,
                         ExprFactory::integer(-16),
                         ExprFactory::rational(7204, 64),
                         ExprFactory::rational(-232832, 512),
                         ExprFactory::rational(4726206, 4096),
                         ExprFactory::rational(-61697920, 32768),
                         ExprFactory::rational(505838740, 262144),
                         ExprFactory::rational(-2381353600, 2097152),
                         ExprFactory::rational(4928741769, 16777216)},
                        8);

    Polynomial P5_EMST({INT_ONE,
                        INT_MINUS_ONE,
                        ExprFactory::rational(15, 36),
                        INT_MINUS_ONE,
                        ExprFactory::rational(603, 1296),
                        ExprFactory::rational(-594, 7776)},
                       5);

    vector<ExprPtr> rootsAry;

    size_t degree = 10;

    Polynomial P3_EMST2({ExprFactory::integer(8),
                         ExprFactory::integer(-60),
                         ExprFactory::integer(-810),
                         ExprFactory::integer(-285)},
                        3);

    std::cout << "\n ----Solving of Polynomial OF DEG=" << degree << " --- :\n\n";

    Polynomial HermitePoly = generateHermite(degree);

    HermitePoly.printPolynomial();

    std::cout << "\n ----LOGS --- :\n\n";

    try
    {
        // EMST::rootExtractorEMST(HermitePoly, rootsAry);
        // EMST::rootExtractor(HermitePoly, rootsAry);
    }
    catch (std::runtime_error e)
    {
        std::cout << e.what();
    }

    // EMST::rootExtractor(P6_H, rootsAry);

    std::cout << "\n ----Roots of Polynomial --- :\n\n";

    for (size_t k = 0; k < degree; k++)
    {
        std::cout << "Root " << k + 1 << " : \n";
        rootsAry[k]->printPreorder();
        std::cout << "\n";
    }
    return 0;
}

int mainPrime(cpp_int n)
{
    auto f = factorize(n);

    for (size_t i = 0; i < f.size(); ++i)
    {
        if (i)
            std::cout << " * ";
        std::cout << f[i].first;
        if (f[i].second > 1)
            std::cout << "^" << f[i].second;
    }
    std::cout << '\n';

    return 0;
}

void quintic_formula()
{

    // Polynomial P4 = Polynomial({INT_ONE, intg(-7), intg(19), intg(-25), intg(16), intg(-4)}, 5);

    // ExprPtr p = intg(-55), q = intg(-110), r = ratg(1155, 4), s = intg(473);
    ExprPtr p = INT_ZERO;
    ExprPtr q = INT_ZERO;
    ExprPtr r = intg(15);
    ExprPtr s = intg(12);

    // 6th degree resolvent
    // Base variable powers
    ExprPtr p2 = S(p ^ 2);
    ExprPtr p3 = S(p ^ 3);
    ExprPtr p4 = S(p ^ 4);
    ExprPtr p5 = S(p ^ 5);
    ExprPtr p6 = S(p ^ 6);
    ExprPtr p7 = S(p ^ 7);

    ExprPtr q2 = S(q ^ 2);
    ExprPtr q3 = S(q ^ 3);
    ExprPtr q4 = S(q ^ 4);
    ExprPtr q5 = S(q ^ 5);
    ExprPtr q6 = S(q ^ 6);
    ExprPtr q8 = S(q ^ 8);

    ExprPtr r2 = S(r ^ 2);
    ExprPtr r3 = S(r ^ 3);
    ExprPtr r4 = S(r ^ 4);
    ExprPtr r5 = S(r ^ 5);
    ExprPtr r6 = S(r ^ 6);

    ExprPtr s2 = S(s ^ 2);
    ExprPtr s3 = S(s ^ 3);
    ExprPtr s4 = S(s ^ 4);

    // Frequently repeating grouped terms
    ExprPtr p_q2 = S(p * q2);
    ExprPtr q_s = S(q * s);
    ExprPtr q_r_s = S(q * r * s);

    // Coefficient variables
    ExprPtr a5_C = 8 * r;

    ExprPtr a4_C = (2 * p_q2) - (6 * p2 * r) + (40 * r2) - (50 * q_s);

    ExprPtr a3_C = (-2 * q4) + (21 * p_q2 * r) - (40 * p2 * r2) + (160 * r3) - (15 * p2 * q_s) - (400 * q_r_s) + (125 * p * s2);

    ExprPtr a2_C = (p2 * q4) - (6 * p3 * q2 * r) - (8 * q4 * r) + (9 * p4 * r2) + (76 * p_q2 * r2) - (136 * p2 * r3) + (400 * r4) - (50 * p * q3 * s) + (90 * p2 * q_r_s) - (1400 * q * r2 * s) + (625 * q2 * s2) + (500 * p * r * s2);

    ExprPtr a1_C = (-2 * p * q6) + (19 * p2 * q4 * r) - (51 * p3 * q2 * r2) + (3 * q4 * r2) + (32 * p4 * r3) + (76 * p_q2 * r3) - (256 * p2 * r4) + (512 * r5) - (31 * p3 * q3 * s) - (58 * q5 * s) + (117 * p4 * q_r_s) + (105 * p * q3 * r * s) + (260 * p2 * q * r2 * s) - (2400 * q * r3 * s) - (108 * p5 * s2) - (325 * p2 * q2 * s2) + (525 * p3 * r * s2) + (2750 * q2 * r * s2) - (500 * p * r2 * s2) + (625 * p * q * s3) - (3125 * s4);

    ExprPtr a0_C = (q8) - (13 * p * q6 * r) + (p5 * q2 * r2) + (65 * p2 * q4 * r2) - (4 * p6 * r3) - (128 * p3 * q2 * r3) + (17 * q4 * r3) + (48 * p4 * r4) - (16 * p_q2 * r4) - (192 * p2 * r5) + (256 * r6) - (4 * p5 * q3 * s) - (12 * p2 * q5 * s) + (18 * p6 * q_r_s) + (12 * p3 * q3 * r * s) - (124 * q5 * r * s) + (196 * p4 * q * r2 * s) + (590 * p * q3 * r2 * s) - (160 * p2 * q * r3 * s) - (1600 * q * r4 * s) - (27 * p7 * s2) - (150 * p4 * q2 * s2) - (125 * p * q4 * s2) - (99 * p5 * r * s2) - (725 * p2 * q2 * r * s2) + (1200 * p3 * r2 * s2) + (3250 * q2 * r2 * s2) - (2000 * p * r3 * s2) - (1250 * p * q * r * s3) + (3125 * p2 * s4) - (9375 * r * s4);

    ExprPtr a5_ = S(a5_C);
    ExprPtr a4_ = S(a4_C);
    ExprPtr a3_ = S(a3_C);
    ExprPtr a2_ = S(a2_C);
    ExprPtr a1_ = S(a1_C);
    ExprPtr a0_ = S(a0_C);

    Polynomial resolvent = Polynomial({INT_ONE, a5_, a4_, a3_, a2_, a1_, a0_}, 6);
    resolvent.printPolynomial();

    // ExprPtr normalizationFactor;
    // Polynomial resolvent_Ezouidi = Polynomial(resolvent, normalizationFactor, false);
    // resolvent_Ezouidi.printPolynomial();

    // vector<ExprPtr> rootsAry;
    // EMST::rootExtractor(resolvent, rootsAry, false, true);

    // -385 normalization

    ExprPtr a5 = INT_ZERO;
    ExprPtr a4 = intg(-2858625);
    ExprPtr a3 = intg(-471673125);
    ExprPtr a2 = ratg(8171736890625, 4);
    ExprPtr a1 = ratg(907062794859375, 2);
    ExprPtr a0 = ratg(-869677098584765625, 4);

    ExprPtr b2 = S(a5 / 2);
    ExprPtr b1 = S((a4 / 2) - ((a5 ^ 2) / 8));
    ExprPtr b0 = S((a3 / 2) + ((a5 ^ 3) / 16) - ((a4 * a5) / 4));

    ExprPtr c1sq = S((ratg(5, 64) * (a5 ^ 4)) - ((ratg(3, 8) * a4 * (a5 ^ 2))) + ((a4 ^ 2) / 4) + ((a3 * a5) / 2) - a2);

    ExprPtr c0sq = S((a3 * a4 * ratg(1, 4)) + (a4 * (a5 ^ 3) * ratg(1, 16)) - ((a4 ^ 2) * a5 * ratg(1, 8)) - (a3 * (a5 ^ 2) * ratg(1, 16)) - ((a5 ^ 5) * ratg(1, 128)) - (a1 * ratg(1, 2)));

    c1sq->printPreorder();
    c0sq->printPreorder();

    ExprPtr k11 = b2;
    ExprPtr k12 = S(b1 - powg(c1sq, FRACT_ONE_HALF, false));
    // ExprPtr k13 =

    // Polynomial cubic1 = Polynomial({k11, k12, k13}, 3);
}

int EMST_Tests()
{
    Polynomial p1({
                      INT_ONE,
                      ExprFactory::integer(-20),
                      ExprFactory::integer(190),
                      ExprFactory::integer(-1140),
                      ExprFactory::integer(4845),
                      ExprFactory::integer(-15504),
                      ExprFactory::integer(38760),
                      ExprFactory::integer(-77520),
                      ExprFactory::integer(125970),
                      ExprFactory::integer(-167960),
                      ExprFactory::integer(184756),
                      ExprFactory::integer(-167960),
                      ExprFactory::integer(125970),
                      ExprFactory::integer(-77520),
                      ExprFactory::integer(38760),
                      ExprFactory::integer(-15504),
                      ExprFactory::integer(4845),
                      ExprFactory::integer(-1140),
                      ExprFactory::integer(190),
                      ExprFactory::integer(-20),
                      INT_ONE,
                  },
                  20);

    Polynomial p2({INT_ONE,
                   ExprFactory::integer(18),
                   ExprFactory::integer(153),
                   ExprFactory::integer(816),
                   ExprFactory::integer(3060),
                   ExprFactory::integer(8568),
                   ExprFactory::integer(18564),
                   ExprFactory::integer(31824),
                   ExprFactory::integer(43758), INT_ONE,
                   ExprFactory::integer(-393813),
                   ExprFactory::integer(-1718460),
                   ExprFactory::integer(-4065432),
                   ExprFactory::integer(-6117426),
                   ExprFactory::integer(-6122934),
                   ExprFactory::integer(-4083180),
                   ExprFactory::integer(-1750131),
                   ExprFactory::integer(-437553),
                   ExprFactory::multiply({ExprFactory::integer(-5), ExprFactory::power(ExprFactory::integer(6), ExprFactory::rational(1, 2))})},
                  18);

    Polynomial p3({ExprFactory::integer(16),
                   INT_ZERO,
                   ExprFactory::integer(-48),
                   INT_ZERO,
                   ExprFactory::integer(12)},
                  4);

    Polynomial p4({INT_ONE,
                   ExprFactory::integer(2),
                   INT_ONE},
                  2);

    Polynomial p5({INT_ONE,
                   ExprFactory::integer(-12),
                   ExprFactory::rational(2185, 36),
                   ExprFactory::rational(-35760, 216),
                   ExprFactory::rational(332754, 1296),
                   ExprFactory::rational(-1668528, 7776),
                   ExprFactory::rational(3520800, 46656)},
                  6);

    // Note input should be l0, l1 values fix it.
    Polynomial p7({INT_ONE,
                   ExprFactory::integer(14),
                   ExprFactory::integer(84),
                   ExprFactory::integer(623),
                   ExprFactory::integer(3304),
                   ExprFactory::integer(8904),
                   ExprFactory::integer(11424),
                   ExprFactory::integer(5616)},
                  7);
    /*
        std::cout << "\nOriginal Polynomial :\n";

        p1.printPolynomial();

        std::cout << "\nApplied Ezouidi :\n";
    */
    try
    {
        // std::shared_ptr<Polynomial> ps2 = p7.normalize(1);
        // ps2->printPolynomial();

        RootValidator();
    }
    catch (std::runtime_error e)
    {
        std::cout << std::endl
                  << "ERROR" << std::endl;
        std::cout << e.what();
    }
    catch (std::exception e)
    {
        std::cout << std::endl
                  << "ERROR" << std::endl;

        std::cout << e.what();
    }

    // int R = 14; // 20sec

    int N = 7;

    vector<ExprPtr> skNodes{
        INT_ZERO, INT_ZERO, INT_ZERO, INT_ZERO, ExprFactory::integer(-343),
        INT_ZERO, INT_ZERO};

    // vector<ExprPtr> skNodes = {INT_ONE, INT_ONE, INT_ZERO};

    vector<ExprPtr> structuralRemainderNodes;

    structuralRemainderNodes.resize(skNodes.size());

    // std::cout << "\n --- H FUNCTION (r) : ---\n";

    // EMST::HFunction(R, skNodes)->printPreorder();

    // std::cout << "\n --- Delta FUNCTION (r) : ---\n";

    // EMST::DeltaFunction(7, 4, 4)->printPreorder(); // Delta(m , 4)

    // EMST::DiscriminantOperator(N, 1, skNodes)->printPreorder();

    // EMST::StructrualRemainder(N, skNodes, structuralRemainderNodes);

    // std::cout << "\n --- StructrualRemainder FUNCTION (N) : ---\n";

    // for (const auto &node : structuralRemainderNodes)
    // {
    // node->printPreorder();
    // }

    return 0;
}

int main()
{
    // simplifierTest();
    // RootExtractorDriver();
    // PolyTester();
    // PolyTester();
    quintic_formula();

    return 0;
}