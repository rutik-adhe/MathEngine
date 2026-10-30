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

// Factory

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

// Function Pointers

ExprPtr (*intg)(cpp_int) = ExprFactory::integer;
ExprPtr (*ratg)(cpp_int, cpp_int) = ExprFactory::rational;
ExprPtr (*powg)(ExprPtr, ExprPtr, bool) = ExprFactory::power;
ExprPtr (*addg)(const vector<ExprPtr> &, bool) = ExprFactory::add;
ExprPtr (*mulg)(const vector<ExprPtr> &, bool) = ExprFactory::multiply;

// Commonly used Values

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

// Operator Overloading

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

// Simplifier

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

// Main

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

    // ExprPtr expr22 = EMST::HFunction(R, skNodes);

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

int main()
{
    // simplifierTest();
    // RootExtractorDriver();
    // PolyTester();
    // PolyTester();
    simplifierTest();

    return 0;
}