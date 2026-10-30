#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
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

class Rational {
public:
    long long n; // numerator
    long long d; // denominator, always positive

    Rational(long long num = 0, long long den = 1) : n(num), d(den) {
        normalize();
    }

    void normalize() {
        if (d == 0) {
            throw std::runtime_error("Rational denominator cannot be zero");
        }

        if (d < 0) {
            n = -n;
            d = -d;
        }

        long long g = std::gcd(std::llabs(n), std::llabs(d));
        if (g == 0) {
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

    string str() const {
        if (d == 1) return std::to_string(n);
        return std::to_string(n) + "/" + std::to_string(d);
    }

    friend bool operator==(const Rational& a, const Rational& b) {
        return a.n == b.n && a.d == b.d;
    }

    friend bool operator!=(const Rational& a, const Rational& b) {
        return !(a == b);
    }

    friend Rational operator+(const Rational& a, const Rational& b) {
        return Rational(a.n * b.d + b.n * a.d, a.d * b.d);
    }

    friend Rational operator-(const Rational& a, const Rational& b) {
        return Rational(a.n * b.d - b.n * a.d, a.d * b.d);
    }

    friend Rational operator*(const Rational& a, const Rational& b) {
        return Rational(a.n * b.n, a.d * b.d);
    }

    friend Rational operator/(const Rational& a, const Rational& b) {
        if (b.n == 0) throw std::runtime_error("Division by zero rational");
        return Rational(a.n * b.d, a.d * b.n);
    }
};

static Rational powRationalInt(const Rational& base, long long exp) {
    if (exp == 0) return Rational(1);
    if (base.isZero()) {
        if (exp < 0) throw std::runtime_error("0 cannot be raised to negative power");
        return Rational(0);
    }

    Rational result(1);
    Rational b = base;
    long long e = std::llabs(exp);

    while (e > 0) {
        if (e & 1LL) result = result * b;
        b = b * b;
        e >>= 1LL;
    }

    if (exp < 0) {
        return Rational(result.d, result.n);
    }
    return result;
}

// exact integer nth root for perfect powers
static bool exactNthRoot(long long x, int n, long long& root) {
    if (n <= 0) return false;
    if (x == 0) {
        root = 0;
        return true;
    }

    bool neg = x < 0;
    long long ax = std::llabs(x);

    if (neg && (n % 2 == 0)) return false;

    long long lo = 0, hi = ax;
    while (lo <= hi) {
        long long mid = lo + (hi - lo) / 2;

        __int128 p = 1;
        for (int i = 0; i < n; ++i) {
            p *= mid;
            if (p > ax) break;
        }

        if (p == ax) {
            root = neg ? -mid : mid;
            return true;
        } else if (p < ax) {
            lo = mid + 1;
        } else {
            hi = mid - 1;
        }
    }

    return false;
}

// ============================================================
// Expression hierarchy
// ============================================================

class Expr {
public:
    enum class Kind {
        Integer,
        Rational,
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

class IntegerNode : public Expr {
public:
    long long value;

    explicit IntegerNode(long long v) : Expr(Kind::Integer), value(v) {}

    shared_ptr<Expr> clone() const override {
        return make_shared<IntegerNode>(value);
    }

    string toString() const override {
        return std::to_string(value);
    }

    void printPreorder(int indent = 0) const override {
        std::cout << string(indent, ' ') << "Integer(" << value << ")\n";
    }
};

class RationalNode : public Expr {
public:
    Rational value;

    explicit RationalNode(const Rational& r) : Expr(Kind::Rational), value(r) {}

    shared_ptr<Expr> clone() const override {
        return make_shared<RationalNode>(value);
    }

    string toString() const override {
        return value.str();
    }

    void printPreorder(int indent = 0) const override {
        std::cout << string(indent, ' ') << "Rational(" << value.str() << ")\n";
    }
};

class SymbolNode : public Expr {
public:
    string name;

    explicit SymbolNode(string n) : Expr(Kind::Symbol), name(std::move(n)) {}

    shared_ptr<Expr> clone() const override {
        return make_shared<SymbolNode>(name);
    }

    string toString() const override {
        return name;
    }

    void printPreorder(int indent = 0) const override {
        std::cout << string(indent, ' ') << "Symbol(" << name << ")\n";
    }
};

class ConstantNode : public Expr {
public:
    string name;

    explicit ConstantNode(string n) : Expr(Kind::Constant), name(std::move(n)) {}

    shared_ptr<Expr> clone() const override {
        return make_shared<ConstantNode>(name);
    }

    string toString() const override {
        return name;
    }

    void printPreorder(int indent = 0) const override {
        std::cout << string(indent, ' ') << "Constant(" << name << ")\n";
    }
};

class AddNode : public Expr {
public:
    vector<ExprPtr> children;

    explicit AddNode(vector<ExprPtr> c) : Expr(Kind::Add), children(std::move(c)) {}

    shared_ptr<Expr> clone() const override {
        vector<ExprPtr> out;
        out.reserve(children.size());
        for (const auto& c : children) out.push_back(c->clone());
        return make_shared<AddNode>(out);
    }

    string toString() const override {
        vector<string> parts;
        for (const auto& c : children) parts.push_back(c->toString());
        std::ostringstream oss;
        oss << "Add(";
        for (size_t i = 0; i < parts.size(); ++i) {
            if (i) oss << ",";
            oss << parts[i];
        }
        oss << ")";
        return oss.str();
    }

    void printPreorder(int indent = 0) const override {
        std::cout << string(indent, ' ') << "Add\n";
        for (const auto& c : children) c->printPreorder(indent + 2);
    }
};

class MultiplyNode : public Expr {
public:
    vector<ExprPtr> children;

    explicit MultiplyNode(vector<ExprPtr> c) : Expr(Kind::Multiply), children(std::move(c)) {}

    shared_ptr<Expr> clone() const override {
        vector<ExprPtr> out;
        out.reserve(children.size());
        for (const auto& c : children) out.push_back(c->clone());
        return make_shared<MultiplyNode>(out);
    }

    string toString() const override {
        vector<string> parts;
        for (const auto& c : children) parts.push_back(c->toString());
        std::ostringstream oss;
        oss << "Mul(";
        for (size_t i = 0; i < parts.size(); ++i) {
            if (i) oss << ",";
            oss << parts[i];
        }
        oss << ")";
        return oss.str();
    }

    void printPreorder(int indent = 0) const override {
        std::cout << string(indent, ' ') << "Multiply\n";
        for (const auto& c : children) c->printPreorder(indent + 2);
    }
};

class PowerNode : public Expr {
public:
    ExprPtr base;
    ExprPtr exponent;

    PowerNode(ExprPtr b, ExprPtr e) : Expr(Kind::Power), base(std::move(b)), exponent(std::move(e)) {}

    shared_ptr<Expr> clone() const override {
        return make_shared<PowerNode>(base->clone(), exponent->clone());
    }

    string toString() const override {
        return "Pow(" + base->toString() + "," + exponent->toString() + ")";
    }

    void printPreorder(int indent = 0) const override {
        std::cout << string(indent, ' ') << "Power\n";
        base->printPreorder(indent + 2);
        exponent->printPreorder(indent + 2);
    }
};

class FunctionNode : public Expr {
public:
    string name;
    vector<ExprPtr> args;

    FunctionNode(string n, vector<ExprPtr> a) : Expr(Kind::Function), name(std::move(n)), args(std::move(a)) {}

    shared_ptr<Expr> clone() const override {
        vector<ExprPtr> out;
        out.reserve(args.size());
        for (const auto& a : args) out.push_back(a->clone());
        return make_shared<FunctionNode>(name, out);
    }

    string toString() const override {
        std::ostringstream oss;
        oss << "Fn(" << name << "(";
        for (size_t i = 0; i < args.size(); ++i) {
            if (i) oss << ",";
            oss << args[i]->toString();
        }
        oss << "))";
        return oss.str();
    }

    void printPreorder(int indent = 0) const override {
        std::cout << string(indent, ' ') << "Function(" << name << ")\n";
        for (const auto& a : args) a->printPreorder(indent + 2);
    }
};

// ============================================================
// Linear factory helpers
// ============================================================

class ExprFactory {
public:
    static ExprPtr integer(long long v) {
        return make_shared<IntegerNode>(v);
    }

    static ExprPtr rational(long long n, long long d) {
        Rational r(n, d);
        if (r.d == 1) return integer(r.n);
        return make_shared<RationalNode>(r);
    }

    static ExprPtr rational(const Rational& r) {
        if (r.d == 1) return integer(r.n);
        return make_shared<RationalNode>(r);
    }

    static ExprPtr symbol(const string& name) {
        return make_shared<SymbolNode>(name);
    }

    static ExprPtr constant(const string& name) {
        return make_shared<ConstantNode>(name);
    }

    static ExprPtr add(const vector<ExprPtr>& children) {
        if (children.size() == 1) return children[0];
        return make_shared<AddNode>(children);
    }

    static ExprPtr multiply(const vector<ExprPtr>& children) {
        if (children.size() == 1) return children[0];
        return make_shared<MultiplyNode>(children);
    }

    static ExprPtr power(ExprPtr base, ExprPtr exponent) {
        return make_shared<PowerNode>(std::move(base), std::move(exponent));
    }

    static ExprPtr function(const string& name, const vector<ExprPtr>& args) {
        return make_shared<FunctionNode>(name, args);
    }
};

// ============================================================
// Simplifier
// ============================================================

class Simplifier {
public:
    static ExprPtr simplify(const ExprPtr& e) {
        if (!e) return nullptr;

        switch (e->getKind()) {
            case Expr::Kind::Integer:
            case Expr::Kind::Rational:
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
    static bool isNumeric(const ExprPtr& e, Rational& out) {
        if (!e) return false;

        if (e->getKind() == Expr::Kind::Integer) {
            out = Rational(std::static_pointer_cast<IntegerNode>(e)->value);
            return true;
        }

        if (e->getKind() == Expr::Kind::Rational) {
            out = std::static_pointer_cast<RationalNode>(e)->value;
            return true;
        }

        return false;
    }

    static ExprPtr numericToExpr(const Rational& r) {
        return ExprFactory::rational(r);
    }

    static string keyOf(const ExprPtr& e) {
        return e ? e->toString() : "null";
    }

    static void sortCanonical(vector<ExprPtr>& items) {
        std::sort(items.begin(), items.end(),
                  [](const ExprPtr& a, const ExprPtr& b) {
                      return a->toString() < b->toString();
                  });
    }

    static std::pair<Rational, ExprPtr> splitCoefficient(const ExprPtr& term) {
        Rational r;
        if (isNumeric(term, r)) {
            return {r, nullptr};
        }

        if (term->getKind() == Expr::Kind::Multiply) {
            auto m = std::static_pointer_cast<MultiplyNode>(term);
            Rational coeff(1);
            vector<ExprPtr> rest;

            for (const auto& f : m->children) {
                Rational fr;
                if (isNumeric(f, fr)) {
                    coeff = coeff * fr;
                } else {
                    rest.push_back(f);
                }
            }

            if (rest.empty()) return {coeff, nullptr};
            if (rest.size() == 1) return {coeff, rest[0]};
            return {coeff, ExprFactory::multiply(rest)};
        }

        return {Rational(1), term};
    }

    static ExprPtr simplifyFunction(const std::shared_ptr<FunctionNode>& f) {
        vector<ExprPtr> args;
        args.reserve(f->args.size());
        for (const auto& a : f->args) args.push_back(simplify(a));
        return make_shared<FunctionNode>(f->name, args);
    }

    static ExprPtr simplifyPower(const std::shared_ptr<PowerNode>& p) {
        ExprPtr base = simplify(p->base);
        ExprPtr exp  = simplify(p->exponent);

        Rational expNum;
        bool expIsNum = isNumeric(exp, expNum);

        // x^1 = x
        if (expIsNum && expNum.isOne()) {
            return base;
        }

        // x^0 = 1  (except 0^0, left as-is)
        if (expIsNum && expNum.isZero()) {
            Rational baseNum;
            if (!(isNumeric(base, baseNum) && baseNum.isZero())) {
                return ExprFactory::integer(1);
            }
        }

        // 1^x = 1
        Rational baseNum;
        if (isNumeric(base, baseNum)) {
            if (baseNum.isOne()) return ExprFactory::integer(1);
            if (baseNum.isZero()) {
                if (expIsNum && expNum.n > 0) return ExprFactory::integer(0);
            }
        }

        // numeric power evaluation
        if (expIsNum) {
            Rational b;
            if (isNumeric(base, b)) {
                if (expNum.isInteger()) {
                    return numericToExpr(powRationalInt(b, expNum.n));
                }

                // exact nth roots for perfect powers
                if (std::llabs(expNum.n) == 1 && expNum.d > 1) {
                    int n = static_cast<int>(expNum.d);
                    long long root = 0;

                    if (b.d == 1) {
                        if (exactNthRoot(b.n, n, root)) {
                            if (expNum.n == 1) {
                                return ExprFactory::integer(root);
                            } else {
                                return ExprFactory::rational(Rational(1, root));
                            }
                        }
                    }
                }
            }
        }

        // (a^m)^n => a^(m*n)
        if (base->getKind() == Expr::Kind::Power) {
            auto inner = std::static_pointer_cast<PowerNode>(base);
            Rational innerExp;
            if (isNumeric(inner->exponent, innerExp) && expIsNum) {
                Rational combined = innerExp * expNum;
                return simplify(ExprFactory::power(inner->base, numericToExpr(combined)));
            }
        }

        return ExprFactory::power(base, exp);
    }

    static ExprPtr simplifyMultiply(const std::shared_ptr<MultiplyNode>& m) {
        // Step 1: simplify children and flatten nested multiplies
        vector<ExprPtr> factors;
        for (const auto& c : m->children) {
            ExprPtr s = simplify(c);
            if (s->getKind() == Expr::Kind::Multiply) {
                auto mm = std::static_pointer_cast<MultiplyNode>(s);
                for (const auto& k : mm->children) factors.push_back(k);
            } else {
                factors.push_back(s);
            }
        }

        // Step 2: combine coefficient and repeated powers
        Rational coeff(1);
        std::map<string, std::pair<ExprPtr, Rational>> powerMap;

        for (const auto& f : factors) {
            Rational r;
            if (isNumeric(f, r)) {
                coeff = coeff * r;
                continue;
            }

            if (f->getKind() == Expr::Kind::Power) {
                auto p = std::static_pointer_cast<PowerNode>(f);
                Rational expNum;
                if (isNumeric(p->exponent, expNum)) {
                    string key = keyOf(p->base);
                    auto it = powerMap.find(key);
                    if (it == powerMap.end()) {
                        powerMap[key] = {p->base, expNum};
                    } else {
                        it->second.second = it->second.second + expNum;
                    }
                    continue;
                }
            }

            string key = keyOf(f);
            auto it = powerMap.find(key);
            if (it == powerMap.end()) {
                powerMap[key] = {f, Rational(1)};
            } else {
                it->second.second = it->second.second + Rational(1);
            }
        }

        // zero annihilates product
        if (coeff.isZero()) return ExprFactory::integer(0);

        vector<ExprPtr> out;

        // rebuild powers
        for (auto& kv : powerMap) {
            ExprPtr base = kv.second.first;
            Rational exp = kv.second.second;

            if (exp.isZero()) continue;

            if (exp.isOne()) {
                out.push_back(base);
            } else {
                out.push_back(simplify(ExprFactory::power(base, numericToExpr(exp))));
            }
        }

        sortCanonical(out);

        // attach coefficient
        if (!coeff.isOne() || out.empty()) {
            if (out.empty()) {
                return numericToExpr(coeff);
            }
            out.insert(out.begin(), numericToExpr(coeff));
        }

        if (out.size() == 1) return out[0];
        return ExprFactory::multiply(out);
    }

    static ExprPtr simplifyAdd(const std::shared_ptr<AddNode>& a) {
        // Step 1: simplify children and flatten nested adds
        vector<ExprPtr> terms;
        for (const auto& c : a->children) {
            ExprPtr s = simplify(c);
            if (s->getKind() == Expr::Kind::Add) {
                auto aa = std::static_pointer_cast<AddNode>(s);
                for (const auto& k : aa->children) terms.push_back(k);
            } else {
                terms.push_back(s);
            }
        }

        // Step 2: combine numeric terms and like terms
        Rational numericSum(0);
        std::map<string, std::pair<Rational, ExprPtr>> groups;

        for (const auto& t : terms) {
            Rational r;
            if (isNumeric(t, r)) {
                numericSum = numericSum + r;
                continue;
            }

            auto [coeff, rest] = splitCoefficient(t);

            if (!rest) {
                numericSum = numericSum + coeff;
                continue;
            }

            string key = keyOf(rest);
            auto it = groups.find(key);
            if (it == groups.end()) {
                groups[key] = {coeff, rest};
            } else {
                it->second.first = it->second.first + coeff;
            }
        }

        // Step 3: rebuild
        vector<ExprPtr> out;

        if (!numericSum.isZero()) {
            out.push_back(numericToExpr(numericSum));
        }

        for (auto& kv : groups) {
            Rational coeff = kv.second.first;
            ExprPtr rest = kv.second.second;

            if (coeff.isZero()) continue;

            if (coeff.isOne()) {
                out.push_back(rest);
            } else {
                out.push_back(simplify(ExprFactory::multiply({numericToExpr(coeff), rest})));
            }
        }

        sortCanonical(out);

        if (out.empty()) return ExprFactory::integer(0);
        if (out.size() == 1) return out[0];
        return ExprFactory::add(out);
    }
};

// ============================================================
// Demo / main
// ============================================================

int main() {
    // Build:
    // (3 + 2i)^2 + (x + x) + (6*x*y) + (2*x*y) + sqrt(9)
    //
    // sqrt(9) is written as 9^(1/2)

    ExprPtr i = ExprFactory::constant("i");
    ExprPtr x = ExprFactory::symbol("x");
    ExprPtr y = ExprFactory::symbol("y");

    ExprPtr complexPart = ExprFactory::power(
        ExprFactory::add({
            ExprFactory::integer(3),
            ExprFactory::multiply({ExprFactory::integer(2), i})
        }),
        ExprFactory::integer(2)
    );

    ExprPtr xPlusX = ExprFactory::add({x, x});
    ExprPtr sixxy = ExprFactory::multiply({ExprFactory::integer(6), x, y});
    ExprPtr twoxy = ExprFactory::multiply({ExprFactory::integer(2), x, y});
    ExprPtr sqrt9 = ExprFactory::power(ExprFactory::integer(9), ExprFactory::rational(1, 2));

    ExprPtr expr = ExprFactory::add({complexPart, xPlusX, sixxy, twoxy, sqrt9});

    std::cout << "Original tree (preorder):\n";
    expr->printPreorder();

    ExprPtr simp = Simplifier::simplify(expr);

    std::cout << "\nSimplified tree (preorder):\n";
    simp->printPreorder();

    std::cout << "\nCanonical string form:\n";
    std::cout << simp->toString() << "\n";

    return 0;
}