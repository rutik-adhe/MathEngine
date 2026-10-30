#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_int.hpp>

using namespace std;
using boost::multiprecision::cpp_int;

static const uint32_t POW10[10] = {
    1u, 10u, 100u, 1000u, 10000u,
    100000u, 1000000u, 10000000u, 100000000u, 1000000000u
};

cpp_int parse_positive_decimal(const string& raw) {
    string s;
    s.reserve(raw.size());

    for (unsigned char c : raw) {
        if (!isspace(c)) s.push_back(char(c));
    }

    if (s.empty()) {
        throw runtime_error("Empty number");
    }

    size_t first_non_zero = s.find_first_not_of('0');
    if (first_non_zero == string::npos) {
        return 0;
    }
    s = s.substr(first_non_zero);

    cpp_int x = 0;

    for (size_t i = 0; i < s.size(); ) {
        size_t len = min<size_t>(9, s.size() - i);
        uint32_t chunk = 0;

        for (size_t j = 0; j < len; ++j) {
            char c = s[i + j];
            if (c < '0' || c > '9') {
                throw runtime_error("Invalid digit in input");
            }
            chunk = chunk * 10 + uint32_t(c - '0');
        }

        x *= POW10[len];
        x += chunk;
        i += len;
    }

    return x;
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string op, a_str, b_str;
    if (!(cin >> op >> a_str >> b_str)) {
        cerr << "Expected: <op> <a> <b>\n";
        return 1;
    }

    string out_file = "result.txt";
    if (argc >= 2) {
        out_file = argv[1];
    }

    try {
        cpp_int a = parse_positive_decimal(a_str);
        cpp_int b = parse_positive_decimal(b_str);
        cpp_int r;

        if (op == "+") {
            r = a + b;
        } else if (op == "-") {
            if (a < b) {
                throw runtime_error("Subtraction would become negative");
            }
            r = a - b;
        } else if (op == "*") {
            r = a * b;
        } else if (op == "/") {
            if (b == 0) {
                throw runtime_error("Division by zero");
            }
            r = a / b;   // integer quotient
        } else {
            throw runtime_error("Operation must be one of +, -, *, /");
        }

        ofstream out(out_file);
        if (!out) {
            throw runtime_error("Could not open output file");
        }

        out << r << '\n';
        cout << "Result written to " << out_file << '\n';
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << '\n';
        return 1;
    }

    return 0;
}