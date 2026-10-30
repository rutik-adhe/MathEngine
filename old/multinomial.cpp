#include <iostream>
#include <vector>
#include <string>

using namespace std;

class MultinomialGenerator
{
    int n;
    vector<string> vars;

public:
    MultinomialGenerator(int degree,
                         vector<string> variables)
        : n(degree),
          vars(move(variables))
    {
    }

    void generate()
    {
        vector<int> exp(vars.size());

        dfs(0, n, 1ULL, exp);
    }

private:
    void dfs(int pos,
             int remaining,
             unsigned long long coeff,
             vector<int>& exp)
    {
        // m is total number of variables
        int m = vars.size();

        if(pos == m - 1)
        {
            exp[pos] = remaining;

            printTerm(coeff, exp);
            return;
        }

        unsigned long long currentCoeff = coeff;

        exp[pos] = remaining;

        dfs(pos + 1,
            0,
            currentCoeff,
            exp);

        for(int k = remaining - 1; k >= 0; --k)
        {
            currentCoeff =
                currentCoeff *
                (unsigned long long)(k + 1)
                /
                (unsigned long long)(remaining - k);

            exp[pos] = k;

            dfs(pos + 1,
                remaining - k,
                currentCoeff,
                exp);
        }
    }

    void printTerm(unsigned long long coeff,
                   const vector<int>& exp)
    {
        cout << coeff;

        for(size_t i = 0; i < exp.size(); ++i)
        {
            if(exp[i] == 0)
                continue;

            cout << "*" << vars[i];

            if(exp[i] > 1)
                cout << "^" << exp[i];
        }

        cout << '\n';
    }
};

int main()
{
    MultinomialGenerator gen(
        2,
        {"a","b"}
    );

    gen.generate();
}