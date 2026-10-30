struct Frame
{
    int pos;
    int remaining;
    unsigned long long coeff;
    std::vector<int> exp;
};

class MultinomialGenerator
{
    int n;
    std::vector<ExprPtr> vars;

public:
    MultinomialGenerator(
        int degree,
        std::vector<ExprPtr> variables)
        : n(degree),
          vars(std::move(variables))
    {
    }

    ExprPtr generate()
    {
        mulNodes.clear();

        if (vars.empty())
            return ExprFactory::integer(0);

        std::stack<Frame> st;

        st.push({0,
                 n,
                 1ULL,
                 std::vector<int>(vars.size(), 0)});

        while (!st.empty())
        {
            Frame cur = std::move(st.top());
            st.pop();

            int m = vars.size();

            if (cur.pos == m - 1)
            {
                cur.exp[cur.pos] = cur.remaining;

                mulNodes.push_back(
                    printTerm(
                        cur.coeff,
                        cur.exp));
                continue;
            }

            unsigned long long currentCoeff =
                cur.coeff;

            auto expCopy = cur.exp;

            expCopy[cur.pos] =
                cur.remaining;

            st.push({cur.pos + 1,
                     0,
                     currentCoeff,
                     expCopy});

            for (int k = cur.remaining - 1;
                 k >= 0;
                 --k)
            {
                currentCoeff =
                    currentCoeff *
                    (unsigned long long)(k + 1) /
                    (unsigned long long)(cur.remaining - k);

                expCopy = cur.exp;

                expCopy[cur.pos] = k;

                st.push({cur.pos + 1,
                         cur.remaining - k,
                         currentCoeff,
                         expCopy});
            }
        }

        if (mulNodes.empty())
            return ExprFactory::integer(0);

        if (mulNodes.size() == 1)
            return mulNodes[0];

        return ExprFactory::add(mulNodes);
    }

private:
    std::vector<ExprPtr> mulNodes;

    ExprPtr printTerm(
        unsigned long long coeff,
        const std::vector<int> &exp)
    {
        std::vector<ExprPtr> nodes;

        nodes.push_back(
            ExprFactory::integer(coeff));

        for (size_t i = 0;
             i < exp.size();
             ++i)
        {
            if (exp[i] == 0)
                continue;

            if (exp[i] == 1)
            {
                nodes.push_back(vars[i]);
            }
            else
            {
                nodes.push_back(
                    ExprFactory::power(
                        vars[i],
                        ExprFactory::integer(exp[i])));
            }
        }

        return ExprFactory::multiply(nodes);
    }
};