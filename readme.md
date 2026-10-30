## C++ MathEngine

Some Features and Details of my Project

[(1) Classes](#(1)-Classes)

[(2) Programming Language Concepts Used](#(2)-Programming-Language-Concepts-Used)

[(3) Fast Exponentiation](#(3)-Fast-Exponentiation)

[(4) Prime Factorization Algorithm](#(4)-Prime-Factorization-Algorithm)

[(5) Multinomial Theorem Algorithm](#(5)-Multinomial-Theorem-Algorithm)

[(6) Expression Simplification Procedure](#(6)-Expression-Simplification-Procedure)

[(7) Test Cases](#(7)-Test-Cases)

---

### (1) Classes

<table>
  <tr>
    <th colspan="3">Expression</th>
  <tr>
    <td>Data Expr Nodes</td>
    <td>Operator Expr Nodes</td>
  </tr>
  <tr>
    <td>Integer</td>
    <td>Add</td>
  </tr>
  <tr>
    <td>Rational</td>
    <td>Multiply</td>
  </tr>
  <tr>
    <td>Symbol</td>
    <td>Power</td>
  </tr>
  <tr>
    <td>Constant</td>
    <td>Function</td>
  </tr>
  
</table>

### (2) Programming Language Concepts Used

1. Object Oriented Programming : Inheritance, Polymorphism, Encapsulation, Abstraction
2. C++ Smart Pointers
3. Error Handling
4. Templates
5. Input / Output Streams
6. Maps and Dictionaries
7. Big Integer Library - CppInt
8. Dynamic Memory Allocation

### (3) Fast Exponentiation

$$ a^b $$

1. Initialize `result = 1` and `e = |exp|`.
2. While `e > 0`:
    1. If `e` is odd, update `result = result * base`.
    2. Square the base: `base = base * base`.
    3. Divide `e` by `2` using a right shift.
3. Return `result`.

### (4) Prime Factorization Algorithm

An integer $N$ as a product of prime powers

$$ N = \prod_i p_i^q $$

Prime Factorization – Miller–Rabin + Parallel Pollard's Rho Algorithm

1. Take an integer `n` and handle special cases (`0`, negative numbers and `1`).
2. Recursively factorize `n`:
    1. Check whether `n` is probably prime using the Miller–Rabin test.
    2. If probably prime, append `n` to the factor list.
    3. Otherwise, find a non-trivial factor `d` using Parallel Pollard's Rho.
    4. Recursively factorize `d` and `n / d`.

3. Miller–Rabin Primality Test:
    1. Check divisibility by predefined small primes.
    2. Express `n - 1 = d × 2^s`, where `d` is odd.
    3. Perform modular exponentiation for each witness base.
    4. Use repeated squaring to check primality conditions.
    5. Return whether `n` is probably prime.

4. Pollard's Rho Factorization:
    1. Initialize random `x` and `c`, and take `y = x`.
    2. Define `f(x) = (x² + c) mod n`.
    3. While no divisor is obtained:
        1. Update `x = f(x)` [Tortoise Move].
        2. Update `y = f(f(y))` [Hare Move].
        3. Calculate `d = GCD(|x-y|, n)`.
    4. If `d = n`, restart with new random values.
    5. Otherwise, return the non-trivial factor `d`.

5. Parallel Pollard's Rho:
    1. Create multiple threads with independent random seeds.
    2. Run independent Pollard's Rho searches concurrently.
    3. Store the first non-trivial factor obtained.
    4. Stop further searches using an atomic flag.
    5. Join all threads and return the discovered factor.

6. Sort the obtained prime factors and combine identical factors by counting their exponents.

7. Print the prime factorization in the form:
   `n = p1^e1 * p2^e2 * ... * pk^ek`.


### (5) Multinomial Theorem Algorithm

The Algorithm Generates Expansion of
$$ (x_1 + x_2 + x_3 + ... + x_m)^n $$

Each Generated term has

$$ \frac{n!}{e_1!e_2!...e_m!} x_1^{e_1}x_2^{e_2}...x_m^{e_m} $$

where

$$ e_1 + e_2 + ... + e_m = n, e_i \ge 0 $$

Algorithm
1. Initialize the generator with degree `n` and a list of variables `vars`.
2. Initialize an exponent vector `exp` of size equal to the number of variables.
3. Call `dfs(0, n, 1, exp)` to generate all multinomial terms.
4. In the recursive DFS function:
    1. Calculate `m` as the total number of variables.
    2. If `pos` is the last variable:
        1. Assign `remaining` to `exp[pos]`.
        2. Generate the term using `printTerm(coeff, exp)`.
        3. Append the generated term to `mulNodes`.
        4. Return.
    3. Initialize `currentCoeff` equal to `coeff`.
    4. Assign `remaining` to `exp[pos]`.
    5. Recursively call `dfs(pos + 1, 0, currentCoeff, exp)`.
    6. For each `k` from `remaining - 1` down to `0`:
        1. Update the coefficient:
           `currentCoeff = currentCoeff * (k + 1) / (remaining - k)`.
        2. Assign `k` to `exp[pos]`.
        3. Recursively call `dfs(pos + 1, remaining - k, currentCoeff, exp)`.
5. In the `printTerm` function:
    1. Initialize an empty list of expression nodes.
    2. Add the coefficient as an integer expression.
    3. For each variable:
        1. If its exponent is `0`, skip the variable.
        2. If its exponent is `1`, add the variable directly.
        3. Otherwise, add the variable raised to its corresponding exponent.
    4. Multiply all expression nodes to construct the multinomial term.
6. After DFS completes:
    1. If `mulNodes` is empty, return integer `0`.
    2. If `mulNodes` contains exactly one term, return that term.
    3. Otherwise, return the sum of all terms in `mulNodes`.


### (6) Expression Simplification Procedure

Simplify Add Nodes

1. Remove Nested Add Nodes 
    
    $(a+ (b + c)) = a+b+c$
2. Perform addtion of all Integers and Rational in Expr
  
    $5 + \frac{1}{2} + x = \frac{11}{2}$

3. Group all coefficients of like terms and add them

    $2x^2 + \frac{7}{2}x^2+yx+xy=\frac{11}{2}x^2+2xy$

Simplify Multiply Nodes

1. Remove nested Multiply Nodes
   
   $(a)\cdot(bc\cdot(d)) = a\cdot b \cdot c\cdot d$
2. Multiply all Integers and Rationals to generate coefficients

    $2 \cdot 4 \cdot \frac{7}{2}xy^2 = 28xy^2$
3. Sum powers of all like symbols

    $x \cdot xy \cdot yz = x^2y^2z$
4. Distributative Law
    $(a+b)\cdot(c+d)=ac+ad+bc+bd$

Simplify Power Nodes

1. Addition of Exponents for same base and extraction of Integer/Rational of Power Nodes

    $2^3 \cdot \sqrt{4} = 16$

2. Multinomial Theorem Application for integer exponents

    $(x+2)^2=x^2+4x+4$

### (7) Test Cases

[1] $\sqrt{144}=12$

[2] $20+(-2)+\frac{1}{2}+(2i)+i=\frac{37}{2}+3i$

[3] $(1+2i)\cdot (1+3i)=-5+5i$

[4] $0^0 = {error}$

[5] $(\frac{12}{24})^{-2}=4$

[6] $i^4=1$

[7] $i^2=-1$

[8] $1^{2536} = 1$

[9] $144^{5/4} = 288\sqrt{3}$

[10] $(-2)^{-1/2}=\frac{1}{\sqrt{2}}i$

[11] $(1+\sqrt{2})^{(5/2)}=(3+2\sqrt{2})\sqrt{1+\sqrt{2}}$

[12] $\frac{3}{\sqrt{6}} \cdot 2= \sqrt{6}$