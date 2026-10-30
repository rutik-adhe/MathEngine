#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <cctype>
#include <stdexcept>

// Define all possible fundamental symbols in your language
enum class TokenType {
    IDENTIFIER, NUMBER, 
    PLUS, MINUS, STAR, SLASH, POWER, // Math operators
    EQUALS, COMMA,                   // Assignment & separators
    LPAREN, RPAREN, LBRACKET, RBRACKET, // Grouping & indexing
    PRINT, SIMPLIFY, SOLVE,          // Keywords
    END_OF_FILE
};

struct Token {
    TokenType type;
    std::string value; // Stores the actual string, e.g., "42", "x", "+"
};

class Lexer {
    std::string source;
    size_t pos = 0;

public:
    Lexer(const std::string& src) : source(src) {}

    Token getNextToken() {
        while (pos < source.length() && std::isspace(source[pos])) pos++;
        if (pos >= source.length()) return {TokenType::END_OF_FILE, ""};

        char current = source[pos];

        // Tokenise Numbers
        if (std::isdigit(current)) {
            std::string num;
            while (pos < source.length() && (std::isdigit(source[pos]) || source[pos] == '.')) {
                num += source[pos++];
            }
            return {TokenType::NUMBER, num};
        }

        // Tokenise Identifiers and Keywords
        if (std::isalpha(current)) {
            std::string id;
            while (pos < source.length() && std::isalnum(source[pos])) {
                id += source[pos++];
            }
            if (id == "print") return {TokenType::PRINT, id};
            if (id == "simplify") return {TokenType::SIMPLIFY, id};
            if (id == "solve") return {TokenType::SOLVE, id};
            return {TokenType::IDENTIFIER, id};
        }

        // Tokenise Single-character Operators
        pos++;
        switch (current) {
            case '=': return {TokenType::EQUALS, "="};
            case '+': return {TokenType::PLUS, "+"};
            case '-': return {TokenType::MINUS, "-"};
            case '*': return {TokenType::STAR, "*"};
            case '/': return {TokenType::SLASH, "/"};
            case '^': return {TokenType::POWER, "^"};
            case ',': return {TokenType::COMMA, ","};
            case '(': return {TokenType::LPAREN, "("};
            case ')': return {TokenType::RPAREN, ")"};
            case '[': return {TokenType::LBRACKET, "["};
            case ']': return {TokenType::RBRACKET, "]"};
        }
        throw std::runtime_error(std::string("Unknown character: ") + current);
    }
};

// Environment stores variable states between statements
using Environment = std::map<std::string, double>;

// Base AST Node
struct ASTNode {
    virtual ~ASTNode() = default;
    virtual double evaluate(Environment& env) = 0; 
};

// Leaf Node: A raw number (e.g., 3, 5)
struct NumberNode : public ASTNode {
    double value;
    NumberNode(double val) : value(val) {}
    double evaluate(Environment& env) override { return value; }
};

// Leaf Node: A variable lookup (e.g., a, b, x)
struct VariableNode : public ASTNode {
    std::string name;
    VariableNode(const std::string& n) : name(n) {}
    double evaluate(Environment& env) override {
        if (env.find(name) == env.end()) throw std::runtime_error("Undefined variable: " + name);
        return env[name];
    }
};

// Branch Node: Arithmetic operations (e.g., 3 + 5)
struct BinaryOpNode : public ASTNode {
    std::unique_ptr<ASTNode> left;
    Token op;
    std::unique_ptr<ASTNode> right;
    
    BinaryOpNode(std::unique_ptr<ASTNode> l, Token o, std::unique_ptr<ASTNode> r)
        : left(std::move(l)), op(o), right(std::move(r)) {}

    double evaluate(Environment& env) override {
        double l_val = left->evaluate(env);
        double r_val = right->evaluate(env);
        if (op.type == TokenType::PLUS) return l_val + r_val;
        // Handle MINUS, STAR, SLASH...
        return 0;
    }
};

// Branch Node: Function execution mapped to internal C++ logic
struct FunctionCallNode : public ASTNode {
    std::string func_name;
    std::vector<std::unique_ptr<ASTNode>> arguments;
    int array_index; // For [0] syntax in solve()
    
    FunctionCallNode(const std::string& name, int idx = -1) : func_name(name), array_index(idx) {}

    double evaluate(Environment& env) override {
        if (func_name == "simplify") {
            // Evaluate the argument node (e.g., 3 + 5 becomes 8)
            return arguments[0]->evaluate(env); 
        } 
        else if (func_name == "solve") {
            // Here you would serialise the argument AST back to a string and pass it 
            // to a symbolic engine (like SymEngine) to extract roots.
            // For now, simulate returning the first root of x^2 + 2x + 1
            std::vector<double> mock_roots = {-1.0, -1.0};
            return mock_roots[array_index];
        }
        throw std::runtime_error("Unknown function");
    }
};

// Statement Node: Assigns a value to the Environment
struct AssignmentNode : public ASTNode {
    std::string var_name;
    std::unique_ptr<ASTNode> expression;

    AssignmentNode(const std::string& name, std::unique_ptr<ASTNode> expr)
        : var_name(name), expression(std::move(expr)) {}

    double evaluate(Environment& env) override {
        double val = expression->evaluate(env);
        env[var_name] = val; // Store in state
        return val;
    }
};

class Parser {
    Lexer lexer;
    Token current_token;

    void advance() { current_token = lexer.getNextToken(); }
    void match(TokenType type) {
        if (current_token.type == type) advance();
        else throw std::runtime_error("Syntax error: Unexpected token " + current_token.value);
    }

public:
    Parser(const std::string& source) : lexer(source) {
        advance(); // Load first token
    }

    // Parses a full line of code
    std::unique_ptr<ASTNode> parseStatement() {
        if (current_token.type == TokenType::IDENTIFIER) {
            std::string var_name = current_token.value;
            advance();
            match(TokenType::EQUALS);
            
            // Parse what comes after the equals sign
            auto expr = parseExpression(); 
            return std::make_unique<AssignmentNode>(var_name, std::move(expr));
        }
        else if (current_token.type == TokenType::PRINT) {
            advance();
            match(TokenType::LPAREN);
            auto expr = parseExpression();
            match(TokenType::RPAREN);
            // PrintNode implementation omitted for brevity, but follows same pattern
            return nullptr; 
        }
        throw std::runtime_error("Unknown statement");
    }

    // Parses mathematical expressions (handles operator precedence)
    std::unique_ptr<ASTNode> parseExpression() {
        // A full implementation parses Terms and Factors to respect BEDMAS/BODMAS.
        // For this skeleton, we handle basic parsing.
        if (current_token.type == TokenType::SIMPLIFY) {
            auto func_node = std::make_unique<FunctionCallNode>("simplify");
            advance();
            match(TokenType::LPAREN);
            // In a complete parser, this would recursively call parseExpression()
            func_node->arguments.push_back(std::make_unique<NumberNode>(8)); 
            match(TokenType::RPAREN);
            return func_node;
        }
        return nullptr;
    }
};

int main() {
    std::string code_line = "a = simplify(3 + 5)";
    Environment global_env;

    try {
        Parser parser(code_line);
        std::unique_ptr<ASTNode> ast = parser.parseStatement();
        
        // Walking the tree evaluates the math and updates global_env
        ast->evaluate(global_env); 
        
        std::cout << "Variable 'a' stored as: " << global_env["a"] << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Compilation Error: " << e.what() << std::endl;
    }

    return 0;
}