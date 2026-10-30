#include <fstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
#include <cstdint>

using std::string;
using std::vector;
using std::shared_ptr;
using std::make_shared;
using std::ifstream;
using std::ofstream;
using std::istream;
using std::ostream;

namespace ExprArchive
{
    static void writeU8(ostream& out, std::uint8_t v)
    {
        out.write(reinterpret_cast<const char*>(&v), sizeof(v));
    }

    static void writeU64(ostream& out, std::uint64_t v)
    {
        out.write(reinterpret_cast<const char*>(&v), sizeof(v));
    }

    static std::uint8_t readU8(istream& in)
    {
        std::uint8_t v = 0;
        in.read(reinterpret_cast<char*>(&v), sizeof(v));
        if (!in)
        {
            throw std::runtime_error("Unexpected end of file while reading u8");
        }
        return v;
    }

    static std::uint64_t readU64(istream& in)
    {
        std::uint64_t v = 0;
        in.read(reinterpret_cast<char*>(&v), sizeof(v));
        if (!in)
        {
            throw std::runtime_error("Unexpected end of file while reading u64");
        }
        return v;
    }

    static void writeString(ostream& out, const string& s)
    {
        std::uint64_t n = static_cast<std::uint64_t>(s.size());
        writeU64(out, n);
        if (n > 0)
        {
            out.write(s.data(), static_cast<std::streamsize>(n));
        }
    }

    static string readString(istream& in)
    {
        std::uint64_t n = readU64(in);
        string s;
        s.resize(static_cast<size_t>(n));
        if (n > 0)
        {
            in.read(&s[0], static_cast<std::streamsize>(n));
            if (!in)
            {
                throw std::runtime_error("Unexpected end of file while reading string");
            }
        }
        return s;
    }

    class CanonicalRegistry
    {
    public:
        static void registerCanonical(const string& name, const ExprPtr& expr)
        {
            if (!expr)
            {
                throw std::runtime_error("Cannot register null canonical expression: " + name);
            }

            auto& nm = byName();
            auto& ptr = byPtr();

            if (nm.find(name) != nm.end())
            {
                throw std::runtime_error("Duplicate canonical name: " + name);
            }

            nm.emplace(name, expr);
            ptr.emplace(expr.get(), name);
        }

        static const string* nameOf(const ExprPtr& expr)
        {
            if (!expr)
            {
                return nullptr;
            }

            auto& ptr = byPtr();
            auto it = ptr.find(expr.get());
            if (it == ptr.end())
            {
                return nullptr;
            }
            return &it->second;
        }

        static bool isCanonical(const ExprPtr& expr)
        {
            return nameOf(expr) != nullptr;
        }

        static ExprPtr resolve(const string& name)
        {
            auto& nm = byName();
            auto it = nm.find(name);
            if (it == nm.end())
            {
                throw std::runtime_error("Unknown canonical constant name in file: " + name);
            }
            return it->second;
        }

    private:
        static std::unordered_map<string, ExprPtr>& byName()
        {
            static std::unordered_map<string, ExprPtr> table;
            return table;
        }

        static std::unordered_map<const Expr*, string>& byPtr()
        {
            static std::unordered_map<const Expr*, string> table;
            return table;
        }
    };

    static std::size_t childCount(const ExprPtr& node)
    {
        switch (node->getKind())
        {
            case Expr::Kind::Add:
                return std::dynamic_pointer_cast<AddNode>(node)->children.size();

            case Expr::Kind::Multiply:
                return std::dynamic_pointer_cast<MultiplyNode>(node)->children.size();

            case Expr::Kind::Power:
                return 2;

            case Expr::Kind::Function:
                return std::dynamic_pointer_cast<FunctionNode>(node)->args.size();

            default:
                return 0;
        }
    }

    static ExprPtr childAt(const ExprPtr& node, std::size_t index)
    {
        switch (node->getKind())
        {
            case Expr::Kind::Add:
            {
                auto p = std::dynamic_pointer_cast<AddNode>(node);
                return p->children.at(index);
            }

            case Expr::Kind::Multiply:
            {
                auto p = std::dynamic_pointer_cast<MultiplyNode>(node);
                return p->children.at(index);
            }

            case Expr::Kind::Power:
            {
                auto p = std::dynamic_pointer_cast<PowerNode>(node);
                if (index == 0)
                {
                    return p->base;
                }
                if (index == 1)
                {
                    return p->exponent;
                }
                throw std::runtime_error("PowerNode child index out of range");
            }

            case Expr::Kind::Function:
            {
                auto p = std::dynamic_pointer_cast<FunctionNode>(node);
                return p->args.at(index);
            }

            default:
                throw std::runtime_error("Leaf node has no children");
        }
    }

    static void collectReachableNodes(
        const ExprPtr& root,
        vector<ExprPtr>& ordered,
        std::unordered_map<const Expr*, std::uint64_t>& idOf)
    {
        if (!root)
        {
            return;
        }

        struct Frame
        {
            ExprPtr node;
            std::size_t nextChild;
        };

        // 0 = unvisited, 1 = active, 2 = done
        std::unordered_map<const Expr*, std::uint8_t> state;
        vector<Frame> st;
        st.push_back({root, 0});

        while (!st.empty())
        {
            Frame& fr = st.back();
            const Expr* raw = fr.node.get();
            std::uint8_t& s = state[raw];

            if (s == 0)
            {
                s = 1;
                idOf[raw] = static_cast<std::uint64_t>(ordered.size());
                ordered.push_back(fr.node);

                // Canonical nodes are stored by name only; do not expand them.
                if (CanonicalRegistry::isCanonical(fr.node))
                {
                    s = 2;
                    st.pop_back();
                    continue;
                }
            }

            std::size_t cnt = childCount(fr.node);

            if (fr.nextChild < cnt)
            {
                ExprPtr ch = childAt(fr.node, fr.nextChild);

                if (!ch)
                {
                    throw std::runtime_error("Null child encountered in expression graph");
                }

                std::uint8_t childState = state[ch.get()];
                if (childState == 1)
                {
                    throw std::runtime_error("Cycle detected in expression graph");
                }

                ++fr.nextChild;

                if (childState == 0)
                {
                    st.push_back({ch, 0});
                }

                continue;
            }

            s = 2;
            st.pop_back();
        }
    }

    static void writeNodeRecord(
        ostream& out,
        const ExprPtr& node,
        const std::unordered_map<const Expr*, std::uint64_t>& idOf)
    {
        const string* canonName = CanonicalRegistry::nameOf(node);

        if (canonName)
        {
            writeU8(out, 0); // canonical record
            writeString(out, *canonName);
            return;
        }

        writeU8(out, 1); // concrete record
        writeU8(out, static_cast<std::uint8_t>(node->getKind()));
        writeU8(out, node->isSimplified() ? 1 : 0);

        switch (node->getKind())
        {
            case Expr::Kind::Integer:
            {
                auto p = std::dynamic_pointer_cast<IntegerNode>(node);
                writeString(out, p->value.str());
                break;
            }

            case Expr::Kind::Rational:
            {
                auto p = std::dynamic_pointer_cast<RationalNode>(node);
                writeString(out, p->value.str());
                break;
            }

            case Expr::Kind::Symbol:
            {
                auto p = std::dynamic_pointer_cast<SymbolNode>(node);
                writeString(out, p->name);
                break;
            }

            case Expr::Kind::Constant:
            {
                auto p = std::dynamic_pointer_cast<ConstantNode>(node);
                writeString(out, p->name);
                break;
            }

            case Expr::Kind::Complex:
            {
                auto p = std::dynamic_pointer_cast<ComplexNode>(node);
                writeString(out, p->value.str());
                break;
            }

            case Expr::Kind::Add:
            case Expr::Kind::Multiply:
            case Expr::Kind::Function:
            {
                std::size_t cnt = childCount(node);
                writeU64(out, static_cast<std::uint64_t>(cnt));
                for (std::size_t i = 0; i < cnt; ++i)
                {
                    ExprPtr ch = childAt(node, i);
                    auto it = idOf.find(ch.get());
                    if (it == idOf.end())
                    {
                        throw std::runtime_error("Internal error: child missing from ID map");
                    }
                    writeU64(out, it->second);
                }
                break;
            }

            case Expr::Kind::Power:
            {
                ExprPtr b = childAt(node, 0);
                ExprPtr e = childAt(node, 1);

                auto itb = idOf.find(b.get());
                auto ite = idOf.find(e.get());

                if (itb == idOf.end() || ite == idOf.end())
                {
                    throw std::runtime_error("Internal error: Power child missing from ID map");
                }

                writeU64(out, itb->second);
                writeU64(out, ite->second);
                break;
            }

            default:
                throw std::runtime_error("Unknown expression kind while saving");
        }
    }

    struct LoadedSlot
    {
        ExprPtr node;
        bool canonical = false;
        Expr::Kind kind = Expr::Kind::Integer;
        vector<std::uint64_t> childIds;
    };

    static void attachChildrenFromIds(
        const LoadedSlot& slot,
        vector<ExprPtr>& table)
    {
        if (slot.canonical)
        {
            return;
        }

        switch (slot.kind)
        {
            case Expr::Kind::Add:
            {
                auto p = std::dynamic_pointer_cast<AddNode>(slot.node);
                p->children.clear();
                for (std::uint64_t id : slot.childIds)
                {
                    if (id >= table.size())
                    {
                        throw std::runtime_error("Invalid child ID in file");
                    }
                    p->children.push_back(table[static_cast<size_t>(id)]);
                }
                break;
            }

            case Expr::Kind::Multiply:
            {
                auto p = std::dynamic_pointer_cast<MultiplyNode>(slot.node);
                p->children.clear();
                for (std::uint64_t id : slot.childIds)
                {
                    if (id >= table.size())
                    {
                        throw std::runtime_error("Invalid child ID in file");
                    }
                    p->children.push_back(table[static_cast<size_t>(id)]);
                }
                break;
            }

            case Expr::Kind::Power:
            {
                auto p = std::dynamic_pointer_cast<PowerNode>(slot.node);
                if (slot.childIds.size() != 2)
                {
                    throw std::runtime_error("Power node must have exactly 2 children");
                }
                p->base = table[static_cast<size_t>(slot.childIds[0])];
                p->exponent = table[static_cast<size_t>(slot.childIds[1])];
                break;
            }

            case Expr::Kind::Function:
            {
                auto p = std::dynamic_pointer_cast<FunctionNode>(slot.node);
                p->args.clear();
                for (std::uint64_t id : slot.childIds)
                {
                    if (id >= table.size())
                    {
                        throw std::runtime_error("Invalid child ID in file");
                    }
                    p->args.push_back(table[static_cast<size_t>(id)]);
                }
                break;
            }

            default:
                break;
        }
    }

    void saveExprToFile(const ExprPtr& root, const string& fileName)
    {
        ofstream out(fileName, std::ios::binary);
        if (!out)
        {
            throw std::runtime_error("Cannot open file for writing: " + fileName);
        }

        // Header
        out.write("EXPR", 4);
        writeU8(out, 2); // version 2 = graph + canonical registry

        if (!root)
        {
            writeU8(out, 0);
            return;
        }

        writeU8(out, 1);

        vector<ExprPtr> ordered;
        std::unordered_map<const Expr*, std::uint64_t> idOf;
        collectReachableNodes(root, ordered, idOf);

        auto itRoot = idOf.find(root.get());
        if (itRoot == idOf.end())
        {
            throw std::runtime_error("Root was not assigned an ID");
        }

        writeU64(out, itRoot->second);
        writeU64(out, static_cast<std::uint64_t>(ordered.size()));

        for (const auto& node : ordered)
        {
            writeNodeRecord(out, node, idOf);
        }
    }

    ExprPtr loadExprFromFile(const string& fileName)
    {
        ifstream in(fileName, std::ios::binary);
        if (!in)
        {
            throw std::runtime_error("Cannot open file for reading: " + fileName);
        }

        char magic[4];
        in.read(magic, 4);
        if (!in || magic[0] != 'E' || magic[1] != 'X' || magic[2] != 'P' || magic[3] != 'R')
        {
            throw std::runtime_error("Invalid file format");
        }

        std::uint8_t version = readU8(in);
        if (version != 2)
        {
            throw std::runtime_error("Unsupported Expr file version");
        }

        std::uint8_t hasRoot = readU8(in);
        if (hasRoot == 0)
        {
            return nullptr;
        }

        std::uint64_t rootId = readU64(in);
        std::uint64_t count = readU64(in);

        if (count == 0)
        {
            throw std::runtime_error("Corrupt file: non-null root but zero nodes");
        }

        vector<LoadedSlot> slots(static_cast<size_t>(count));
        vector<ExprPtr> table(static_cast<size_t>(count));

        for (std::uint64_t i = 0; i < count; ++i)
        {
            std::uint8_t recordType = readU8(in);

            if (recordType == 0)
            {
                string name = readString(in);
                ExprPtr expr = CanonicalRegistry::resolve(name);
                slots[static_cast<size_t>(i)].node = expr;
                slots[static_cast<size_t>(i)].canonical = true;
                table[static_cast<size_t>(i)] = expr;
                continue;
            }

            if (recordType != 1)
            {
                throw std::runtime_error("Unknown record type in file");
            }

            std::uint8_t kindByte = readU8(in);
            std::uint8_t simplifiedByte = readU8(in);

            Expr::Kind kind = static_cast<Expr::Kind>(kindByte);
            bool simplified = (simplifiedByte != 0);

            LoadedSlot slot;
            slot.kind = kind;

            switch (kind)
            {
                case Expr::Kind::Integer:
                {
                    string s = readString(in);
                    cpp_int v(s);
                    slot.node = make_shared<IntegerNode>(v);
                    break;
                }

                case Expr::Kind::Rational:
                {
                    string s = readString(in);

                    // Adjust this line if your Rational type uses another parser.
                    Rational r(s);

                    slot.node = make_shared<RationalNode>(r);
                    break;
                }

                case Expr::Kind::Symbol:
                {
                    string name = readString(in);
                    slot.node = make_shared<SymbolNode>(name);
                    break;
                }

                case Expr::Kind::Constant:
                {
                    string name = readString(in);
                    slot.node = make_shared<ConstantNode>(name);
                    break;
                }

                case Expr::Kind::Complex:
                {
                    string s = readString(in);

                    // Adjust this line if your ComplexExact type uses another parser.
                    ComplexExact c(s);

                    slot.node = make_shared<ComplexNode>(c);
                    break;
                }

                case Expr::Kind::Add:
                case Expr::Kind::Multiply:
                case Expr::Kind::Function:
                {
                    std::uint64_t n = readU64(in);
                    slot.childIds.resize(static_cast<size_t>(n));
                    for (std::uint64_t j = 0; j < n; ++j)
                    {
                        slot.childIds[static_cast<size_t>(j)] = readU64(in);
                    }

                    if (kind == Expr::Kind::Add)
                    {
                        vector<ExprPtr> empty;
                        slot.node = make_shared<AddNode>(std::move(empty), simplified);
                    }
                    else if (kind == Expr::Kind::Multiply)
                    {
                        vector<ExprPtr> empty;
                        slot.node = make_shared<MultiplyNode>(std::move(empty), simplified);
                    }
                    else
                    {
                        string fnName = ""; // read above? not yet.
                        // Function needs its name before child count in the save format.
                        // This branch should not happen because Function is handled separately below.
                        throw std::runtime_error("Internal format mismatch for Function node");
                    }

                    break;
                }

                case Expr::Kind::Power:
                {
                    std::uint64_t baseId = readU64(in);
                    std::uint64_t expId = readU64(in);
                    slot.childIds.push_back(baseId);
                    slot.childIds.push_back(expId);

                    slot.node = make_shared<PowerNode>(nullptr, nullptr, simplified);
                    break;
                }

                default:
                    throw std::runtime_error("Unknown node kind while loading");
            }

            // Function node must be handled separately because it has a name.
            if (kind == Expr::Kind::Function)
            {
                string fnName = readString(in);
                std::uint64_t n = readU64(in);
                slot.childIds.resize(static_cast<size_t>(n));
                for (std::uint64_t j = 0; j < n; ++j)
                {
                    slot.childIds[static_cast<size_t>(j)] = readU64(in);
                }

                vector<ExprPtr> empty;
                slot.node = make_shared<FunctionNode>(fnName, std::move(empty));
            }

            slots[static_cast<size_t>(i)] = std::move(slot);
            table[static_cast<size_t>(i)] = slots[static_cast<size_t>(i)].node;
        }

        for (std::uint64_t i = 0; i < count; ++i)
        {
            attachChildrenFromIds(slots[static_cast<size_t>(i)], table);
        }

        if (rootId >= count)
        {
            throw std::runtime_error("Invalid root ID in file");
        }

        return table[static_cast<size_t>(rootId)];
    }
}