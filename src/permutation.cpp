#include "permutation.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <set>
#include <stdexcept>

namespace
{
    // Orders labels by absolute value, with -k before k.
    bool LabelLess(int a, int b)
    {
        const int absA = std::abs(a);
        const int absB = std::abs(b);
        if (absA != absB) {
            return absA < absB;
        }
        return a < b;
    }
}

Permutation Permutation::FromCycles(const std::vector<std::vector<int>>& cycles)
{
    Permutation p;
    for (const std::vector<int>& cycle : cycles) {
        for (size_t i = 0; i < cycle.size(); ++i) {
            const int from = cycle[i];
            const int to = cycle[(i + 1) % cycle.size()];
            if (from == 0) {
                throw std::invalid_argument("Permutation labels must be nonzero");
            }
            if (from == to) {
                continue; // singleton cycle: fixed point, not stored
            }
            if (!p.map.emplace(from, to).second) {
                throw std::invalid_argument("Label appears more than once in cycles");
            }
        }
    }
    return p;
}

std::vector<std::vector<int>> ParseCycleNotation(const std::string& input)
{
    // Optionally strip one pair of square brackets wrapping the whole
    // expression (SnapPy PD code style): recognized when the bracket matching
    // the leading '[' is the final character and more brackets occur inside.
    std::string s = input;
    {
        const size_t first = s.find_first_not_of(" \t\r\n");
        const size_t last = s.find_last_not_of(" \t\r\n");
        if (first != std::string::npos && s[first] == '[') {
            int depth = 0;
            size_t match = std::string::npos;
            for (size_t i = first; i <= last; ++i) {
                if (s[i] == '[' || s[i] == '(') {
                    ++depth;
                }
                else if (s[i] == ']' || s[i] == ')') {
                    --depth;
                    if (depth == 0) {
                        match = i;
                        break;
                    }
                }
            }
            const std::string inner =
                (match == last) ? s.substr(first + 1, last - first - 1) : std::string();
            if (match == last &&
                (inner.find('(') != std::string::npos || inner.find('[') != std::string::npos)) {
                s = inner;
            }
        }
    }

    std::vector<std::vector<int>> cycles;
    std::vector<int> current;
    char expectedCloser = '\0'; // nonzero while inside a cycle

    size_t i = 0;
    while (i < s.size()) {
        const char c = s[i];
        if (std::isspace(static_cast<unsigned char>(c)) || c == ',') {
            ++i;
        }
        else if (c == '(' || c == '[') {
            if (expectedCloser != '\0') {
                throw std::invalid_argument("Nested bracket in cycle notation");
            }
            expectedCloser = (c == '(') ? ')' : ']';
            current.clear();
            ++i;
        }
        else if (c == ')' || c == ']') {
            if (expectedCloser == '\0') {
                throw std::invalid_argument("Unmatched closing bracket in cycle notation");
            }
            if (c != expectedCloser) {
                throw std::invalid_argument("Mismatched brackets in cycle notation");
            }
            cycles.push_back(current);
            expectedCloser = '\0';
            ++i;
        }
        else if (c == '-' || std::isdigit(static_cast<unsigned char>(c))) {
            if (expectedCloser == '\0') {
                throw std::invalid_argument("Number outside of a cycle in cycle notation");
            }
            size_t consumed = 0;
            const int value = std::stoi(s.substr(i), &consumed);
            current.push_back(value);
            i += consumed;
        }
        else {
            throw std::invalid_argument(std::string("Unexpected character '") + c + "' in cycle notation");
        }
    }
    if (expectedCloser != '\0') {
        throw std::invalid_argument("Missing closing bracket in cycle notation");
    }

    return cycles;
}

Permutation Permutation::FromString(const std::string& input)
{
    std::vector<std::vector<int>> cycles = ParseCycleNotation(input);

    // PD-style duplicate labels: each edge label should occur once as +k and
    // once as -k; when a label occurs twice with the same sign (e.g. SnapPy PD
    // codes, which are all positive), the second occurrence is negated.
    std::set<int> used;
    for (std::vector<int>& cycle : cycles) {
        if (cycle.empty()) {
            continue;
        }
        for (int& value : cycle) {
            if (used.count(value) != 0 && used.count(-value) == 0) {
                value = -value;
            }
            used.insert(value);
        }
    }
    return FromCycles(cycles);
}

int Permutation::Apply(int h) const
{
    const auto it = map.find(h);
    return (it != map.end()) ? it->second : h;
}

int Permutation::MaxAbsLabel() const
{
    int maxAbs = 0;
    for (const auto& entry : map) {
        maxAbs = std::max(maxAbs, std::abs(entry.first));
    }
    return maxAbs;
}

Permutation Permutation::Inverse() const
{
    Permutation inv;
    for (const auto& entry : map) {
        inv.map.emplace(entry.second, entry.first);
    }
    return inv;
}

Permutation Permutation::operator*(const Permutation& rhs) const
{
    Permutation result;
    for (const auto& entry : map) {
        const int image = rhs.Apply(entry.second);
        if (image != entry.first) {
            result.map.emplace(entry.first, image);
        }
    }
    for (const auto& entry : rhs.map) {
        if (map.count(entry.first) == 0) {
            result.map.emplace(entry.first, entry.second);
        }
    }
    return result;
}

bool Permutation::operator==(const Permutation& rhs) const
{
    return map == rhs.map;
}

bool Permutation::operator!=(const Permutation& rhs) const
{
    return !(*this == rhs);
}

std::vector<std::vector<int>> Permutation::Cycles() const
{
    std::vector<int> labels;
    labels.reserve(map.size());
    for (const auto& entry : map) {
        labels.push_back(entry.first);
    }
    std::sort(labels.begin(), labels.end(), LabelLess);

    std::vector<std::vector<int>> cycles;
    std::set<int> visited;
    for (const int start : labels) {
        if (visited.count(start) != 0) {
            continue;
        }
        std::vector<int> cycle;
        int h = start;
        do {
            cycle.push_back(h);
            visited.insert(h);
            h = Apply(h);
        } while (h != start);
        cycles.push_back(std::move(cycle));
    }
    return cycles;
}

std::string Permutation::ToString() const
{
    const std::vector<std::vector<int>> cycles = Cycles();
    if (cycles.empty()) {
        return "()";
    }
    std::string out;
    for (const std::vector<int>& cycle : cycles) {
        out += '(';
        for (size_t i = 0; i < cycle.size(); ++i) {
            if (i > 0) {
                out += ',';
            }
            out += std::to_string(cycle[i]);
        }
        out += ')';
    }
    return out;
}

Permutation MakeEdgeInvolution(int n)
{
    std::vector<std::vector<int>> cycles;
    cycles.reserve(static_cast<size_t>(n));
    for (int k = 1; k <= n; ++k) {
        cycles.push_back({-k, k});
    }
    return Permutation::FromCycles(cycles);
}
