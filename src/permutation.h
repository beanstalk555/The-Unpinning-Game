#pragma once

#include <string>
#include <unordered_map>
#include <vector>

// A permutation of half-edge labels (nonzero signed integers).
// Acts as the identity on any label it does not explicitly map.
class Permutation
{
public:
    Permutation() = default; // identity

    static Permutation FromCycles(const std::vector<std::vector<int>>& cycles);
    // Parses cycle notation, e.g. "(6,1,-7,-2)(11,2,-12,-3)".
    static Permutation FromString(const std::string& s);

    int Apply(int h) const;
    // Largest absolute value among moved labels (0 for the identity).
    int MaxAbsLabel() const;
    Permutation Inverse() const;
    // Left-to-right composition: (a * b)(x) == b.Apply(a.Apply(x)).
    // With this convention the face permutation is phi = (sigma * epsilon).Inverse().
    Permutation operator*(const Permutation& rhs) const;
    bool operator==(const Permutation& rhs) const;
    bool operator!=(const Permutation& rhs) const;

    // Canonical cycle decomposition: fixed points are omitted, each cycle
    // starts at its smallest label (ordered by absolute value, -k before k),
    // and cycles are sorted by their first label.
    std::vector<std::vector<int>> Cycles() const;
    std::string ToString() const;

private:
    // Invariant: never stores fixed points (h -> h), so two Permutations are
    // equal exactly when their maps are equal.
    std::unordered_map<int, int> map;
};

// The standard edge involution (-1,1)(-2,2)...(-n,n), i.e. h -> -h.
Permutation MakeEdgeInvolution(int n);

// Runs self-tests, printing results to stdout. Returns true if all passed.
bool RunPermutationTests();
