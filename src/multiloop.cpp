#include "multiloop.h"

#include <cstdlib>
#include <functional>
#include <map>
#include <set>

namespace
{
    // Orbits of a permutation acting on all half-edges +-1..+-n, including
    // fixed points. Orbits are listed in order of their smallest label
    // (by absolute value, -k before k), starting at that label.
    std::vector<std::vector<int>> OrbitsOnHalfEdges(const Permutation& p, int n)
    {
        std::vector<std::vector<int>> orbits;
        std::set<int> visited;
        for (int k = 1; k <= n; ++k) {
            for (int sign = -1; sign <= 1; sign += 2) {
                const int start = sign * k;
                if (visited.count(start) != 0) {
                    continue;
                }
                std::vector<int> orbit;
                int h = start;
                do {
                    orbit.push_back(h);
                    visited.insert(h);
                    h = p.Apply(h);
                } while (h != start);
                orbits.push_back(std::move(orbit));
            }
        }
        return orbits;
    }
}

Multiloop::Multiloop(const Permutation& vertexPermutation)
    : numEdges(vertexPermutation.MaxAbsLabel())
    , sigma(vertexPermutation)
    , epsilon(MakeEdgeInvolution(numEdges))
    , phi((sigma * epsilon).Inverse())
    , rho(sigma * sigma * epsilon)
{
    faces = OrbitsOnHalfEdges(phi, numEdges);

    // Strands: union edge labels along each rho orbit. The two orbits of a
    // strand (one per direction) carry the same labels, so they merge.
    std::vector<int> parent(numEdges + 1);
    for (int k = 0; k <= numEdges; ++k) {
        parent[k] = k;
    }
    std::function<int(int)> find = [&](int a) {
        while (parent[a] != a) {
            parent[a] = parent[parent[a]];
            a = parent[a];
        }
        return a;
    };
    for (const std::vector<int>& orbit : OrbitsOnHalfEdges(rho, numEdges)) {
        const int root = find(std::abs(orbit.front()));
        for (const int h : orbit) {
            parent[find(std::abs(h))] = root;
        }
    }

    strandOfEdge.assign(numEdges + 1, -1);
    numStrands = 0;
    std::map<int, int> idOfRoot;
    for (int k = 1; k <= numEdges; ++k) {
        const int root = find(k);
        if (idOfRoot.count(root) == 0) {
            idOfRoot[root] = numStrands++;
        }
        strandOfEdge[k] = idOfRoot[root];
    }
}

std::string ValidateMultiloopSigma(const Permutation& sigma)
{
    const std::vector<std::vector<int>> cycles = sigma.Cycles();
    if (cycles.empty()) {
        return "sigma is empty";
    }
    for (const std::vector<int>& cycle : cycles) {
        if (cycle.size() != 4) {
            return "every crossing must be 4-valent: all sigma cycles must have length 4";
        }
    }

    const int n = sigma.MaxAbsLabel();
    std::set<int> labels;
    for (const std::vector<int>& cycle : cycles) {
        for (const int h : cycle) {
            labels.insert(h);
        }
    }
    for (int k = 1; k <= n; ++k) {
        if (labels.count(k) == 0 || labels.count(-k) == 0) {
            return "half-edge labels must be exactly +-1..+-" + std::to_string(n) +
                   " (label " + std::to_string(k) + " is incomplete)";
        }
    }

    // Connectivity: half-edges at one vertex are joined by sigma, the two
    // halves of an edge by negation.
    std::set<int> reached;
    std::vector<int> stack = {1};
    while (!stack.empty()) {
        const int h = stack.back();
        stack.pop_back();
        if (!reached.insert(h).second) {
            continue;
        }
        stack.push_back(sigma.Apply(h));
        stack.push_back(-h);
    }
    if (static_cast<int>(reached.size()) != 2 * n) {
        return "The multiloop is not connected.";
    }

    const Permutation phi = (sigma * MakeEdgeInvolution(n)).Inverse();
    int numFaces = 0;
    {
        std::set<int> visited;
        for (int k = 1; k <= n; ++k) {
            for (int sign = -1; sign <= 1; sign += 2) {
                int h = sign * k;
                if (visited.count(h) != 0) {
                    continue;
                }
                ++numFaces;
                const int start = h;
                do {
                    visited.insert(h);
                    h = phi.Apply(h);
                } while (h != start);
            }
        }
    }
    const int numVertices = static_cast<int>(cycles.size());
    if (numVertices - n + numFaces != 2) {
        return "sigma does not describe a planar diagram (Euler characteristic is not 2)";
    }
    return "";
}
