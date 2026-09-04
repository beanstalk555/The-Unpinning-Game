#include "multiloop.h"

#include <cstdlib>
#include <functional>
#include <map>
#include <set>
#include <sstream>
#include <unordered_set>

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

    void ApplyPdDuplicateFix(std::vector<std::vector<int>>& cycles)
    {
        std::set<int> used;
        for (std::vector<int>& cycle : cycles) {
            for (int& value : cycle) {
                if (used.count(value) != 0 && used.count(-value) == 0) {
                    value = -value;
                }
                used.insert(value);
            }
        }
    }

    std::vector<int> CollectComponent(const Permutation& sigma, int start, std::set<int>& seen)
    {
        std::vector<int> out;
        std::vector<int> stack = {start};
        while (!stack.empty()) {
            const int h = stack.back();
            stack.pop_back();
            if (!seen.insert(h).second) {
                continue;
            }
            out.push_back(h);
            stack.push_back(sigma.Apply(h));
            stack.push_back(-h);
        }
        return out;
    }

    std::string CycleToString(const std::vector<int>& cycle)
    {
        std::ostringstream out;
        out << "(";
        for (size_t i = 0; i < cycle.size(); ++i) {
            if (i > 0) {
                out << ",";
            }
            out << cycle[i];
        }
        out << ")";
        return out.str();
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

SigmaNotation ParseSigmaNotation(const std::string& s)
{
    std::vector<std::vector<int>> cycles = ParseCycleNotation(s);
    ApplyPdDuplicateFix(cycles);

    std::vector<std::vector<int>> nonempty;
    nonempty.reserve(cycles.size());
    for (const std::vector<int>& cycle : cycles) {
        if (!cycle.empty()) {
            nonempty.push_back(cycle);
        }
    }

    SigmaNotation result;
    result.sigma = Permutation::FromCycles(nonempty);

    std::set<int> emitted;
    for (const std::vector<int>& cycle : cycles) {
        if (cycle.empty()) {
            MultiloopComponent circle;
            circle.kind = MultiloopPart::Circle;
            result.components.push_back(std::move(circle));
            continue;
        }
        const int start = cycle.front();
        if (emitted.count(start) != 0) {
            continue;
        }
        MultiloopComponent part;
        part.kind = MultiloopPart::Crossings;
        part.halfEdges = CollectComponent(result.sigma, start, emitted);
        result.components.push_back(std::move(part));
    }
    return result;
}

std::string FormatSigmaNotation(const Permutation& sigma,
                                const std::vector<MultiloopComponent>& components)
{
    if (components.empty()) {
        return sigma.ToString();
    }
    std::string out;
    const std::vector<std::vector<int>> cycles = sigma.Cycles();
    for (const MultiloopComponent& component : components) {
        if (component.kind == MultiloopPart::Circle) {
            out += "()";
            continue;
        }
        std::unordered_set<int> hs(component.halfEdges.begin(), component.halfEdges.end());
        for (const std::vector<int>& cycle : cycles) {
            if (!cycle.empty() && hs.count(cycle.front()) != 0) {
                out += CycleToString(cycle);
            }
        }
    }
    return out;
}

DenseRelabel RelabelComponent(const Permutation& sigma, const std::vector<int>& halfEdges)
{
    DenseRelabel relabel;
    std::set<int> absLabels;
    for (const int h : halfEdges) {
        if (h != 0) {
            absLabels.insert(std::abs(h));
        }
    }
    int next = 1;
    for (const int a : absLabels) {
        relabel.globalToLocal[a] = next;
        relabel.globalToLocal[-a] = -next;
        relabel.localToGlobal[next] = a;
        relabel.localToGlobal[-next] = -a;
        ++next;
    }

    std::vector<std::vector<int>> localCycles;
    for (const std::vector<int>& cycle : sigma.Cycles()) {
        if (cycle.empty()) {
            continue;
        }
        const auto it = relabel.globalToLocal.find(cycle.front());
        if (it == relabel.globalToLocal.end()) {
            continue;
        }
        std::vector<int> local;
        local.reserve(cycle.size());
        for (const int h : cycle) {
            local.push_back(relabel.globalToLocal.at(h));
        }
        localCycles.push_back(std::move(local));
    }
    relabel.localSigma = Permutation::FromCycles(localCycles);
    return relabel;
}

int FaceIndexContainingHalfEdge(const Multiloop& loop, int halfEdge)
{
    for (int f = 0; f < static_cast<int>(loop.faces.size()); ++f) {
        for (const int h : loop.faces[static_cast<size_t>(f)]) {
            if (h == halfEdge) {
                return f;
            }
        }
    }
    return -1;
}

int FirstFaceOnComponent(const Multiloop& loop, const std::vector<int>& halfEdges)
{
    std::unordered_set<int> hs(halfEdges.begin(), halfEdges.end());
    for (int f = 0; f < static_cast<int>(loop.faces.size()); ++f) {
        for (const int h : loop.faces[static_cast<size_t>(f)]) {
            if (hs.count(h) != 0) {
                return f;
            }
        }
    }
    return -1;
}

void ApplyDefaultExteriors(const Multiloop& loop, std::vector<MultiloopComponent>& components)
{
    for (MultiloopComponent& component : components) {
        if (component.kind != MultiloopPart::Crossings) {
            component.exteriorFace = -1;
            continue;
        }
        if (component.exteriorFace < 0 ||
            component.exteriorFace >= static_cast<int>(loop.faces.size()) ||
            !ComponentOwnsFace(component, loop, component.exteriorFace)) {
            component.exteriorFace = FirstFaceOnComponent(loop, component.halfEdges);
        }
    }
}

std::vector<int> ExteriorFaceIndices(const std::vector<MultiloopComponent>& components)
{
    std::vector<int> faces;
    for (const MultiloopComponent& component : components) {
        if (component.kind == MultiloopPart::Crossings && component.exteriorFace >= 0) {
            faces.push_back(component.exteriorFace);
        }
    }
    return faces;
}

bool ComponentOwnsFace(const MultiloopComponent& component, const Multiloop& loop, int face)
{
    if (component.kind != MultiloopPart::Crossings || face < 0 ||
        face >= static_cast<int>(loop.faces.size())) {
        return false;
    }
    std::unordered_set<int> hs(component.halfEdges.begin(), component.halfEdges.end());
    for (const int h : loop.faces[static_cast<size_t>(face)]) {
        if (hs.count(h) != 0) {
            return true;
        }
    }
    return false;
}

std::string ValidateMultiloopSigma(const Permutation& sigma)
{
    const std::vector<std::vector<int>> cycles = sigma.Cycles();
    if (cycles.empty()) {
        return "";
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

    const Permutation phi = (sigma * MakeEdgeInvolution(n)).Inverse();
    std::set<int> seen;
    for (int k = 1; k <= n; ++k) {
        if (seen.count(k) != 0) {
            continue;
        }
        const std::vector<int> halfEdges = CollectComponent(sigma, k, seen);
        std::unordered_set<int> hset(halfEdges.begin(), halfEdges.end());

        int numVertices = 0;
        for (const std::vector<int>& cycle : cycles) {
            if (!cycle.empty() && hset.count(cycle.front()) != 0) {
                ++numVertices;
            }
        }

        std::set<int> edges;
        for (const int h : halfEdges) {
            edges.insert(std::abs(h));
        }
        const int numEdges = static_cast<int>(edges.size());

        int numFaces = 0;
        std::set<int> visited;
        for (const int start : halfEdges) {
            if (visited.count(start) != 0) {
                continue;
            }
            ++numFaces;
            int h = start;
            do {
                visited.insert(h);
                h = phi.Apply(h);
            } while (h != start);
        }

        if (numVertices - numEdges + numFaces != 2) {
            return "sigma does not describe a planar diagram (Euler characteristic is not 2)";
        }
    }
    return "";
}
