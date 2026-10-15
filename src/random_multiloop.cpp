#include "random_multiloop.h"

#include "move_rewrite.h"
#include "multiloop.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <random>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace
{
    struct P2
    {
        double x;
        double y;
    };

    // For a counterclockwise triangle (a, b, c), is p strictly inside its
    // circumcircle?
    bool InCircumcircle(const P2& a, const P2& b, const P2& c, const P2& p)
    {
        const double ax = a.x - p.x;
        const double ay = a.y - p.y;
        const double bx = b.x - p.x;
        const double by = b.y - p.y;
        const double cx = c.x - p.x;
        const double cy = c.y - p.y;
        const double det = (ax * ax + ay * ay) * (bx * cy - cx * by) -
                           (bx * bx + by * by) * (ax * cy - cx * ay) +
                           (cx * cx + cy * cy) * (ax * by - bx * ay);
        return det > 0.0;
    }

    // Bowyer-Watson incremental Delaunay triangulation. Returns CCW triangles
    // over point indices 0..n-1.
    std::vector<std::array<int, 3>> DelaunayTriangulate(const std::vector<P2>& points)
    {
        const int n = static_cast<int>(points.size());
        std::vector<P2> p = points;
        p.push_back({-1000.0, -1000.0});
        p.push_back({1000.0, -1000.0});
        p.push_back({0.5, 1000.0});

        std::vector<std::array<int, 3>> triangles = {{n, n + 1, n + 2}};
        for (int i = 0; i < n; ++i) {
            std::vector<std::array<int, 3>> bad;
            std::vector<std::array<int, 3>> good;
            for (const std::array<int, 3>& t : triangles) {
                if (InCircumcircle(p[t[0]], p[t[1]], p[t[2]], p[i])) {
                    bad.push_back(t);
                }
                else {
                    good.push_back(t);
                }
            }
            std::set<std::pair<int, int>> directed;
            for (const std::array<int, 3>& t : bad) {
                directed.insert({t[0], t[1]});
                directed.insert({t[1], t[2]});
                directed.insert({t[2], t[0]});
            }
            triangles = std::move(good);
            for (const std::pair<int, int>& e : directed) {
                if (directed.count({e.second, e.first}) == 0) {
                    triangles.push_back({e.first, e.second, i});
                }
            }
        }

        std::vector<std::array<int, 3>> result;
        for (const std::array<int, 3>& t : triangles) {
            if (t[0] < n && t[1] < n && t[2] < n) {
                result.push_back(t);
            }
        }
        return result;
    }

    bool IsConnected(int numVertices, const std::vector<std::pair<int, int>>& edges,
                     const std::vector<char>& alive)
    {
        std::vector<std::vector<int>> adjacency(numVertices);
        for (size_t e = 0; e < edges.size(); ++e) {
            if (alive[e]) {
                adjacency[edges[e].first].push_back(edges[e].second);
                adjacency[edges[e].second].push_back(edges[e].first);
            }
        }
        std::vector<char> seen(numVertices, 0);
        std::vector<int> stack = {0};
        seen[0] = 1;
        int count = 0;
        while (!stack.empty()) {
            const int v = stack.back();
            stack.pop_back();
            ++count;
            for (const int w : adjacency[v]) {
                if (!seen[w]) {
                    seen[w] = 1;
                    stack.push_back(w);
                }
            }
        }
        return count == numVertices;
    }

    // Bigon hook: replaces the 4-cycle at index ci, rotated so it reads
    // (p, q, r, s), by the two 4-cycles (p, q, x, y) and (r, s, -y, -x) where
    // x = n + 1 and y = n + 2 are new edge labels forming a bigon between the
    // two new crossings.
    void BigonHook(std::vector<std::vector<int>>& cycles, size_t ci, int rotation, int numEdges)
    {
        const std::vector<int> c = cycles[ci];
        const int p = c[rotation % 4];
        const int q = c[(rotation + 1) % 4];
        const int r = c[(rotation + 2) % 4];
        const int s = c[(rotation + 3) % 4];
        const int x = numEdges + 1;
        const int y = numEdges + 2;
        cycles[ci] = {p, q, x, y};
        cycles.push_back({r, s, -y, -x});
    }

    // Repeatedly hooks crossings where two distinct strands meet, joining
    // them, until the multiloop is a single strand.
    void ForceOneStrand(Permutation& sigma, std::mt19937& rng)
    {
        while (true) {
            const Multiloop loop(sigma);
            if (loop.numStrands <= 1) {
                return;
            }
            std::vector<std::vector<int>> cycles = sigma.Cycles();
            std::vector<size_t> candidates;
            for (size_t i = 0; i < cycles.size(); ++i) {
                if (loop.strandOfEdge[std::abs(cycles[i][0])] !=
                    loop.strandOfEdge[std::abs(cycles[i][1])]) {
                    candidates.push_back(i);
                }
            }
            // Connectivity guarantees two strands meet somewhere.
            const size_t pick = candidates[rng() % candidates.size()];
            // (a,b,c,d) -> (x,y,b,c)(a,-y,-x,d): split between (b,c) and (d,a).
            BigonHook(cycles, pick, 1, loop.numEdges);
            sigma = Permutation::FromCycles(cycles);
        }
    }

    // Repeatedly hooks crossings where a strand meets itself until every
    // strand is a simple closed curve. Each hook removes the chosen
    // self-intersection (splitting its strand in two) and creates none.
    void ForceSimple(Permutation& sigma, std::mt19937& rng)
    {
        while (true) {
            const Multiloop loop(sigma);
            std::vector<std::vector<int>> cycles = sigma.Cycles();
            std::vector<size_t> candidates;
            for (size_t i = 0; i < cycles.size(); ++i) {
                if (loop.strandOfEdge[std::abs(cycles[i][0])] ==
                    loop.strandOfEdge[std::abs(cycles[i][1])]) {
                    candidates.push_back(i);
                }
            }
            if (candidates.empty()) {
                return;
            }
            const size_t pick = candidates[rng() % candidates.size()];
            const std::vector<int>& c = cycles[pick];
            // Follow the strand out along d's edge until it re-enters this
            // crossing: at a self-intersection (a,b,c,d) the first return is
            // at a or c, determining the pairing of consecutive labels.
            int h = c[3];
            int arrival = 0;
            while (true) {
                const int arr = -h; // half-edge where the travelled edge ends
                if (arr == c[0] || arr == c[1] || arr == c[2] || arr == c[3]) {
                    arrival = arr;
                    break;
                }
                h = sigma.Apply(sigma.Apply(arr)); // straight ahead
            }
            // d returns to a: hook (a,b,c,d) -> (a,b,x,y)(c,d,-y,-x).
            // d returns to c: the same with the cycle read as (b,c,d,a).
            const int rotation = (arrival == c[0]) ? 0 : 1;
            BigonHook(cycles, pick, rotation, loop.numEdges);
            sigma = Permutation::FromCycles(cycles);
        }
    }

    int RemapExteriorHalfEdge(const Multiloop& loop, const RewriteResult& rw, int halfEdge)
    {
        const int direct = rw.MapHalfEdge(halfEdge);
        if (direct != 0) {
            return direct;
        }
        for (const std::vector<int>& face : loop.faces) {
            bool onFace = false;
            for (const int h : face) {
                if (h == halfEdge) {
                    onFace = true;
                    break;
                }
            }
            if (!onFace) {
                continue;
            }
            for (const int h : face) {
                const int mapped = rw.MapHalfEdge(h);
                if (mapped != 0) {
                    return mapped;
                }
            }
        }
        return 0;
    }

    bool HalfEdgeOnFace(const Multiloop& loop, int halfEdge)
    {
        for (const std::vector<int>& face : loop.faces) {
            for (const int h : face) {
                if (h == halfEdge) {
                    return true;
                }
            }
        }
        return false;
    }

    // One adjacent smoothing of cycle (a, b, c, d). rotation 0 pairs (a, b)
    // with (c, d); rotation 1 pairs (b, c) with (d, a). A pair whose two
    // labels are the ends of one edge is deleted, and the other pair is
    // spliced, matching an R1.
    bool TryAdjacentSmoothing(const Multiloop& loop, int crossing, int rotation, RewriteResult* out)
    {
        const std::vector<std::vector<int>> cycles = loop.sigma.Cycles();
        if (crossing < 0 || crossing >= static_cast<int>(cycles.size()) ||
            cycles[static_cast<size_t>(crossing)].size() != 4) {
            return false;
        }
        const std::vector<int>& cycle = cycles[static_cast<size_t>(crossing)];
        const int h0 = cycle[static_cast<size_t>(rotation % 4)];
        const int h1 = cycle[static_cast<size_t>((rotation + 1) % 4)];
        const int h2 = cycle[static_cast<size_t>((rotation + 2) % 4)];
        const int h3 = cycle[static_cast<size_t>((rotation + 3) % 4)];
        const int pairs[2][2] = {{h0, h1}, {h2, h3}};
        std::vector<int> deletedEdges;
        std::vector<SplicePair> splices;
        for (int p = 0; p < 2; ++p) {
            const int u = pairs[p][0];
            const int v = pairs[p][1];
            if (std::abs(u) == std::abs(v)) {
                deletedEdges.push_back(std::abs(u));
            }
            else {
                splices.push_back(SplicePair{u, v});
            }
        }
        RewriteResult rw;
        if (!DeleteCrossingsAndSplice(loop, {crossing}, deletedEdges, splices, &rw)) {
            return false;
        }
        if (!rw.circleArcs.empty() || !ValidateMultiloopSigma(rw.sigma).empty()) {
            return false;
        }
        try {
            const Multiloop next(rw.sigma);
            const int crossings = static_cast<int>(next.sigma.Cycles().size());
            if (next.numStrands != 1 || crossings != static_cast<int>(cycles.size()) - 1 ||
                crossings <= 0) {
                return false;
            }
        }
        catch (const std::exception&) {
            return false;
        }
        *out = std::move(rw);
        return true;
    }

    // Deletes one uniformly chosen crossing by the adjacent smoothing that
    // keeps a single strand. Returns false when that smoothing is not unique
    // or the exterior face cannot be followed.
    bool ResolveOneCrossing(RandomMultiloop* sample, std::mt19937& rng)
    {
        const Multiloop loop(sample->sigma);
        const std::vector<std::vector<int>> cycles = loop.sigma.Cycles();
        if (cycles.empty()) {
            return false;
        }
        const int crossing = static_cast<int>(rng() % cycles.size());
        RewriteResult kept;
        int winners = 0;
        for (int rotation = 0; rotation < 2; ++rotation) {
            RewriteResult rw;
            if (!TryAdjacentSmoothing(loop, crossing, rotation, &rw)) {
                continue;
            }
            kept = std::move(rw);
            ++winners;
        }
        if (winners != 1) {
            return false;
        }
        const int exterior = RemapExteriorHalfEdge(loop, kept, sample->exteriorHalfEdge);
        if (exterior == 0) {
            return false;
        }
        const Multiloop next(kept.sigma);
        if (!HalfEdgeOnFace(next, exterior)) {
            return false;
        }
        sample->sigma = kept.sigma;
        sample->exteriorHalfEdge = exterior;
        return true;
    }
}

RandomMultiloop GenerateRandomMultiloop(int minFaces, double density, bool forceOneStrand,
                                        bool forceSimple, unsigned seed)
{
    const int numPoints = std::max(3, minFaces);
    density = std::max(0.0, std::min(1.0, density));
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    for (int attempt = 0; attempt < 5; ++attempt) {
        // 1. Random point cloud in the unit square.
        std::vector<P2> points(numPoints);
        for (P2& p : points) {
            p = {uniform(rng), uniform(rng)};
        }

        // 2. Delaunay triangulation, reduced to its edge set.
        const std::vector<std::array<int, 3>> triangles = DelaunayTriangulate(points);
        std::set<std::pair<int, int>> edgeSet;
        for (const std::array<int, 3>& t : triangles) {
            for (int i = 0; i < 3; ++i) {
                const int a = t[i];
                const int b = t[(i + 1) % 3];
                edgeSet.insert({std::min(a, b), std::max(a, b)});
            }
        }
        std::vector<std::pair<int, int>> edges(edgeSet.begin(), edgeSet.end());
        const int fullEdges = static_cast<int>(edges.size());
        if (fullEdges < numPoints - 1 || !IsConnected(numPoints, edges, std::vector<char>(fullEdges, 1))) {
            continue; // degenerate triangulation; retry with new points
        }

        // 3. Random edge deletion down to the density target, keeping the
        // graph connected. Bridges stay bridges as edges disappear, so one
        // shuffled pass suffices.
        const int treeEdges = numPoints - 1;
        const int target = treeEdges +
            static_cast<int>(std::lround(density * static_cast<double>(fullEdges - treeEdges)));
        std::vector<int> order(fullEdges);
        for (int i = 0; i < fullEdges; ++i) {
            order[i] = i;
        }
        std::shuffle(order.begin(), order.end(), rng);
        std::vector<char> alive(fullEdges, 1);
        int aliveCount = fullEdges;
        for (const int e : order) {
            if (aliveCount <= target) {
                break;
            }
            alive[e] = 0;
            if (IsConnected(numPoints, edges, alive)) {
                --aliveCount;
            }
            else {
                alive[e] = 1;
            }
        }

        std::vector<std::pair<int, int>> kept;
        for (int e = 0; e < fullEdges; ++e) {
            if (alive[e]) {
                kept.push_back(edges[e]);
            }
        }
        const int numEdges = static_cast<int>(kept.size());

        // 4. Combinatorial embedding: darts 2e (u -> v) and 2e + 1 (v -> u),
        // with the counterclockwise rotation at each vertex from the geometry.
        const auto tailOf = [&](int d) { return (d & 1) ? kept[d >> 1].second : kept[d >> 1].first; };
        const auto headOf = [&](int d) { return (d & 1) ? kept[d >> 1].first : kept[d >> 1].second; };
        std::vector<std::vector<int>> dartsAt(numPoints);
        for (int d = 0; d < 2 * numEdges; ++d) {
            dartsAt[tailOf(d)].push_back(d);
        }
        std::vector<int> rotNext(2 * numEdges);
        for (int v = 0; v < numPoints; ++v) {
            std::vector<int>& darts = dartsAt[v];
            std::sort(darts.begin(), darts.end(), [&](int a, int b) {
                const double angA = std::atan2(points[headOf(a)].y - points[v].y,
                                               points[headOf(a)].x - points[v].x);
                const double angB = std::atan2(points[headOf(b)].y - points[v].y,
                                               points[headOf(b)].x - points[v].x);
                return angA < angB;
            });
            for (size_t i = 0; i < darts.size(); ++i) {
                rotNext[darts[i]] = darts[(i + 1) % darts.size()];
            }
        }

        // The dart preceding d in its face walk (faces kept on the left).
        const auto prevInFace = [&](int d) { return rotNext[d] ^ 1; };

        // 5. Medial graph as a multiloop: one 4-valent vertex per kept edge,
        // one medial edge per dart d, labeled d + 1, with the +(d + 1) end on
        // the medial vertex of edge(d). Counterclockwise around the medial
        // vertex of e with dart d = 2e:
        //   (+(d+1), -(prev(d)+1), +(rev(d)+1), -(prev(rev(d))+1)).
        std::vector<std::vector<int>> cycles;
        for (int e = 0; e < numEdges; ++e) {
            const int d = 2 * e;
            cycles.push_back({d + 1, -(prevInFace(d) + 1), (d ^ 1) + 1, -(prevInFace(d ^ 1) + 1)});
        }

        RandomMultiloop result;
        result.sigma = Permutation::FromCycles(cycles);

        // 6. The unbounded region: at the bottommost vertex u, the outer face
        // lies between the most-clockwise and most-counterclockwise darts;
        // its walk contains the reverse of the smallest-angle dart, and the
        // corresponding multiloop face contains that dart's +(d + 1) label.
        int bottom = 0;
        for (int v = 1; v < numPoints; ++v) {
            if (points[v].y < points[bottom].y ||
                (points[v].y == points[bottom].y && points[v].x < points[bottom].x)) {
                bottom = v;
            }
        }
        int minDart = -1;
        double minAngle = 0.0;
        for (const int d : dartsAt[bottom]) {
            double ang = std::atan2(points[headOf(d)].y - points[bottom].y,
                                    points[headOf(d)].x - points[bottom].x);
            if (ang < 0.0) {
                ang += 6.283185307179586;
            }
            if (minDart == -1 || ang < minAngle) {
                minDart = d;
                minAngle = ang;
            }
        }
        result.exteriorHalfEdge = (minDart ^ 1) + 1;

        // 7. Optional bigon-hook post-processing. Hooks never remove labels,
        // so the exterior half-edge stays on the unbounded face.
        if (forceOneStrand) {
            ForceOneStrand(result.sigma, rng);
        }
        else if (forceSimple) {
            ForceSimple(result.sigma, rng);
        }
        return result;
    }
    throw std::runtime_error("Random multiloop generation failed (degenerate point cloud)");
}

RandomMultiloop MakeOddOneStrandMultiloop(unsigned seed)
{
    std::mt19937 rng(seed);
    // Same request as the sandbox default: 10 points, mesh density 0.8.
    constexpr int kPoints = 10;
    constexpr double kDensity = 0.8;
    constexpr int kAttempts = 48;
    for (int attempt = 0; attempt < kAttempts; ++attempt) {
        RandomMultiloop sample;
        try {
            sample = GenerateRandomMultiloop(kPoints, kDensity, true, false, rng());
        }
        catch (const std::exception&) {
            continue;
        }
        {
            const Multiloop loop(sample.sigma);
            const int crossings = static_cast<int>(loop.sigma.Cycles().size());
            if (loop.numStrands != 1 || crossings <= 0) {
                continue;
            }
            if (crossings % 2 == 0 && !ResolveOneCrossing(&sample, rng)) {
                continue;
            }
        }
        const Multiloop played(sample.sigma);
        const int crossings = static_cast<int>(played.sigma.Cycles().size());
        if (played.numStrands != 1 || crossings <= 0 || crossings % 2 == 0 ||
            !HalfEdgeOnFace(played, sample.exteriorHalfEdge)) {
            continue;
        }
        return sample;
    }
    throw std::runtime_error("Random multiloop generation failed (no odd one-strand diagram)");
}
