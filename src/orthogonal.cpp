// C++ port of SnapPy / spherogram's orthogonal.py (https://github.com/3-manifolds/Spherogram/blob/3062478dd69a52a16e55b967dd46dfe940af9ea7/spherogram_src/links/orthogonal.py),
// adapted to the permutation representation of a multiloop. 
//  
// References:
//
// Tamassia, On embedding a graph in the grid with the minimum number of
// bends. Siam J. Comput. 16 (1987) http://dx.doi.org/10.1137/0216030
//
// and
//
// Bridgeman et. al., Turn-Regularity and Planar Orthogonal Drawings
// ftp://ftp.cs.brown.edu/pub/techreports/99/cs99-04.pdf
// A more concise summary of the algorithm is contained in
//
// Hashemi and Tahmasbi, A better heuristic for area-compaction of orthogonal
// representations.  https://dx.doi.org/10.1016/j.amc.2005.03.007
//
// All crossings are 4-valent, so the flow network has no vertex arcs and
// every crossing corner is a 90-degree angle.

#include "orthogonal.h"

#include "geometry.h"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <deque>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

namespace
{
    constexpr int kInf = std::numeric_limits<int>::max() / 4;

    void Require(bool condition, const char* message)
    {
        if (!condition) {
            throw std::runtime_error(std::string("Grid layout failed: ") + message);
        }
    }

    // ----------------------------------------------------------------
    // Minimum cost maximum flow (successive shortest augmenting paths).
    // ----------------------------------------------------------------

    class MinCostFlow
    {
    public:
        explicit MinCostFlow(int numNodes) : graph(numNodes) {}

        int AddEdge(int from, int to, int capacity, int cost)
        {
            graph[from].push_back(static_cast<int>(edges.size()));
            edges.push_back({to, capacity, cost, 0});
            graph[to].push_back(static_cast<int>(edges.size()));
            edges.push_back({from, 0, -cost, 0});
            return static_cast<int>(edges.size()) - 2;
        }

        int FlowOn(int edgeId) const { return edges[edgeId].flow; }

        int Solve(int source, int sink)
        {
            const int n = static_cast<int>(graph.size());
            int total = 0;
            for (;;) {
                std::vector<int> dist(n, kInf);
                std::vector<int> prevEdge(n, -1);
                std::vector<char> inQueue(n, 0);
                std::deque<int> queue;
                dist[source] = 0;
                queue.push_back(source);
                while (!queue.empty()) {
                    const int u = queue.front();
                    queue.pop_front();
                    inQueue[u] = 0;
                    for (const int id : graph[u]) {
                        const Edge& e = edges[id];
                        if (e.capacity - e.flow > 0 && dist[u] + e.cost < dist[e.to]) {
                            dist[e.to] = dist[u] + e.cost;
                            prevEdge[e.to] = id;
                            if (!inQueue[e.to]) {
                                inQueue[e.to] = 1;
                                queue.push_back(e.to);
                            }
                        }
                    }
                }
                if (prevEdge[sink] == -1) {
                    break;
                }
                int augment = kInf;
                for (int v = sink; v != source;) {
                    const Edge& e = edges[prevEdge[v]];
                    augment = std::min(augment, e.capacity - e.flow);
                    v = edges[prevEdge[v] ^ 1].to;
                }
                for (int v = sink; v != source;) {
                    edges[prevEdge[v]].flow += augment;
                    edges[prevEdge[v] ^ 1].flow -= augment;
                    v = edges[prevEdge[v] ^ 1].to;
                }
                total += augment;
            }
            return total;
        }

    private:
        struct Edge
        {
            int to;
            int capacity;
            int cost;
            int flow;
        };
        std::vector<Edge> edges;
        std::vector<std::vector<int>> graph;
    };

    // ----------------------------------------------------------------
    // Directions and compass slots.
    // ----------------------------------------------------------------

    // Matches the reference implementation: dirs = [left, up, right, down];
    // a turn of t advances the index by t (mod 4).
    enum Dir { kLeft = 0, kUp = 1, kRight = 2, kDown = 3 };

    // Compass slots around a vertex in counterclockwise link order, matching
    // the reference sort key (edge.head == v, edge is vertical).
    enum Slot { kWest = 0, kSouth = 1, kEast = 2, kNorth = 3 };

    // ----------------------------------------------------------------
    // Orthogonal representation: a planar graph whose edges are each
    // horizontal (oriented rightward) or vertical (oriented upward).
    // ----------------------------------------------------------------

    struct OEdge
    {
        int tail;
        int head;
        bool horizontal;
        bool dummy;
    };

    struct OFaceElem
    {
        int edge;
        int vertex; // the endpoint giving the clockwise orientation
    };

    struct OFace
    {
        std::vector<OFaceElem> elems;
        std::vector<int> turns; // turns[i] is the turn at elems[i]'s vertex
        bool exterior = false;
    };

    struct SwitchInfo
    {
        int index;
        bool isSink;
        int turn;
    };

    // partial_sums + kitty_corner from the reference code.
    bool FindKittyCorner(const std::vector<int>& turns, int* r0, int* r1)
    {
        std::vector<int> rotations(turns.size() + 1, 0);
        for (size_t i = 0; i < turns.size(); ++i) {
            rotations[i + 1] = rotations[i] + turns[i];
        }
        std::vector<int> reflex;
        for (size_t i = 0; i < turns.size(); ++i) {
            if (turns[i] == -1) {
                reflex.push_back(static_cast<int>(i));
            }
        }
        for (size_t a = 0; a < reflex.size(); ++a) {
            for (size_t b = a + 1; b < reflex.size(); ++b) {
                if (rotations[reflex[b]] - rotations[reflex[a]] == 2) {
                    *r0 = reflex[a];
                    *r1 = reflex[b];
                    return true;
                }
            }
        }
        return false;
    }

    std::vector<std::pair<int, int>> SaturateFace(std::vector<SwitchInfo> info)
    {
        for (size_t i = 0; i < info.size(); ++i) {
            if (info[i].turn == -1) {
                std::rotate(info.begin(), info.begin() + i, info.end());
                break;
            }
        }
        for (size_t i = 0; i + 2 < info.size(); ++i) {
            const SwitchInfo x = info[i];
            const SwitchInfo y = info[i + 1];
            const SwitchInfo z = info[i + 2];
            if (x.turn == -1 && y.turn == 1 && z.turn == 1) {
                SwitchInfo a = x;
                SwitchInfo b = z;
                if (!x.isSink) {
                    std::swap(a, b);
                }
                std::vector<SwitchInfo> remaining(info.begin(), info.begin() + i);
                remaining.push_back({z.index, z.isSink, 1});
                remaining.insert(remaining.end(), info.begin() + i + 3, info.end());
                std::vector<std::pair<int, int>> rest = SaturateFace(std::move(remaining));
                rest.insert(rest.begin(), {a.index, b.index});
                return rest;
            }
        }
        return {};
    }

    struct DagEdge
    {
        int tail;
        int head;
        bool dummy;
    };

    // topological_numbering from the reference code: a layered numbering of a
    // DAG followed by local moves that shorten non-dummy edges.
    std::vector<int> TopologicalNumbering(int numNodes, const std::vector<DagEdge>& dagEdges)
    {
        std::vector<std::vector<int>> out(numNodes);
        std::vector<std::vector<int>> in(numNodes);
        std::vector<int> inDegree(numNodes, 0);
        for (size_t i = 0; i < dagEdges.size(); ++i) {
            out[dagEdges[i].tail].push_back(static_cast<int>(i));
            in[dagEdges[i].head].push_back(static_cast<int>(i));
            ++inDegree[dagEdges[i].head];
        }

        std::vector<int> numbering(numNodes, 0);
        std::vector<int> current;
        for (int v = 0; v < numNodes; ++v) {
            if (inDegree[v] == 0) {
                current.push_back(v);
            }
        }
        int remaining = numNodes;
        int layer = 0;
        while (remaining > 0) {
            Require(!current.empty(), "constraint graph has a cycle");
            std::vector<int> next;
            for (const int v : current) {
                numbering[v] = layer;
                --remaining;
                for (const int id : out[v]) {
                    if (--inDegree[dagEdges[id].head] == 0) {
                        next.push_back(dagEdges[id].head);
                    }
                }
            }
            current = std::move(next);
            ++layer;
        }

        for (int pass = 0; pass < 1000; ++pass) {
            bool success = false;
            for (int v = 0; v < numNodes; ++v) {
                int below = 0;
                for (const int id : in[v]) {
                    if (!dagEdges[id].dummy) {
                        ++below;
                    }
                }
                int above = 0;
                for (const int id : out[v]) {
                    if (!dagEdges[id].dummy) {
                        ++above;
                    }
                }
                if (above == below) {
                    continue;
                }
                int newPos;
                if (above > below) {
                    newPos = kInf;
                    for (const int id : out[v]) {
                        newPos = std::min(newPos, numbering[dagEdges[id].head]);
                    }
                    --newPos;
                }
                else {
                    newPos = -kInf;
                    for (const int id : in[v]) {
                        newPos = std::max(newPos, numbering[dagEdges[id].tail]);
                    }
                    ++newPos;
                }
                if (newPos != numbering[v]) {
                    numbering[v] = newPos;
                    success = true;
                }
            }
            if (!success) {
                break;
            }
        }
        return numbering;
    }

    class OrthoRep
    {
    public:
        explicit OrthoRep(int numVertices)
            : slots(numVertices, std::array<int, 4>{{-1, -1, -1, -1}})
        {
        }

        int AddEdge(int tail, int head, bool horizontal, bool dummy)
        {
            const int id = static_cast<int>(edges.size());
            const int tailSlot = horizontal ? kEast : kNorth;
            const int headSlot = horizontal ? kWest : kSouth;
            Require(slots[tail][tailSlot] == -1 && slots[head][headSlot] == -1,
                    "two edges leave a vertex in the same direction");
            slots[tail][tailSlot] = id;
            slots[head][headSlot] = id;
            edges.push_back({tail, head, horizontal, dummy});
            return id;
        }

        int OtherEnd(int edgeId, int v) const
        {
            const OEdge& e = edges[edgeId];
            return e.tail == v ? e.head : e.tail;
        }

        // The next edge counterclockwise around v.
        int NextEdgeAtVertex(int edgeId, int v) const
        {
            int current = -1;
            for (int s = 0; s < 4; ++s) {
                if (slots[v][s] == edgeId) {
                    current = s;
                }
            }
            Require(current != -1, "edge missing from vertex link");
            for (int step = 1; step <= 4; ++step) {
                const int s = (current + step) % 4;
                if (slots[v][s] != -1) {
                    return slots[v][s];
                }
            }
            Require(false, "isolated edge in vertex link");
            return -1;
        }

        OFace TraceFace(int edgeId, int v) const
        {
            OFace face;
            face.elems.push_back({edgeId, v});
            for (;;) {
                const OFaceElem& last = face.elems.back();
                const int e = NextEdgeAtVertex(last.edge, last.vertex);
                const int w = OtherEnd(e, last.vertex);
                if (e == face.elems.front().edge && w == face.elems.front().vertex) {
                    break;
                }
                face.elems.push_back({e, w});
            }

            const size_t n = face.elems.size();
            int rotation = 0;
            for (size_t i = 0; i < n; ++i) {
                const OEdge& e0 = edges[face.elems[i].edge];
                const OEdge& e1 = edges[face.elems[(i + 1) % n].edge];
                const int v0 = face.elems[i].vertex;
                if (e0.horizontal == e1.horizontal) {
                    face.turns.push_back(0);
                }
                else {
                    const bool t = ((e0.tail == v0) ^ (e1.head == v0) ^ e0.horizontal);
                    face.turns.push_back(t ? -1 : 1);
                }
                rotation += face.turns.back();
            }
            Require(std::abs(rotation) == 4, "face rotation is not +-4");
            face.exterior = (rotation == -4);
            return face;
        }

        void BuildFaces()
        {
            faces.clear();
            std::set<std::pair<int, int>> sides;
            for (size_t e = 0; e < edges.size(); ++e) {
                sides.insert({static_cast<int>(e), edges[e].tail});
                sides.insert({static_cast<int>(e), edges[e].head});
            }
            while (!sides.empty()) {
                const std::pair<int, int> seed = *sides.begin();
                OFace face = TraceFace(seed.first, seed.second);
                for (const OFaceElem& elem : face.elems) {
                    sides.erase({elem.edge, elem.vertex});
                }
                faces.push_back(std::move(face));
            }
        }

        // Adds dummy edges between kitty corners until all interior faces
        // are turn-regular.
        void MakeTurnRegular()
        {
            BuildFaces();
            std::vector<OFace> regular;
            std::vector<OFace> irregular;
            int r0 = 0;
            int r1 = 0;
            for (OFace& face : faces) {
                if (face.exterior || !FindKittyCorner(face.turns, &r0, &r1)) {
                    regular.push_back(std::move(face));
                }
                else {
                    irregular.push_back(std::move(face));
                }
            }
            while (!irregular.empty()) {
                OFace face = std::move(irregular.back());
                irregular.pop_back();
                Require(FindKittyCorner(face.turns, &r0, &r1), "lost kitty corner");
                const int v0 = face.elems[r0].vertex;
                const int v1 = face.elems[r1].vertex;
                Require(v0 != v1, "degenerate kitty corner");

                // Both free slots of a reflex corner open into this face, and
                // kitty corners face opposite directions, so some axis pairing
                // of free slots exists.
                int e;
                if (slots[v0][kEast] == -1 && slots[v1][kWest] == -1) {
                    e = AddEdge(v0, v1, true, true);
                }
                else if (slots[v0][kNorth] == -1 && slots[v1][kSouth] == -1) {
                    e = AddEdge(v0, v1, false, true);
                }
                else if (slots[v0][kWest] == -1 && slots[v1][kEast] == -1) {
                    e = AddEdge(v1, v0, true, true);
                }
                else if (slots[v0][kSouth] == -1 && slots[v1][kNorth] == -1) {
                    e = AddEdge(v1, v0, false, true);
                }
                else {
                    Require(false, "cannot insert turn-regularity edge");
                    return;
                }

                for (const int v : {v0, v1}) {
                    OFace newFace = TraceFace(e, v);
                    if (newFace.exterior || !FindKittyCorner(newFace.turns, &r0, &r1)) {
                        regular.push_back(std::move(newFace));
                    }
                    else {
                        irregular.push_back(std::move(newFace));
                    }
                }
            }
            faces = std::move(regular);
        }

        std::vector<std::pair<int, int>> SaturationEdges(bool swapHorizontal) const
        {
            std::vector<std::pair<int, int>> result;
            for (const OFace& face : faces) {
                const size_t n = face.elems.size();
                std::vector<SwitchInfo> info;
                for (size_t i = 0; i < n; ++i) {
                    const OEdge& e0 = edges[face.elems[i].edge];
                    const OEdge& e1 = edges[face.elems[(i + 1) % n].edge];
                    const int v0 = face.elems[i].vertex;
                    const bool swap0 = swapHorizontal && e0.horizontal;
                    const bool swap1 = swapHorizontal && e1.horizontal;
                    const int tail0 = swap0 ? e0.head : e0.tail;
                    const int head0 = swap0 ? e0.tail : e0.head;
                    const int tail1 = swap1 ? e1.head : e1.tail;
                    const int head1 = swap1 ? e1.tail : e1.head;
                    if (tail0 == v0 && tail1 == v0) {
                        info.push_back({static_cast<int>(i), false, face.turns[i]});
                    }
                    else if (head0 == v0 && head1 == v0) {
                        info.push_back({static_cast<int>(i), true, face.turns[i]});
                    }
                }
                for (const std::pair<int, int>& p : SaturateFace(std::move(info))) {
                    result.push_back({face.elems[p.first].vertex, face.elems[p.second].vertex});
                }
            }
            return result;
        }

        // Coordinates from chains of edges of one kind: chains of horizontal
        // edges give y, chains of vertical edges give x.
        std::vector<int> ChainCoordinates(bool horizontalChains) const
        {
            const int numVertices = static_cast<int>(slots.size());
            std::vector<int> parent(numVertices);
            for (int v = 0; v < numVertices; ++v) {
                parent[v] = v;
            }
            const std::function<int(int)> find = [&](int a) {
                while (parent[a] != a) {
                    parent[a] = parent[parent[a]];
                    a = parent[a];
                }
                return a;
            };
            for (const OEdge& e : edges) {
                if (e.horizontal == horizontalChains) {
                    parent[find(e.tail)] = find(e.head);
                }
            }

            std::vector<int> chainId(numVertices, -1);
            int numChains = 0;
            for (int v = 0; v < numVertices; ++v) {
                const int root = find(v);
                if (chainId[root] == -1) {
                    chainId[root] = numChains++;
                }
            }
            const auto chainOf = [&](int v) { return chainId[find(v)]; };

            std::vector<DagEdge> dagEdges;
            for (const OEdge& e : edges) {
                if (e.horizontal != horizontalChains) {
                    Require(chainOf(e.tail) != chainOf(e.head), "edge within its own chain");
                    dagEdges.push_back({chainOf(e.tail), chainOf(e.head), e.dummy});
                }
            }
            for (const std::pair<int, int>& p : SaturationEdges(false)) {
                if (chainOf(p.first) != chainOf(p.second)) {
                    dagEdges.push_back({chainOf(p.first), chainOf(p.second), true});
                }
            }
            for (const std::pair<int, int>& p : SaturationEdges(true)) {
                int u = p.first;
                int v = p.second;
                if (!horizontalChains) {
                    std::swap(u, v);
                }
                if (chainOf(u) != chainOf(v)) {
                    dagEdges.push_back({chainOf(u), chainOf(v), true});
                }
            }

            const std::vector<int> numbering = TopologicalNumbering(numChains, dagEdges);
            std::vector<int> coords(numVertices);
            for (int v = 0; v < numVertices; ++v) {
                coords[v] = numbering[chainOf(v)];
            }
            return coords;
        }

        std::vector<std::array<int, 4>> slots; // per-vertex W/S/E/N edge ids
        std::vector<OEdge> edges;
        std::vector<OFace> faces;
    };
}

GridLayout ComputeGridLayout(const Multiloop& loop, int exteriorFace)
{
    const int n = loop.numEdges;
    Require(n > 0, "empty multiloop");
    const std::vector<std::vector<int>>& faces = loop.faces;
    const int numFaces = static_cast<int>(faces.size());
    Require(exteriorFace >= 0 && exteriorFace < numFaces, "invalid face index");

    std::unordered_map<int, int> faceOf;
    for (int f = 0; f < numFaces; ++f) {
        for (const int h : faces[f]) {
            faceOf[h] = f;
        }
    }
    const std::vector<std::vector<int>> crossings = loop.sigma.Cycles();
    const int numCrossings = static_cast<int>(crossings.size());
    std::unordered_map<int, int> crossingOf;
    for (int c = 0; c < numCrossings; ++c) {
        for (const int h : crossings[c]) {
            crossingOf[h] = c;
        }
    }
    for (int k = 1; k <= n; ++k) {
        Require(faceOf.count(k) && faceOf.count(-k) && crossingOf.count(k) && crossingOf.count(-k),
                "incomplete half-edge labels");
    }

    // ---- Tamassia's flow network: minimize the number of bends ----

    const int source = numFaces;
    const int sink = numFaces + 1;
    MinCostFlow flow(numFaces + 2);
    int sourceDemand = 0;
    int sinkDemand = 0;
    for (int f = 0; f < numFaces; ++f) {
        const int k = static_cast<int>(faces[f].size());
        const int sourceCap = (f == exteriorFace) ? 0 : std::max(4 - k, 0);
        const int sinkCap = (f == exteriorFace) ? k + 4 : std::max(k - 4, 0);
        if (sourceCap > 0) {
            flow.AddEdge(source, f, sourceCap, 0);
            sourceDemand += sourceCap;
        }
        if (sinkCap > 0) {
            flow.AddEdge(f, sink, sinkCap, 0);
            sinkDemand += sinkCap;
        }
    }
    Require(sourceDemand == sinkDemand, "flow demands do not balance (diagram is not planar)");

    std::map<std::pair<int, int>, int> pairFlowEdge;
    for (int f = 0; f < numFaces; ++f) {
        for (const int h : faces[f]) {
            const int g = faceOf[-h];
            if (g != f && pairFlowEdge.count({f, g}) == 0) {
                pairFlowEdge[{f, g}] = flow.AddEdge(f, g, kInf, 1);
            }
        }
    }

    const int sent = flow.Solve(source, sink);
    Require(sent == sourceDemand, "bend flow is infeasible (diagram is not planar)");

    // ---- Convert flow to bends. All flow between a face pair runs through
    // one shared edge; bendTurns[k] is recorded along the +k direction. ----

    std::vector<std::vector<int>> bendTurns(n + 1);
    std::set<std::pair<int, int>> handledPairs;
    for (int f = 0; f < numFaces; ++f) {
        for (const int h : faces[f]) {
            const int g = faceOf[-h];
            if (g == f) {
                continue;
            }
            const std::pair<int, int> key = {std::min(f, g), std::max(f, g)};
            if (!handledPairs.insert(key).second) {
                continue;
            }
            int wOut = flow.FlowOn(pairFlowEdge[{f, g}]);
            int wIn = flow.FlowOn(pairFlowEdge[{g, f}]);
            const int cancel = std::min(wOut, wIn);
            wOut -= cancel;
            wIn -= cancel;
            if (wOut + wIn == 0) {
                continue;
            }
            // Turns along h as traversed by face f: convex (+1) bends for
            // outgoing flow, reflex (-1) for incoming.
            std::vector<int> turns(wOut, 1);
            turns.insert(turns.end(), wIn, -1);
            const int k = std::abs(h);
            if (h < 0) {
                std::reverse(turns.begin(), turns.end());
                for (int& t : turns) {
                    t = -t;
                }
            }
            Require(bendTurns[k].empty(), "edge bent twice");
            bendTurns[k] = std::move(turns);
        }
    }

    // ---- Subdivide edges at bends; build segments. ----

    int numVertices = numCrossings;
    std::vector<std::vector<int>> bendVerts(n + 1);
    for (int k = 1; k <= n; ++k) {
        for (size_t i = 0; i < bendTurns[k].size(); ++i) {
            bendVerts[k].push_back(numVertices++);
        }
    }

    struct Seg
    {
        int u;
        int v;
    };
    std::vector<Seg> segs;
    std::vector<std::vector<int>> edgeSegs(n + 1);
    for (int k = 1; k <= n; ++k) {
        int prev = crossingOf[k];
        for (const int bend : bendVerts[k]) {
            edgeSegs[k].push_back(static_cast<int>(segs.size()));
            segs.push_back({prev, bend});
            prev = bend;
        }
        edgeSegs[k].push_back(static_cast<int>(segs.size()));
        segs.push_back({prev, crossingOf[-k]});
    }

    // ---- Face boundary walks over segments, with turns. ----

    struct WalkElem
    {
        int seg;
        bool forward;
    };
    std::vector<std::vector<WalkElem>> walkElems(numFaces);
    std::vector<std::vector<int>> walkTurns(numFaces);
    for (int f = 0; f < numFaces; ++f) {
        for (const int h : faces[f]) {
            const int k = std::abs(h);
            const std::vector<int>& ss = edgeSegs[k];
            const std::vector<int>& bt = bendTurns[k];
            const int m = static_cast<int>(bt.size());
            if (h > 0) {
                for (int i = 0; i <= m; ++i) {
                    walkElems[f].push_back({ss[i], true});
                    walkTurns[f].push_back(i < m ? bt[i] : 1);
                }
            }
            else {
                for (int i = m; i >= 0; --i) {
                    walkElems[f].push_back({ss[i], false});
                    walkTurns[f].push_back(i > 0 ? -bt[i - 1] : 1);
                }
            }
        }
        int rotation = 0;
        for (const int t : walkTurns[f]) {
            rotation += t;
        }
        Require(rotation == (f == exteriorFace ? -4 : 4), "face rotation after bending");
    }

    // ---- Propagate edge directions across faces. ----

    std::vector<int> segDir(segs.size(), -1); // direction along the canonical (+k) orientation
    const auto canonical = [](int dir, bool forward) {
        return forward ? dir : (dir + 2) % 4;
    };
    std::vector<std::vector<int>> segFaces(segs.size());
    for (int f = 0; f < numFaces; ++f) {
        for (const WalkElem& elem : walkElems[f]) {
            segFaces[elem.seg].push_back(f);
        }
    }

    segDir[walkElems[0][0].seg] = canonical(kRight, walkElems[0][0].forward);
    std::deque<int> faceQueue = {0};
    std::vector<char> queued(numFaces, 0);
    queued[0] = 1;
    int processedFaces = 0;
    while (!faceQueue.empty()) {
        const int f = faceQueue.front();
        faceQueue.pop_front();
        ++processedFaces;
        const std::vector<WalkElem>& elems = walkElems[f];
        const int size = static_cast<int>(elems.size());
        int start = -1;
        int dir = -1;
        for (int i = 0; i < size; ++i) {
            if (segDir[elems[i].seg] != -1) {
                start = i;
                dir = canonical(segDir[elems[i].seg], elems[i].forward);
                break;
            }
        }
        Require(start != -1, "face with no oriented edge");
        for (int j = 0; j < size; ++j) {
            const int i = (start + j) % size;
            const int canon = canonical(dir, elems[i].forward);
            if (segDir[elems[i].seg] == -1) {
                segDir[elems[i].seg] = canon;
            }
            else {
                Require(segDir[elems[i].seg] == canon, "inconsistent edge orientation");
            }
            for (const int g : segFaces[elems[i].seg]) {
                if (!queued[g]) {
                    queued[g] = 1;
                    faceQueue.push_back(g);
                }
            }
            dir = ((dir + walkTurns[f][i]) % 4 + 4) % 4;
        }
    }
    Require(processedFaces == numFaces, "diagram is disconnected");

    // ---- Build the orthogonal representation and compact it. ----

    OrthoRep rep(numVertices);
    for (size_t s = 0; s < segs.size(); ++s) {
        switch (segDir[s]) {
            case kRight: rep.AddEdge(segs[s].u, segs[s].v, true, false); break;
            case kLeft: rep.AddEdge(segs[s].v, segs[s].u, true, false); break;
            case kUp: rep.AddEdge(segs[s].u, segs[s].v, false, false); break;
            case kDown: rep.AddEdge(segs[s].v, segs[s].u, false, false); break;
            default: Require(false, "unoriented segment");
        }
    }
    rep.MakeTurnRegular();
    std::vector<int> ys = rep.ChainCoordinates(true);
    std::vector<int> xs = rep.ChainCoordinates(false);

    int minX = kInf;
    int minY = kInf;
    for (int v = 0; v < numVertices; ++v) {
        minX = std::min(minX, xs[v]);
        minY = std::min(minY, ys[v]);
    }
    for (int v = 0; v < numVertices; ++v) {
        xs[v] -= minX;
        ys[v] -= minY;
    }

    // Mirror horizontally so that sigma cycles wind counterclockwise when the
    // drawing is rendered with y pointing up.
    int maxX = 0;
    for (int v = 0; v < numVertices; ++v) {
        maxX = std::max(maxX, xs[v]);
    }
    for (int v = 0; v < numVertices; ++v) {
        xs[v] = maxX - xs[v];
    }

    std::set<std::pair<int, int>> seenPositions;
    for (int v = 0; v < numVertices; ++v) {
        Require(seenPositions.insert({xs[v], ys[v]}).second, "two vertices share a grid point");
    }

    // ---- Assemble the result. ----

    GridLayout layout;
    for (int v = 0; v < numVertices; ++v) {
        layout.width = std::max(layout.width, xs[v]);
        layout.height = std::max(layout.height, ys[v]);
    }
    layout.crossingPositions.resize(numCrossings);
    for (int c = 0; c < numCrossings; ++c) {
        layout.crossingPositions[c] = {xs[c], ys[c]};
    }
    layout.edgePaths.resize(n);
    for (int k = 1; k <= n; ++k) {
        std::vector<GridPoint>& path = layout.edgePaths[k - 1];
        path.push_back({xs[crossingOf[k]], ys[crossingOf[k]]});
        for (const int bend : bendVerts[k]) {
            path.push_back({xs[bend], ys[bend]});
        }
        path.push_back({xs[crossingOf[-k]], ys[crossingOf[-k]]});
    }
    return layout;
}

std::vector<Vec2d> GridFacePolygon(const Multiloop& loop, const GridLayout& layout, int face)
{
    std::vector<Vec2d> poly;
    if (face < 0 || face >= static_cast<int>(loop.faces.size())) {
        return poly;
    }
    std::vector<GridPoint> gridPoly;
    const auto append = [&](const GridPoint& p) {
        if (!gridPoly.empty() && gridPoly.back().x == p.x && gridPoly.back().y == p.y) {
            return;
        }
        gridPoly.push_back(p);
    };
    for (const int h : loop.faces[face]) {
        const int k = std::abs(h);
        if (k < 1 || k > loop.numEdges) {
            continue;
        }
        const std::vector<GridPoint>& path = layout.edgePaths[static_cast<size_t>(k - 1)];
        if (path.empty()) {
            continue;
        }
        if (h > 0) {
            for (size_t i = 0; i < path.size(); ++i) {
                append(path[i]);
            }
        }
        else {
            for (int i = static_cast<int>(path.size()) - 1; i >= 0; --i) {
                append(path[static_cast<size_t>(i)]);
            }
        }
    }
    if (gridPoly.size() >= 2 && gridPoly.front().x == gridPoly.back().x &&
        gridPoly.front().y == gridPoly.back().y) {
        gridPoly.pop_back();
    }
    poly.reserve(gridPoly.size());
    for (const GridPoint& p : gridPoly) {
        poly.push_back(Vec2d{static_cast<double>(p.x), static_cast<double>(p.y)});
    }
    return poly;
}

double GridFaceArea(const Multiloop& loop, const GridLayout& layout, int face)
{
    return ShoelaceArea(GridFacePolygon(loop, layout, face));
}

namespace
{
    constexpr int kPackGap = 2;

    GridPoint OffsetPoint(GridPoint p, int dx, int dy)
    {
        p.x += dx;
        p.y += dy;
        return p;
    }

    int LocalExteriorFace(const Multiloop& globalLoop, const Multiloop& localLoop,
                          const DenseRelabel& relabel, const MultiloopComponent& component)
    {
        int globalFace = component.exteriorFace;
        if (globalFace < 0 || globalFace >= static_cast<int>(globalLoop.faces.size()) ||
            !ComponentOwnsFace(component, globalLoop, globalFace)) {
            globalFace = FirstFaceOnComponent(globalLoop, component.halfEdges);
        }
        if (globalFace < 0 || globalFace >= static_cast<int>(globalLoop.faces.size())) {
            return 0;
        }
        for (const int h : globalLoop.faces[static_cast<size_t>(globalFace)]) {
            const auto it = relabel.globalToLocal.find(h);
            if (it == relabel.globalToLocal.end()) {
                continue;
            }
            const int localFace = FaceIndexContainingHalfEdge(localLoop, it->second);
            if (localFace >= 0) {
                return localFace;
            }
        }
        return 0;
    }
}

PackedGridLayout ComputePackedGridLayout(const Multiloop& loop,
                                         const std::vector<MultiloopComponent>& componentsIn)
{
    std::vector<MultiloopComponent> components = componentsIn;
    if (components.empty() && loop.numEdges > 0) {
        MultiloopComponent part;
        part.kind = MultiloopPart::Crossings;
        part.halfEdges.reserve(static_cast<size_t>(2 * loop.numEdges));
        for (int k = 1; k <= loop.numEdges; ++k) {
            part.halfEdges.push_back(k);
            part.halfEdges.push_back(-k);
        }
        part.exteriorFace = 0;
        components.push_back(std::move(part));
    }

    PackedGridLayout packed;
    packed.componentOrigin.assign(components.size(), {});
    packed.componentWidth.assign(components.size(), 0);
    packed.componentHeight.assign(components.size(), 0);
    if (loop.numEdges > 0) {
        packed.layout.edgePaths.resize(static_cast<size_t>(loop.numEdges));
        packed.layout.crossingPositions.resize(loop.sigma.Cycles().size());
    }

    std::unordered_map<int, int> globalCrossingOf;
    const std::vector<std::vector<int>> globalCycles = loop.sigma.Cycles();
    for (int c = 0; c < static_cast<int>(globalCycles.size()); ++c) {
        for (const int h : globalCycles[static_cast<size_t>(c)]) {
            globalCrossingOf[h] = c;
        }
    }

    int cursor = 0;
    int maxHeight = 0;
    for (size_t i = 0; i < components.size(); ++i) {
        if (i > 0) {
            cursor += kPackGap;
        }
        const MultiloopComponent& component = components[i];
        int width = 1;
        int height = 1;
        if (component.kind == MultiloopPart::Circle) {
            packed.componentOrigin[i] = {cursor, 0};
        }
        else {
            const DenseRelabel relabel = RelabelComponent(loop.sigma, component.halfEdges);
            if (relabel.localSigma.MaxAbsLabel() <= 0) {
                throw std::runtime_error("crossing component has no edges");
            }
            const Multiloop localLoop(relabel.localSigma);
            const int localExterior = LocalExteriorFace(loop, localLoop, relabel, component);
            GridLayout local = ComputeGridLayout(localLoop, localExterior);
            width = local.width;
            height = local.height;
            packed.componentOrigin[i] = {cursor, 0};
            for (int k = 1; k <= localLoop.numEdges; ++k) {
                const auto git = relabel.localToGlobal.find(k);
                if (git == relabel.localToGlobal.end()) {
                    continue;
                }
                const int globalK = git->second;
                if (globalK < 1 || globalK > loop.numEdges) {
                    continue;
                }
                std::vector<GridPoint> path = local.edgePaths[static_cast<size_t>(k - 1)];
                for (GridPoint& p : path) {
                    p = OffsetPoint(p, cursor, 0);
                }
                packed.layout.edgePaths[static_cast<size_t>(globalK - 1)] = std::move(path);
            }
            const std::vector<std::vector<int>> localCycles = localLoop.sigma.Cycles();
            for (int lc = 0; lc < static_cast<int>(localCycles.size()); ++lc) {
                if (localCycles[static_cast<size_t>(lc)].empty()) {
                    continue;
                }
                const int lh = localCycles[static_cast<size_t>(lc)].front();
                const auto git = relabel.localToGlobal.find(lh);
                if (git == relabel.localToGlobal.end()) {
                    continue;
                }
                const auto cit = globalCrossingOf.find(git->second);
                if (cit == globalCrossingOf.end()) {
                    continue;
                }
                packed.layout.crossingPositions[static_cast<size_t>(cit->second)] =
                    OffsetPoint(local.crossingPositions[static_cast<size_t>(lc)], cursor, 0);
            }
        }
        packed.componentWidth[i] = width;
        packed.componentHeight[i] = height;
        cursor += width;
        maxHeight = std::max(maxHeight, height);
    }

    packed.layout.width = (components.empty()) ? 0 : cursor;
    packed.layout.height = maxHeight;
    return packed;
}
