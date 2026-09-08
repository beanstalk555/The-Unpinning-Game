#include "move_rewrite.h"

#include "geometry.h"

#include <algorithm>
#include <cstdlib>
#include <exception>
#include <set>
#include <unordered_map>
#include <unordered_set>

int RewriteResult::MapHalfEdge(int h) const
{
    const auto it = oldToNew.find(h);
    return (it == oldToNew.end()) ? 0 : it->second;
}

int RewriteResult::MapEdge(int k) const
{
    return std::abs(MapHalfEdge(std::abs(k)));
}

bool DeleteCrossingsAndSplice(const Multiloop& old, const std::vector<int>& deletedCrossings,
                              const std::vector<int>& deletedEdges,
                              const std::vector<SplicePair>& splices, RewriteResult* out)
{
    if (out == nullptr) {
        return false;
    }
    const std::vector<std::vector<int>> cycles = old.sigma.Cycles();
    std::unordered_map<int, int> crossingOf;
    for (int c = 0; c < static_cast<int>(cycles.size()); ++c) {
        for (const int h : cycles[static_cast<size_t>(c)]) {
            crossingOf[h] = c;
        }
    }
    std::unordered_set<int> deadCrossing(deletedCrossings.begin(), deletedCrossings.end());
    std::unordered_set<int> deadEdge;
    for (const int k : deletedEdges) {
        deadEdge.insert(std::abs(k));
    }
    for (const int c : deadCrossing) {
        if (c < 0 || c >= static_cast<int>(cycles.size())) {
            return false;
        }
    }
    const auto isDeadCrossingOf = [&](int h) {
        const auto it = crossingOf.find(h);
        return it != crossingOf.end() && deadCrossing.count(it->second) != 0;
    };
    // A deleted edge must have both ends at deleted crossings.
    for (const int k : deadEdge) {
        if (k < 1 || k > old.numEdges || !isDeadCrossingOf(k) || !isDeadCrossingOf(-k)) {
            return false;
        }
    }

    std::unordered_map<int, int> partner;
    for (const SplicePair& sp : splices) {
        if (sp.u == 0 || sp.v == 0 || sp.u == sp.v || deadEdge.count(std::abs(sp.u)) != 0 ||
            deadEdge.count(std::abs(sp.v)) != 0 || !isDeadCrossingOf(sp.u) ||
            !isDeadCrossingOf(sp.v)) {
            return false;
        }
        if (partner.count(sp.u) != 0 || partner.count(sp.v) != 0) {
            return false;
        }
        partner[sp.u] = sp.v;
        partner[sp.v] = sp.u;
    }
    // Every surviving edge end at a deleted crossing needs a partner.
    for (int k = 1; k <= old.numEdges; ++k) {
        if (deadEdge.count(k) != 0) {
            continue;
        }
        for (const int h : {k, -k}) {
            if (isDeadCrossingOf(h) && partner.count(h) == 0) {
                return false;
            }
        }
    }

    std::vector<char> visited(static_cast<size_t>(old.numEdges + 1), 0);
    std::vector<std::vector<int>> arcs;
    std::vector<std::vector<int>> circles;

    // Traverse from `start` (half-edge at the start of the first edge). Stops
    // when the far end sits at a surviving crossing, or when it closes.
    const auto trace = [&](int start, bool* closed) -> std::vector<int> {
        std::vector<int> arc;
        int cur = start;
        const int limit = old.numEdges + 2;
        for (int guard = 0; guard < limit; ++guard) {
            const int k = std::abs(cur);
            if (visited[static_cast<size_t>(k)]) {
                arc.clear(); // revisiting: inconsistent splice data
                return arc;
            }
            visited[static_cast<size_t>(k)] = 1;
            arc.push_back(cur);
            const int far = -cur;
            if (!isDeadCrossingOf(far)) {
                *closed = false;
                return arc;
            }
            const int next = partner.at(far);
            if (next == start) {
                *closed = true;
                return arc;
            }
            cur = next;
        }
        arc.clear();
        return arc;
    };

    for (int c = 0; c < static_cast<int>(cycles.size()); ++c) {
        if (deadCrossing.count(c) != 0) {
            continue;
        }
        for (const int h : cycles[static_cast<size_t>(c)]) {
            if (deadEdge.count(std::abs(h)) != 0) {
                return false; // a surviving crossing on a deleted edge
            }
            if (visited[static_cast<size_t>(std::abs(h))]) {
                continue;
            }
            bool closed = false;
            std::vector<int> arc = trace(h, &closed);
            if (arc.empty() || closed) {
                return false;
            }
            arcs.push_back(std::move(arc));
        }
    }
    for (int k = 1; k <= old.numEdges; ++k) {
        if (deadEdge.count(k) != 0 || visited[static_cast<size_t>(k)]) {
            continue;
        }
        bool closed = false;
        std::vector<int> arc = trace(k, &closed);
        if (arc.empty() || !closed) {
            return false;
        }
        circles.push_back(std::move(arc));
    }

    // Orient each arc along its smallest constituent label.
    struct ArcRec
    {
        int rep;
        std::vector<int> arc;
    };
    std::vector<ArcRec> recs;
    recs.reserve(arcs.size());
    for (std::vector<int>& arc : arcs) {
        int rep = std::abs(arc.front());
        int repSign = arc.front() > 0 ? 1 : -1;
        for (const int h : arc) {
            if (std::abs(h) < rep) {
                rep = std::abs(h);
                repSign = h > 0 ? 1 : -1;
            }
        }
        if (repSign < 0) {
            std::reverse(arc.begin(), arc.end());
            for (int& h : arc) {
                h = -h;
            }
        }
        recs.push_back({rep, std::move(arc)});
    }
    std::sort(recs.begin(), recs.end(),
              [](const ArcRec& a, const ArcRec& b) { return a.rep < b.rep; });

    RewriteResult result;
    result.numNewEdges = static_cast<int>(recs.size());
    result.mergedArcs.assign(static_cast<size_t>(result.numNewEdges + 1), {});
    for (int m = 1; m <= result.numNewEdges; ++m) {
        const std::vector<int>& arc = recs[static_cast<size_t>(m - 1)].arc;
        for (const int h : arc) {
            result.oldToNew[h] = m;
            result.oldToNew[-h] = -m;
        }
        result.mergedArcs[static_cast<size_t>(m)] = arc;
    }
    result.deletedEdges.assign(deadEdge.begin(), deadEdge.end());
    std::sort(result.deletedEdges.begin(), result.deletedEdges.end());
    for (size_t ci = 0; ci < circles.size(); ++ci) {
        const std::vector<int>& arc = circles[ci];
        for (const int h : arc) {
            result.circleOfHalfEdge[h] = static_cast<int>(ci);
            result.circleOfHalfEdge[-h] = static_cast<int>(ci);
        }
        const int k = std::abs(arc.front());
        result.circleStrand.push_back(
            (k >= 1 && k < static_cast<int>(old.strandOfEdge.size())) ? old.strandOfEdge[k] : -1);
    }
    result.circleArcs = std::move(circles);

    std::vector<std::vector<int>> newCycles;
    for (int c = 0; c < static_cast<int>(cycles.size()); ++c) {
        if (deadCrossing.count(c) != 0) {
            continue;
        }
        std::vector<int> cycle;
        for (const int h : cycles[static_cast<size_t>(c)]) {
            const auto it = result.oldToNew.find(h);
            if (it == result.oldToNew.end()) {
                return false;
            }
            cycle.push_back(it->second);
        }
        newCycles.push_back(std::move(cycle));
    }
    try {
        result.sigma = Permutation::FromCycles(newCycles);
    }
    catch (const std::exception&) {
        return false;
    }
    if (!ValidateMultiloopSigma(result.sigma).empty()) {
        return false;
    }
    if (result.sigma.MaxAbsLabel() != result.numNewEdges) {
        return false;
    }
    *out = std::move(result);
    return true;
}

void RemapStrandsThroughRewrite(const Multiloop& oldLoop, const RewriteResult& rw,
                                Multiloop* newLoop)
{
    if (newLoop == nullptr) {
        return;
    }
    struct Cand
    {
        int overlap;
        int newStrand;
        int oldStrand;
    };
    std::vector<Cand> cands;
    for (int ns = 0; ns < newLoop->numStrands; ++ns) {
        std::unordered_map<int, int> votes;
        for (int m = 1; m <= newLoop->numEdges; ++m) {
            if (newLoop->strandOfEdge[static_cast<size_t>(m)] != ns ||
                static_cast<size_t>(m) >= rw.mergedArcs.size()) {
                continue;
            }
            for (const int h : rw.mergedArcs[static_cast<size_t>(m)]) {
                const int k = std::abs(h);
                if (k >= 1 && k < static_cast<int>(oldLoop.strandOfEdge.size())) {
                    ++votes[oldLoop.strandOfEdge[static_cast<size_t>(k)]];
                }
            }
        }
        for (const auto& kv : votes) {
            cands.push_back({kv.second, ns, kv.first});
        }
    }
    std::sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) {
        if (a.overlap != b.overlap) {
            return a.overlap > b.overlap;
        }
        if (a.newStrand != b.newStrand) {
            return a.newStrand < b.newStrand;
        }
        return a.oldStrand < b.oldStrand;
    });
    std::vector<int> newToOld(static_cast<size_t>(newLoop->numStrands), -1);
    std::unordered_set<int> taken;
    for (const Cand& c : cands) {
        if (newToOld[static_cast<size_t>(c.newStrand)] >= 0 || taken.count(c.oldStrand) != 0) {
            continue;
        }
        newToOld[static_cast<size_t>(c.newStrand)] = c.oldStrand;
        taken.insert(c.oldStrand);
    }
    int fresh = std::max(oldLoop.numStrands, 0);
    for (int ns = 0; ns < newLoop->numStrands; ++ns) {
        if (newToOld[static_cast<size_t>(ns)] >= 0) {
            continue;
        }
        int pick = -1;
        for (int os = 0; os < oldLoop.numStrands; ++os) {
            if (taken.count(os) == 0) {
                pick = os;
                break;
            }
        }
        if (pick < 0) {
            pick = fresh++;
        }
        newToOld[static_cast<size_t>(ns)] = pick;
        taken.insert(pick);
    }
    int maxId = -1;
    for (int m = 1; m <= newLoop->numEdges; ++m) {
        int& s = newLoop->strandOfEdge[static_cast<size_t>(m)];
        if (s >= 0 && s < static_cast<int>(newToOld.size())) {
            s = newToOld[static_cast<size_t>(s)];
        }
        maxId = std::max(maxId, s);
    }
    newLoop->numStrands = std::max(oldLoop.numStrands, maxId + 1);
}

int CircleClosingOldFace(const Multiloop& oldLoop, const RewriteResult& rw, int oldFace)
{
    if (oldFace < 0 || oldFace >= static_cast<int>(oldLoop.faces.size())) {
        return -1;
    }
    std::set<int> circles;
    for (const int h : oldLoop.faces[static_cast<size_t>(oldFace)]) {
        const auto it = rw.circleOfHalfEdge.find(h);
        if (it != rw.circleOfHalfEdge.end()) {
            circles.insert(it->second);
        }
    }
    return circles.size() == 1 ? *circles.begin() : -1;
}

std::vector<int> MappedFaceHalfEdges(const Multiloop& oldLoop, const RewriteResult& rw,
                                     int oldFace)
{
    std::vector<int> out;
    if (oldFace < 0 || oldFace >= static_cast<int>(oldLoop.faces.size())) {
        return out;
    }
    for (const int h : oldLoop.faces[static_cast<size_t>(oldFace)]) {
        if (rw.circleOfHalfEdge.count(h) != 0) {
            continue;
        }
        const int m = rw.MapHalfEdge(h);
        if (m != 0) {
            out.push_back(m);
        }
    }
    return out;
}

std::vector<int> MatchFacesByOverlap(const Multiloop& oldLoop, const Multiloop& newLoop,
                                     const RewriteResult& rw, const std::vector<int>& oldFaces,
                                     const std::vector<int>& forbidden)
{
    struct Cand
    {
        int overlap;
        int oldIdx;
        int newFace;
    };
    std::unordered_set<int> forbid(forbidden.begin(), forbidden.end());
    std::vector<Cand> cands;
    for (int oi = 0; oi < static_cast<int>(oldFaces.size()); ++oi) {
        const std::vector<int> mapped = MappedFaceHalfEdges(oldLoop, rw, oldFaces[static_cast<size_t>(oi)]);
        if (mapped.empty()) {
            continue;
        }
        std::unordered_set<int> set(mapped.begin(), mapped.end());
        for (int nf = 0; nf < static_cast<int>(newLoop.faces.size()); ++nf) {
            if (forbid.count(nf) != 0) {
                continue;
            }
            int ov = 0;
            for (const int h : newLoop.faces[static_cast<size_t>(nf)]) {
                if (set.count(h) != 0) {
                    ++ov;
                }
            }
            if (ov > 0) {
                cands.push_back({ov, oi, nf});
            }
        }
    }
    std::sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) {
        if (a.overlap != b.overlap) {
            return a.overlap > b.overlap;
        }
        if (a.oldIdx != b.oldIdx) {
            return a.oldIdx < b.oldIdx;
        }
        return a.newFace < b.newFace;
    });
    std::vector<int> assignment(oldFaces.size(), -1);
    std::vector<char> newTaken(newLoop.faces.size(), 0);
    for (const Cand& c : cands) {
        if (assignment[static_cast<size_t>(c.oldIdx)] >= 0 || newTaken[static_cast<size_t>(c.newFace)]) {
            continue;
        }
        assignment[static_cast<size_t>(c.oldIdx)] = c.newFace;
        newTaken[static_cast<size_t>(c.newFace)] = 1;
    }
    return assignment;
}

namespace
{
    // Undirected step along an edge bead's links.
    int OtherNeighbor(const Bead& b, int prev)
    {
        return (b.neighbors[0] == prev) ? b.neighbors[1] : b.neighbors[0];
    }
}

std::vector<int> CircleRingThrough(const Necklace& neck, int bead)
{
    std::vector<int> ring;
    const int n = static_cast<int>(neck.beads.size());
    if (bead < 0 || bead >= n) {
        return ring;
    }
    const Bead& start = neck.beads[static_cast<size_t>(bead)];
    if (start.isVertex || start.removed) {
        return ring;
    }
    int prev = -1;
    int cur = bead;
    for (int guard = 0; guard < n + 2; ++guard) {
        const Bead& b = neck.beads[static_cast<size_t>(cur)];
        if (b.isVertex || b.removed) {
            ring.clear();
            return ring;
        }
        ring.push_back(cur);
        int nxt = -1;
        if (prev < 0) {
            nxt = (b.neighbors[1] >= 0) ? b.neighbors[1] : b.neighbors[0];
        }
        else {
            nxt = OtherNeighbor(b, prev);
        }
        if (nxt < 0 || nxt >= n) {
            ring.clear();
            return ring;
        }
        prev = cur;
        cur = nxt;
        if (cur == bead) {
            return ring;
        }
    }
    ring.clear();
    return ring;
}

bool ApplyRewriteToNecklace(Necklace& neck, const RewriteResult& rw, const Multiloop& newLoop)
{
    const int n = static_cast<int>(neck.beads.size());
    const std::vector<std::vector<int>> newCycles = newLoop.sigma.Cycles();
    std::unordered_map<int, int> newCrossingOf;
    for (int c = 0; c < static_cast<int>(newCycles.size()); ++c) {
        for (const int h : newCycles[static_cast<size_t>(c)]) {
            newCrossingOf[h] = c;
        }
    }

    // ---- Phase 1: compute everything. ----
    struct VertexPlan
    {
        int bead;
        std::vector<int> newCycle;
        int crossing;
    };
    std::vector<VertexPlan> vertexPlans;
    // New signed half-edge -> (vertex bead, slot).
    std::unordered_map<int, std::pair<int, int>> beadOfNewHalf;
    for (int i = 0; i < n; ++i) {
        const Bead& b = neck.beads[static_cast<size_t>(i)];
        if (!b.isVertex || b.removed) {
            continue;
        }
        if (b.sigmaCycle.size() != 4) {
            return false;
        }
        VertexPlan plan;
        plan.bead = i;
        for (size_t s = 0; s < 4; ++s) {
            const int m = rw.MapHalfEdge(b.sigmaCycle[s]);
            if (m == 0) {
                return false; // vertex bead of a deleted crossing was not converted
            }
            plan.newCycle.push_back(m);
            beadOfNewHalf[m] = {i, static_cast<int>(s)};
        }
        const auto it = newCrossingOf.find(plan.newCycle[0]);
        if (it == newCrossingOf.end()) {
            return false;
        }
        plan.crossing = it->second;
        vertexPlans.push_back(std::move(plan));
    }
    if (static_cast<int>(vertexPlans.size()) != static_cast<int>(newCycles.size())) {
        return false;
    }

    std::vector<std::vector<int>> chains(static_cast<size_t>(newLoop.numEdges + 1));
    std::vector<char> claimed(static_cast<size_t>(n), 0);
    for (int m = 1; m <= newLoop.numEdges; ++m) {
        const auto tail = beadOfNewHalf.find(m);
        const auto head = beadOfNewHalf.find(-m);
        if (tail == beadOfNewHalf.end() || head == beadOfNewHalf.end()) {
            return false;
        }
        const int vb = tail->second.first;
        const int vEnd = head->second.first;
        std::vector<int> chain;
        chain.push_back(vb);
        int prev = vb;
        int cur = neck.beads[static_cast<size_t>(vb)].neighborOfHalfEdge[tail->second.second];
        int guard = 0;
        while (cur >= 0 && cur < n && !neck.beads[static_cast<size_t>(cur)].isVertex) {
            const Bead& b = neck.beads[static_cast<size_t>(cur)];
            if (b.removed || claimed[static_cast<size_t>(cur)]) {
                return false;
            }
            claimed[static_cast<size_t>(cur)] = 1;
            chain.push_back(cur);
            const int nxt = OtherNeighbor(b, prev);
            prev = cur;
            cur = nxt;
            if (++guard > n + 8) {
                return false;
            }
        }
        if (cur != vEnd) {
            return false;
        }
        if (neck.beads[static_cast<size_t>(vEnd)].neighborOfHalfEdge[head->second.second] != prev) {
            return false;
        }
        chain.push_back(vEnd);
        chains[static_cast<size_t>(m)] = std::move(chain);
    }

    std::vector<std::vector<int>> rings;
    for (size_t ci = 0; ci < rw.circleArcs.size(); ++ci) {
        int seed = -1;
        for (const int h : rw.circleArcs[ci]) {
            const int k = std::abs(h);
            for (int i = 0; i < n && seed < 0; ++i) {
                const Bead& b = neck.beads[static_cast<size_t>(i)];
                if (!b.isVertex && !b.removed && b.edge == k && !claimed[static_cast<size_t>(i)]) {
                    seed = i;
                }
            }
            if (seed >= 0) {
                break;
            }
        }
        if (seed < 0) {
            return false;
        }
        std::vector<int> ring = CircleRingThrough(neck, seed);
        if (ring.size() < 3) {
            return false;
        }
        for (const int b : ring) {
            if (claimed[static_cast<size_t>(b)]) {
                return false;
            }
            claimed[static_cast<size_t>(b)] = 1;
        }
        rings.push_back(std::move(ring));
    }

    // ---- Phase 2: mutate. ----
    for (const VertexPlan& plan : vertexPlans) {
        Bead& b = neck.beads[static_cast<size_t>(plan.bead)];
        b.sigmaCycle = plan.newCycle;
        b.crossing = plan.crossing;
    }
    neck.crossingBead.assign(newCycles.size(), -1);
    for (const VertexPlan& plan : vertexPlans) {
        neck.crossingBead[static_cast<size_t>(plan.crossing)] = plan.bead;
    }
    for (int m = 1; m <= newLoop.numEdges; ++m) {
        const std::vector<int>& chain = chains[static_cast<size_t>(m)];
        for (size_t t = 1; t + 1 < chain.size(); ++t) {
            Bead& b = neck.beads[static_cast<size_t>(chain[t])];
            b.edge = m;
            b.neighbors[0] = chain[t - 1];
            b.neighbors[1] = chain[t + 1];
        }
    }
    for (size_t ci = 0; ci < rings.size(); ++ci) {
        const std::vector<int>& ring = rings[ci];
        const int strand = (ci < rw.circleStrand.size()) ? rw.circleStrand[ci] : -1;
        for (size_t t = 0; t < ring.size(); ++t) {
            Bead& b = neck.beads[static_cast<size_t>(ring[t])];
            b.edge = 0;
            b.strand = strand;
            b.neighbors[0] = ring[(t + ring.size() - 1) % ring.size()];
            b.neighbors[1] = ring[(t + 1) % ring.size()];
        }
    }
    // Pre-existing circle rings (edge 0) stay as they are; anything else that
    // no walk reached is an orphan of a deleted edge.
    for (int i = 0; i < n; ++i) {
        Bead& b = neck.beads[static_cast<size_t>(i)];
        if (b.isVertex || b.removed || claimed[static_cast<size_t>(i)] || b.edge == 0) {
            continue;
        }
        b.removed = true;
        b.neighbors[0] = -1;
        b.neighbors[1] = -1;
    }
    return true;
}

namespace
{
    struct Piece
    {
        bool circle = false;
        std::vector<int> halfEdges; // new labels (crossing piece)
        int circleIndex = -1;       // index into rw.circleArcs
        int strand = -1;
        int pieceFace = -1;         // new face holding this piece's part of the split face
        bool hasExterior = false;   // owns the image of the old exterior
        int newExterior = -1;
    };

    int FaceOfNewHalf(const Multiloop& newLoop, int h)
    {
        if (h == 0) {
            return -1;
        }
        return FaceIndexContainingHalfEdge(newLoop, h);
    }

    int RingSeedForStrand(const Necklace& neck, int strand)
    {
        for (size_t i = 0; i < neck.beads.size(); ++i) {
            const Bead& b = neck.beads[i];
            if (!b.isVertex && !b.removed && b.edge == 0 && b.strand == strand) {
                return static_cast<int>(i);
            }
        }
        return -1;
    }

    std::vector<Vec2d> BeadsToPolygon(const Necklace& neck, const std::vector<int>& beads)
    {
        std::vector<Vec2d> poly;
        poly.reserve(beads.size());
        for (const int i : beads) {
            if (i < 0 || i >= static_cast<int>(neck.beads.size()) ||
                neck.beads[static_cast<size_t>(i)].removed) {
                continue;
            }
            const BeadPos& p = neck.beads[static_cast<size_t>(i)].pos;
            poly.push_back(Vec2d{static_cast<double>(p.x), static_cast<double>(p.y)});
        }
        return poly;
    }

    std::vector<Vec2d> FacePolygon(const Multiloop& newLoop, const Necklace& neck, int face)
    {
        if (face < 0 || face >= static_cast<int>(newLoop.faces.size())) {
            return {};
        }
        const std::vector<std::vector<int>> cycles = newLoop.sigma.Cycles();
        std::unordered_map<int, int> crossingOf;
        for (int c = 0; c < static_cast<int>(cycles.size()); ++c) {
            for (const int h : cycles[static_cast<size_t>(c)]) {
                crossingOf[h] = c;
            }
        }
        std::vector<int> boundary;
        std::vector<float> dots;
        WalkFaceBoundaryFromPhiCycle(newLoop.faces[static_cast<size_t>(face)], neck, crossingOf,
                                     boundary, dots, nullptr, 0.0f, false);
        return BeadsToPolygon(neck, boundary);
    }

    std::vector<Vec2d> PiecePolygon(const Piece& piece, const Multiloop& newLoop,
                                    const Necklace& neck, const RewriteResult& rw)
    {
        if (piece.circle) {
            const int strand = (piece.circleIndex >= 0 &&
                                piece.circleIndex < static_cast<int>(rw.circleStrand.size()))
                                   ? rw.circleStrand[static_cast<size_t>(piece.circleIndex)]
                                   : piece.strand;
            return BeadsToPolygon(neck, CircleRingThrough(neck, RingSeedForStrand(neck, strand)));
        }
        const int face = piece.pieceFace >= 0 ? piece.pieceFace : piece.newExterior;
        return FacePolygon(newLoop, neck, face);
    }

    Vec2d PolygonCentroid(const std::vector<Vec2d>& poly)
    {
        Vec2d c;
        if (poly.empty()) {
            return c;
        }
        for (const Vec2d& p : poly) {
            c.x += p.x;
            c.y += p.y;
        }
        c.x /= static_cast<double>(poly.size());
        c.y /= static_cast<double>(poly.size());
        return c;
    }

    bool SampleComponentPoint(const MultiloopComponent& c, const Necklace& neck,
                              const RewriteResult& rw, Vec2d* out)
    {
        if (out == nullptr) {
            return false;
        }
        if (c.kind == MultiloopPart::Circle) {
            const std::vector<int> ring = CircleRingThrough(neck, RingSeedForStrand(neck, c.strand));
            const std::vector<Vec2d> poly = BeadsToPolygon(neck, ring);
            if (poly.size() < 3) {
                return false;
            }
            *out = PolygonCentroid(poly);
            return true;
        }
        for (const Bead& b : neck.beads) {
            if (b.removed || !b.isVertex) {
                continue;
            }
            bool hit = false;
            for (const int oldH : c.halfEdges) {
                const int m = rw.MapHalfEdge(oldH);
                if (m == 0) {
                    continue;
                }
                if (std::find(b.sigmaCycle.begin(), b.sigmaCycle.end(), m) != b.sigmaCycle.end()) {
                    hit = true;
                    break;
                }
            }
            if (hit) {
                *out = Vec2d{static_cast<double>(b.pos.x), static_cast<double>(b.pos.y)};
                return true;
            }
        }
        return false;
    }

    void NestInPiece(MultiloopComponent* c, const Piece& piece)
    {
        if (c == nullptr) {
            return;
        }
        if (piece.circle) {
            c->parentFace = -1;
            c->parentCircle = piece.strand;
        }
        else {
            c->parentFace = piece.pieceFace;
            c->parentCircle = -1;
        }
    }

    void MakePad(FaceBead& fb, const std::vector<int>& ring)
    {
        fb.face = -1;
        fb.phiCycle.clear();
        fb.vertexBeads.clear();
        fb.boundaryBeads = ring;
        fb.boundaryDots.clear();
        fb.holeBeads.clear();
        fb.cachedExcessArea = 0.0;
    }
}

bool FinishMoveBookkeeping(const Multiloop& oldLoop, const Multiloop& newLoop,
                           const RewriteResult& rw, Necklace& neck,
                           const std::vector<MultiloopComponent>& oldComponents,
                           const MoveBookkeepingInput& in, MoveComponentResult* out)
{
    if (out == nullptr || in.movedComponent < 0 ||
        in.movedComponent >= static_cast<int>(oldComponents.size())) {
        return false;
    }
    const MultiloopComponent& moved = oldComponents[static_cast<size_t>(in.movedComponent)];
    if (moved.kind != MultiloopPart::Crossings) {
        return false;
    }
    MoveComponentResult result;

    // ---- Pieces of the moved component. ----
    std::vector<Piece> pieces;
    {
        std::unordered_set<int> seen;
        for (const int h : moved.halfEdges) {
            const int m = rw.MapHalfEdge(h);
            if (m == 0 || rw.circleOfHalfEdge.count(h) != 0 || seen.count(m) != 0) {
                continue;
            }
            Piece piece;
            piece.halfEdges = CollectSigmaComponent(newLoop.sigma, m);
            for (const int x : piece.halfEdges) {
                seen.insert(x);
            }
            pieces.push_back(std::move(piece));
        }
        for (size_t ci = 0; ci < rw.circleArcs.size(); ++ci) {
            Piece piece;
            piece.circle = true;
            piece.circleIndex = static_cast<int>(ci);
            piece.strand = rw.circleStrand[ci];
            pieces.push_back(std::move(piece));
        }
    }
    if (pieces.empty() || pieces.size() > 2) {
        return false;
    }
    result.disconnected = pieces.size() == 2;

    const auto pieceOfOldHalf = [&](int h) -> int {
        if (h == 0) {
            return -1;
        }
        const auto circ = rw.circleOfHalfEdge.find(h);
        for (int p = 0; p < static_cast<int>(pieces.size()); ++p) {
            const Piece& piece = pieces[static_cast<size_t>(p)];
            if (piece.circle) {
                if (circ != rw.circleOfHalfEdge.end() && circ->second == piece.circleIndex) {
                    return p;
                }
                continue;
            }
            const int m = rw.MapHalfEdge(h);
            if (m != 0 && std::find(piece.halfEdges.begin(), piece.halfEdges.end(), m) !=
                              piece.halfEdges.end()) {
                return p;
            }
        }
        return -1;
    };

    // ---- Face map. ----
    const int numOldFaces = static_cast<int>(oldLoop.faces.size());
    result.faceMap.assign(static_cast<size_t>(numOldFaces), -1);
    std::vector<int> forbidden;
    int mergedFace = -1;
    if (in.mergeFaceA >= 0 && in.mergeFaceB >= 0) {
        const int ma = rw.MapHalfEdge(in.pieceHalfA);
        const int mb = rw.MapHalfEdge(in.pieceHalfB);
        // When both survivors closed into a circle the merged region is that
        // circle's pad (handled with the face beads below), not a face.
        if (ma != 0 || mb != 0) {
            mergedFace = FaceOfNewHalf(newLoop, ma);
            if (mergedFace < 0 || mergedFace != FaceOfNewHalf(newLoop, mb)) {
                return false;
            }
            result.faceMap[static_cast<size_t>(in.mergeFaceA)] = mergedFace;
            result.faceMap[static_cast<size_t>(in.mergeFaceB)] = mergedFace;
            forbidden.push_back(mergedFace);
        }
    }
    if (in.splitFace >= 0) {
        if (!result.disconnected) {
            return false;
        }
        const int pa = pieceOfOldHalf(in.pieceHalfA);
        const int pb = pieceOfOldHalf(in.pieceHalfB);
        if (pa < 0 || pb < 0 || pa == pb) {
            return false;
        }
        pieces[static_cast<size_t>(pa)].pieceFace = FaceOfNewHalf(newLoop, rw.MapHalfEdge(in.pieceHalfA));
        pieces[static_cast<size_t>(pb)].pieceFace = FaceOfNewHalf(newLoop, rw.MapHalfEdge(in.pieceHalfB));
        for (const Piece& piece : pieces) {
            if (piece.pieceFace >= 0) {
                forbidden.push_back(piece.pieceFace);
            }
        }
    }
    else if (result.disconnected) {
        return false;
    }
    {
        std::vector<int> oldFaces;
        for (int f = 0; f < numOldFaces; ++f) {
            if (f == in.splitFace || f == in.mergeFaceA || f == in.mergeFaceB) {
                continue;
            }
            oldFaces.push_back(f);
        }
        const std::vector<int> matched =
            MatchFacesByOverlap(oldLoop, newLoop, rw, oldFaces, forbidden);
        for (size_t i = 0; i < oldFaces.size(); ++i) {
            result.faceMap[static_cast<size_t>(oldFaces[i])] = matched[i];
        }
    }

    // ---- Which piece owns the old exterior; nesting of the other. ----
    const int oldExterior = moved.exteriorFace;
    int outer = 0;
    int inner = -1;
    bool exteriorIsSplit = false;
    if (oldExterior >= 0 && oldExterior == in.splitFace) {
        exteriorIsSplit = true;
    }
    else if (oldExterior >= 0) {
        int owner = -1;
        for (const int h : oldLoop.faces[static_cast<size_t>(oldExterior)]) {
            owner = pieceOfOldHalf(h);
            if (owner >= 0) {
                break;
            }
        }
        if (owner < 0) {
            return false;
        }
        outer = owner;
        if (result.disconnected) {
            inner = 1 - outer;
        }
    }

    // Polygons of the two new pieces, used to re-parent holes of the split
    // face that sit inside the new inner piece (concentric circles). The
    // inner/outer roles themselves stay combinatorial: the piece that owns
    // the old exterior is outer. A geometric swap here is unreliable while
    // the just-spliced rings still overlap.
    std::vector<std::vector<Vec2d>> piecePoly(pieces.size());
    if (result.disconnected && pieces.size() == 2) {
        for (size_t p = 0; p < pieces.size(); ++p) {
            piecePoly[p] = PiecePolygon(pieces[p], newLoop, neck, rw);
        }
    }

    for (int p = 0; p < static_cast<int>(pieces.size()); ++p) {
        Piece& piece = pieces[static_cast<size_t>(p)];
        if (piece.circle) {
            piece.newExterior = -1;
            continue;
        }
        if (exteriorIsSplit) {
            piece.newExterior = piece.pieceFace;
        }
        else if (p == outer) {
            piece.newExterior = (oldExterior >= 0) ? result.faceMap[static_cast<size_t>(oldExterior)] : -1;
            if (piece.newExterior < 0) {
                piece.newExterior = FirstFaceOnComponent(newLoop, piece.halfEdges);
            }
        }
        else {
            piece.newExterior = piece.pieceFace;
        }
        if (piece.newExterior < 0) {
            return false;
        }
    }

    // Region that the split face's bead (and previously nested children that
    // are not inside the new inner piece) end up in.
    int splitParentFace = -1;
    int splitParentCircle = -1;
    if (result.disconnected && !exteriorIsSplit) {
        const Piece& o = pieces[static_cast<size_t>(outer)];
        if (o.circle) {
            splitParentCircle = o.strand;
        }
        else {
            splitParentFace = o.pieceFace;
        }
    }

    // ---- Exterior faces of every component (to erase beads landing there). ----
    std::unordered_set<int> newExteriors;
    for (const Piece& piece : pieces) {
        if (!piece.circle && piece.newExterior >= 0) {
            newExteriors.insert(piece.newExterior);
        }
    }
    for (size_t i = 0; i < oldComponents.size(); ++i) {
        if (static_cast<int>(i) == in.movedComponent) {
            continue;
        }
        const MultiloopComponent& c = oldComponents[i];
        if (c.kind == MultiloopPart::Crossings && c.exteriorFace >= 0 &&
            c.exteriorFace < numOldFaces) {
            const int nf = result.faceMap[static_cast<size_t>(c.exteriorFace)];
            if (nf >= 0) {
                newExteriors.insert(nf);
            }
        }
    }

    // ---- Face beads. ----
    bool splitHighlighted = false;
    int splitOuterCircle = -1;
    bool splitBeadKeptAsPad = false;
    if (in.splitFace >= 0) {
        for (const FaceBead& fb : neck.faceBeads) {
            if (fb.face == in.splitFace && fb.highlighted) {
                splitHighlighted = true;
                break;
            }
        }
        if (result.disconnected) {
            const Piece& o = pieces[static_cast<size_t>(outer)];
            if (o.circle) {
                splitOuterCircle = o.circleIndex;
            }
        }
    }
    std::vector<size_t> erase;
    std::vector<int> padOfCircle(rw.circleArcs.size(), -1);
    for (size_t fi = 0; fi < neck.faceBeads.size(); ++fi) {
        FaceBead& fb = neck.faceBeads[fi];
        if (fb.face < 0) {
            continue;
        }
        if (fb.face >= numOldFaces) {
            erase.push_back(fi);
            continue;
        }
        int target = -1;
        int padCircle = -1;
        if (fb.face == in.splitFace) {
            if (exteriorIsSplit) {
                erase.push_back(fi);
                continue;
            }
            const Piece& o = pieces[static_cast<size_t>(outer)];
            if (o.circle) {
                padCircle = o.circleIndex;
            }
            else {
                target = o.pieceFace;
            }
        }
        else {
            target = result.faceMap[static_cast<size_t>(fb.face)];
            if (target < 0) {
                padCircle = CircleClosingOldFace(oldLoop, rw, fb.face);
            }
        }
        if (padCircle >= 0) {
            if (padOfCircle[static_cast<size_t>(padCircle)] >= 0) {
                neck.faceBeads[static_cast<size_t>(padOfCircle[static_cast<size_t>(padCircle)])]
                    .highlighted = neck.faceBeads[static_cast<size_t>(
                                       padOfCircle[static_cast<size_t>(padCircle)])]
                                       .highlighted ||
                                   fb.highlighted;
                erase.push_back(fi);
                continue;
            }
            const int seed = RingSeedForStrand(neck, rw.circleStrand[static_cast<size_t>(padCircle)]);
            const std::vector<int> ring = CircleRingThrough(neck, seed);
            if (ring.size() < 3) {
                return false;
            }
            if (fb.face == in.splitFace) {
                splitBeadKeptAsPad = true;
            }
            MakePad(fb, ring);
            padOfCircle[static_cast<size_t>(padCircle)] = static_cast<int>(fi);
            continue;
        }
        if (target < 0 || newExteriors.count(target) != 0) {
            erase.push_back(fi);
            continue;
        }
        RebindFaceBead(fb, target, newLoop, neck);
    }
    // Any circle without a pad gets a fresh one at its ring centroid.
    for (size_t ci = 0; ci < rw.circleArcs.size(); ++ci) {
        if (padOfCircle[ci] >= 0) {
            continue;
        }
        const int seed = RingSeedForStrand(neck, rw.circleStrand[ci]);
        const std::vector<int> ring = CircleRingThrough(neck, seed);
        if (ring.size() < 3) {
            return false;
        }
        int label = 1;
        for (const FaceBead& other : neck.faceBeads) {
            label = std::max(label, other.label + 1);
        }
        FaceBead pad;
        pad.label = label;
        if (splitHighlighted && !splitBeadKeptAsPad &&
            splitOuterCircle == static_cast<int>(ci)) {
            pad.highlighted = true;
        }
        MakePad(pad, ring);
        float cx = 0.0f;
        float cy = 0.0f;
        for (const int b : ring) {
            cx += neck.beads[static_cast<size_t>(b)].pos.x;
            cy += neck.beads[static_cast<size_t>(b)].pos.y;
        }
        pad.pos = {cx / static_cast<float>(ring.size()), cy / static_cast<float>(ring.size())};
        neck.faceBeads.push_back(std::move(pad));
    }
    std::sort(erase.begin(), erase.end());
    for (size_t i = erase.size(); i-- > 0;) {
        neck.faceBeads.erase(neck.faceBeads.begin() + static_cast<long>(erase[i]));
    }
    // Merge pair: two beads now share the merged face.
    if (mergedFace >= 0) {
        int keep = -1;
        int drop = -1;
        for (const FaceBead& fb : neck.faceBeads) {
            if (fb.face != mergedFace) {
                continue;
            }
            if (keep < 0) {
                keep = fb.label;
            }
            else if (fb.label < keep) {
                drop = keep;
                keep = fb.label;
            }
            else {
                drop = fb.label;
            }
        }
        if (drop >= 0) {
            result.mergeKeepLabel = keep;
            result.mergeDropLabel = drop;
        }
    }

    // ---- Components. ----
    for (size_t i = 0; i < oldComponents.size(); ++i) {
        if (static_cast<int>(i) == in.movedComponent) {
            std::vector<int> order;
            order.push_back(outer);
            if (inner >= 0) {
                order.push_back(inner);
            }
            else if (result.disconnected) {
                order.push_back(1 - outer);
            }
            for (size_t oi = 0; oi < order.size(); ++oi) {
                const Piece& piece = pieces[static_cast<size_t>(order[oi])];
                MultiloopComponent c;
                if (piece.circle) {
                    c.kind = MultiloopPart::Circle;
                    c.strand = piece.strand;
                    c.exteriorFace = -1;
                }
                else {
                    c.kind = MultiloopPart::Crossings;
                    c.halfEdges = piece.halfEdges;
                    c.exteriorFace = piece.newExterior;
                    c.strand = moved.strand;
                }
                if (oi == 1 && inner >= 0) {
                    c.parentFace = splitParentFace;
                    c.parentCircle = splitParentCircle;
                }
                else {
                    c.parentFace = (moved.parentFace >= 0 && moved.parentFace < numOldFaces)
                                       ? result.faceMap[static_cast<size_t>(moved.parentFace)]
                                       : -1;
                    c.parentCircle = moved.parentCircle;
                    if (moved.parentFace >= 0 && c.parentFace < 0) {
                        // Parent face was the split face: it lives on in the
                        // outer piece's part.
                        c.parentFace = splitParentFace;
                        c.parentCircle = splitParentCircle;
                    }
                }
                result.components.push_back(std::move(c));
            }
            continue;
        }
        MultiloopComponent c = oldComponents[i];
        if (c.kind == MultiloopPart::Crossings) {
            std::vector<int> mapped;
            for (const int h : c.halfEdges) {
                const int m = rw.MapHalfEdge(h);
                if (m == 0) {
                    return false;
                }
                mapped.push_back(m);
            }
            c.halfEdges = std::move(mapped);
            c.exteriorFace = (c.exteriorFace >= 0 && c.exteriorFace < numOldFaces)
                                 ? result.faceMap[static_cast<size_t>(c.exteriorFace)]
                                 : -1;
            if (c.exteriorFace < 0) {
                c.exteriorFace = FirstFaceOnComponent(newLoop, c.halfEdges);
            }
        }
        if (c.parentFace >= 0) {
            const int nf = (c.parentFace < numOldFaces)
                               ? result.faceMap[static_cast<size_t>(c.parentFace)]
                               : -1;
            if (nf >= 0) {
                c.parentFace = nf;
            }
            else if (c.parentFace == in.splitFace) {
                // Combinatorial default: stay with the outer piece. If this
                // hole sits inside the new inner piece, nest there instead
                // (concentric circles from successive disconnecting R2s).
                bool placed = false;
                if (inner >= 0 && !exteriorIsSplit) {
                    Vec2d sample;
                    if (SampleComponentPoint(c, neck, rw, &sample)) {
                        std::vector<Vec2d> ipoly =
                            (static_cast<size_t>(inner) < piecePoly.size())
                                ? piecePoly[static_cast<size_t>(inner)]
                                : std::vector<Vec2d>{};
                        if (ipoly.size() < 3) {
                            ipoly = PiecePolygon(pieces[static_cast<size_t>(inner)], newLoop,
                                                 neck, rw);
                        }
                        if (ipoly.size() >= 3 && PointInPolygon(sample, ipoly)) {
                            NestInPiece(&c, pieces[static_cast<size_t>(inner)]);
                            placed = true;
                        }
                    }
                }
                if (!placed) {
                    c.parentFace = splitParentFace;
                    c.parentCircle = splitParentCircle;
                }
            }
            else {
                // The parent face may have closed into the interior of a new
                // circle (e.g. an R2 that turns a Hopf shadow into two
                // rings): the component is now nested in that circle.
                const int circle = CircleClosingOldFace(oldLoop, rw, c.parentFace);
                c.parentFace = -1;
                if (circle >= 0 && circle < static_cast<int>(rw.circleStrand.size())) {
                    c.parentCircle = rw.circleStrand[static_cast<size_t>(circle)];
                }
                if (inner >= 0 && !exteriorIsSplit) {
                    const Piece& inn = pieces[static_cast<size_t>(inner)];
                    const Piece& outp = pieces[static_cast<size_t>(outer)];
                    Vec2d sample;
                    if (outp.circle && c.parentCircle == outp.strand &&
                        SampleComponentPoint(c, neck, rw, &sample)) {
                        std::vector<Vec2d> ipoly =
                            (static_cast<size_t>(inner) < piecePoly.size())
                                ? piecePoly[static_cast<size_t>(inner)]
                                : std::vector<Vec2d>{};
                        if (ipoly.size() < 3) {
                            ipoly = PiecePolygon(inn, newLoop, neck, rw);
                        }
                        if (ipoly.size() >= 3 && PointInPolygon(sample, ipoly)) {
                            NestInPiece(&c, inn);
                        }
                    }
                }
            }
        }
        result.components.push_back(std::move(c));
    }

    for (const int k : in.exteriorArcEdges) {
        const int m = rw.MapEdge(k);
        if (m > 0 && std::find(result.exteriorArcs.begin(), result.exteriorArcs.end(), m) ==
                         result.exteriorArcs.end()) {
            result.exteriorArcs.push_back(m);
        }
    }
    *out = std::move(result);
    return true;
}

void RebindFaceBead(FaceBead& fb, int face, const Multiloop& newLoop, const Necklace& neck)
{
    fb.face = face;
    fb.phiCycle.clear();
    fb.vertexBeads.clear();
    fb.boundaryBeads.clear();
    fb.boundaryDots.clear();
    fb.holeBeads.clear();
    fb.cachedExcessArea = 0.0;
    if (face < 0 || face >= static_cast<int>(newLoop.faces.size())) {
        return;
    }
    fb.phiCycle = newLoop.faces[static_cast<size_t>(face)];
    const std::vector<std::vector<int>> cycles = newLoop.sigma.Cycles();
    std::unordered_map<int, int> crossingOf;
    for (int c = 0; c < static_cast<int>(cycles.size()); ++c) {
        for (const int h : cycles[static_cast<size_t>(c)]) {
            crossingOf[h] = c;
        }
    }
    for (const int h : fb.phiCycle) {
        const auto it = crossingOf.find(h);
        if (it != crossingOf.end() && it->second < static_cast<int>(neck.crossingBead.size())) {
            fb.vertexBeads.push_back(neck.crossingBead[static_cast<size_t>(it->second)]);
        }
    }
}
