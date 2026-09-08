#include "r2_move.h"

#include <algorithm>
#include <cstdlib>
#include <exception>
#include <set>

namespace
{
    std::vector<int> WalkSide(const Necklace& neck, int vertexBead, int halfEdge)
    {
        std::vector<int> chain;
        chain.push_back(vertexBead);
        int prev = vertexBead;
        int cur = neck.NeighborToward(vertexBead, halfEdge);
        int guard = 0;
        const int limit = static_cast<int>(neck.beads.size()) + 8;
        while (cur >= 0 && !neck.beads[static_cast<size_t>(cur)].isVertex) {
            if (neck.beads[static_cast<size_t>(cur)].removed) {
                break;
            }
            chain.push_back(cur);
            const Bead& b = neck.beads[static_cast<size_t>(cur)];
            const int nxt = (b.neighbors[0] == prev) ? b.neighbors[1] : b.neighbors[0];
            prev = cur;
            cur = nxt;
            if (++guard > limit) {
                break;
            }
        }
        if (cur >= 0 && neck.beads[static_cast<size_t>(cur)].isVertex) {
            chain.push_back(cur);
        }
        return chain;
    }

    int ComponentContaining(const std::vector<MultiloopComponent>& components, int halfEdge)
    {
        for (size_t i = 0; i < components.size(); ++i) {
            const MultiloopComponent& c = components[i];
            if (c.kind != MultiloopPart::Crossings) {
                continue;
            }
            if (std::find(c.halfEdges.begin(), c.halfEdges.end(), halfEdge) != c.halfEdges.end()) {
                return static_cast<int>(i);
            }
        }
        return -1;
    }

    void ReplaceOnce(Bead& b, int oldNeighbor, int newNeighbor)
    {
        for (int i = 0; i < 2; ++i) {
            if (b.neighbors[i] == oldNeighbor) {
                b.neighbors[i] = newNeighbor;
                return;
            }
        }
    }
}

bool AnalyzeR2Face(const Multiloop& loop, const Necklace* neck, int face, R2FaceInfo* out)
{
    if (out == nullptr || face < 0 || face >= static_cast<int>(loop.faces.size())) {
        return false;
    }
    const std::vector<int>& cycle = loop.faces[static_cast<size_t>(face)];
    if (cycle.size() != 2) {
        return false;
    }
    R2FaceInfo info;
    info.face = face;
    info.h0 = cycle[0];
    info.h1 = cycle[1];
    if (std::abs(info.h0) == std::abs(info.h1)) {
        return false;
    }
    const std::vector<std::vector<int>> crossings = loop.sigma.Cycles();
    for (int c = 0; c < static_cast<int>(crossings.size()); ++c) {
        for (const int h : crossings[static_cast<size_t>(c)]) {
            if (h == info.h0) {
                info.crossing[0] = c;
            }
            if (h == info.h1) {
                info.crossing[1] = c;
            }
        }
    }
    if (info.crossing[0] < 0 || info.crossing[1] < 0 || info.crossing[0] == info.crossing[1]) {
        return false;
    }
    if (loop.sigma.Apply(info.h0) != -info.h1 || loop.sigma.Apply(info.h1) != -info.h0) {
        return false;
    }
    info.a = loop.sigma.Apply(-info.h1);
    info.b = loop.sigma.Apply(info.a);
    if (loop.sigma.Apply(info.b) != info.h0) {
        return false;
    }
    info.ap = loop.sigma.Apply(-info.h0);
    info.bp = loop.sigma.Apply(info.ap);
    if (loop.sigma.Apply(info.bp) != info.h1) {
        return false;
    }
    {
        const std::set<int> labels = {info.h0, -info.h1, info.a, info.b,
                                      info.h1, -info.h0, info.ap, info.bp};
        if (labels.size() != 8) {
            return false;
        }
    }
    info.faceL = FaceIndexContainingHalfEdge(loop, info.a);
    info.faceR = FaceIndexContainingHalfEdge(loop, info.ap);
    if (info.faceL < 0 || info.faceR < 0) {
        return false;
    }
    info.disconnects = info.faceL == info.faceR;

    if (neck != nullptr) {
        for (int i = 0; i < 2; ++i) {
            if (info.crossing[i] >= static_cast<int>(neck->crossingBead.size())) {
                return false;
            }
            info.vertexBead[i] = neck->crossingBead[static_cast<size_t>(info.crossing[i])];
        }
        info.aTerminal = neck->NeighborToward(info.vertexBead[0], info.a);
        info.bTerminal = neck->NeighborToward(info.vertexBead[0], info.b);
        info.apTerminal = neck->NeighborToward(info.vertexBead[1], info.ap);
        info.bpTerminal = neck->NeighborToward(info.vertexBead[1], info.bp);
        std::set<int> terminals;
        for (const int t : {info.aTerminal, info.bTerminal, info.apTerminal, info.bpTerminal}) {
            if (t < 0 || t >= static_cast<int>(neck->beads.size()) ||
                neck->beads[static_cast<size_t>(t)].isVertex ||
                neck->beads[static_cast<size_t>(t)].removed) {
                return false;
            }
            terminals.insert(t);
        }
        if (terminals.size() != 4) {
            return false;
        }
    }
    *out = info;
    return true;
}

int R2OpposingHostHalfEdge(const Multiloop& loop, const R2FaceInfo& info, int exteriorFace)
{
    if (!info.disconnects || info.faceL < 0 || exteriorFace == info.faceL) {
        return 0;
    }
    if (exteriorFace < 0 || exteriorFace >= static_cast<int>(loop.faces.size())) {
        return 0;
    }
    bool outerIsS1 = false;
    bool found = false;
    for (const int h : loop.faces[static_cast<size_t>(exteriorFace)]) {
        const int k = std::abs(h);
        if (k == std::abs(info.a) || k == std::abs(info.bp)) {
            outerIsS1 = true;
            found = true;
            break;
        }
        if (k == std::abs(info.b) || k == std::abs(info.ap)) {
            outerIsS1 = false;
            found = true;
            break;
        }
    }
    if (!found) {
        return 0;
    }
    return outerIsS1 ? info.a : info.ap;
}

std::unique_ptr<Multiloop> FinalizeR2(const Multiloop& oldLoop, Necklace& neck,
                                      const R2FaceInfo& info,
                                      const std::vector<MultiloopComponent>& oldComponents,
                                      MoveComponentResult* out)
{
    // ---- Phase 1: validate and compute. ----
    R2FaceInfo fresh;
    if (!AnalyzeR2Face(oldLoop, &neck, info.face, &fresh) || fresh.h0 != info.h0 ||
        fresh.h1 != info.h1) {
        return nullptr;
    }
    const int moved = ComponentContaining(oldComponents, fresh.h0);
    if (moved < 0) {
        return nullptr;
    }
    RewriteResult rw;
    if (!DeleteCrossingsAndSplice(oldLoop, {fresh.crossing[0], fresh.crossing[1]},
                                  {std::abs(fresh.h0), std::abs(fresh.h1)},
                                  {{fresh.a, fresh.bp}, {fresh.b, fresh.ap}}, &rw)) {
        return nullptr;
    }
    const int X = fresh.vertexBead[0];
    const int Y = fresh.vertexBead[1];
    const std::vector<int> side0 = WalkSide(neck, X, fresh.h0);
    const std::vector<int> side1 = WalkSide(neck, Y, fresh.h1);
    if (side0.size() != 3 || side1.size() != 3 || side0.back() != Y || side1.back() != X) {
        return nullptr;
    }
    const int m0 = side0[1];
    const int m1 = side1[1];
    if (m0 == m1) {
        return nullptr;
    }
    for (const int t : {fresh.aTerminal, fresh.bTerminal, fresh.apTerminal, fresh.bpTerminal}) {
        if (t == m0 || t == m1) {
            return nullptr;
        }
    }
    const auto touches = [&](int bead, int target) {
        const Bead& b = neck.beads[static_cast<size_t>(bead)];
        return b.neighbors[0] == target || b.neighbors[1] == target;
    };
    if (!touches(fresh.aTerminal, X) || !touches(fresh.bTerminal, X) ||
        !touches(fresh.apTerminal, Y) || !touches(fresh.bpTerminal, Y)) {
        return nullptr;
    }
    std::unique_ptr<Multiloop> newLoop;
    try {
        newLoop.reset(new Multiloop(rw.sigma));
    }
    catch (const std::exception&) {
        return nullptr;
    }

    // ---- Phase 2: mutate, restoring the snapshot on failure. ----
    const Necklace backup = neck;
    const auto toEdgeBead = [&](int bead, int tempEdge, int n0, int n1) {
        Bead& v = neck.beads[static_cast<size_t>(bead)];
        v.isVertex = false;
        v.crossing = -1;
        v.sigmaCycle.clear();
        for (int s = 0; s < 4; ++s) {
            v.neighborOfHalfEdge[s] = -1;
        }
        v.edge = tempEdge;
        v.neighbors[0] = n0;
        v.neighbors[1] = n1;
    };
    toEdgeBead(X, std::abs(fresh.a), fresh.aTerminal, m1);
    toEdgeBead(Y, std::abs(fresh.ap), m0, fresh.apTerminal);
    {
        Bead& bm1 = neck.beads[static_cast<size_t>(m1)];
        bm1.neighbors[0] = X;
        bm1.neighbors[1] = fresh.bpTerminal;
        Bead& bm0 = neck.beads[static_cast<size_t>(m0)];
        bm0.neighbors[0] = fresh.bTerminal;
        bm0.neighbors[1] = Y;
    }
    ReplaceOnce(neck.beads[static_cast<size_t>(fresh.bTerminal)], X, m0);
    ReplaceOnce(neck.beads[static_cast<size_t>(fresh.bpTerminal)], Y, m1);

    if (!ApplyRewriteToNecklace(neck, rw, *newLoop)) {
        neck = backup;
        return nullptr;
    }
    RemapStrandsThroughRewrite(oldLoop, rw, newLoop.get());

    MoveBookkeepingInput in;
    in.movedComponent = moved;
    in.pieceHalfA = fresh.a;
    in.pieceHalfB = fresh.ap;
    if (fresh.disconnects) {
        in.splitFace = fresh.faceL;
    }
    else {
        in.mergeFaceA = fresh.faceL;
        in.mergeFaceB = fresh.faceR;
    }
    in.exteriorArcEdges = {std::abs(fresh.a), std::abs(fresh.b), std::abs(fresh.ap),
                           std::abs(fresh.bp)};
    MoveComponentResult result;
    if (!FinishMoveBookkeeping(oldLoop, *newLoop, rw, neck, oldComponents, in, &result)) {
        neck = backup;
        return nullptr;
    }
    if (out != nullptr) {
        *out = std::move(result);
    }
    return newLoop;
}
