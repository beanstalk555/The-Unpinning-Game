#include "r1_move.h"

#include <algorithm>
#include <cstdlib>
#include <exception>
#include <unordered_map>

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
}

bool AnalyzeR1Face(const Multiloop& loop, const Necklace* neck, int face, R1FaceInfo* out)
{
    if (out == nullptr || face < 0 || face >= static_cast<int>(loop.faces.size())) {
        return false;
    }
    const std::vector<int>& cycle = loop.faces[static_cast<size_t>(face)];
    if (cycle.size() != 1) {
        return false;
    }
    R1FaceInfo info;
    info.face = face;
    info.e = cycle[0];
    if (loop.sigma.Apply(info.e) != -info.e) {
        return false;
    }
    info.p = loop.sigma.Apply(-info.e);
    info.q = loop.sigma.Apply(info.p);
    if (loop.sigma.Apply(info.q) != info.e) {
        return false;
    }
    if (std::abs(info.p) == std::abs(info.e) || std::abs(info.q) == std::abs(info.e) ||
        info.p == info.q) {
        return false;
    }
    const std::vector<std::vector<int>> crossings = loop.sigma.Cycles();
    for (int c = 0; c < static_cast<int>(crossings.size()); ++c) {
        for (const int h : crossings[static_cast<size_t>(c)]) {
            if (h == info.e) {
                info.crossing = c;
            }
        }
    }
    if (info.crossing < 0) {
        return false;
    }
    if (neck != nullptr) {
        if (info.crossing >= static_cast<int>(neck->crossingBead.size())) {
            return false;
        }
        info.vertexBead = neck->crossingBead[static_cast<size_t>(info.crossing)];
        info.pTerminal = neck->NeighborToward(info.vertexBead, info.p);
        info.qTerminal = neck->NeighborToward(info.vertexBead, info.q);
        for (const int t : {info.pTerminal, info.qTerminal}) {
            if (t < 0 || t >= static_cast<int>(neck->beads.size()) ||
                neck->beads[static_cast<size_t>(t)].isVertex ||
                neck->beads[static_cast<size_t>(t)].removed) {
                return false;
            }
        }
        if (info.pTerminal == info.qTerminal) {
            return false;
        }
    }
    *out = info;
    return true;
}

std::unique_ptr<Multiloop> FinalizeR1(const Multiloop& oldLoop, Necklace& neck,
                                      const R1FaceInfo& info,
                                      const std::vector<MultiloopComponent>& oldComponents,
                                      MoveComponentResult* out)
{
    // ---- Phase 1: validate and compute. ----
    R1FaceInfo fresh;
    if (!AnalyzeR1Face(oldLoop, &neck, info.face, &fresh) || fresh.e != info.e) {
        return nullptr;
    }
    const int moved = ComponentContaining(oldComponents, fresh.e);
    if (moved < 0) {
        return nullptr;
    }
    RewriteResult rw;
    if (!DeleteCrossingsAndSplice(oldLoop, {fresh.crossing}, {std::abs(fresh.e)},
                                  {{fresh.p, fresh.q}}, &rw)) {
        return nullptr;
    }
    const std::vector<int> loopChain = WalkSide(neck, fresh.vertexBead, fresh.e);
    if (loopChain.size() < 3 || loopChain.back() != fresh.vertexBead) {
        return nullptr;
    }
    for (size_t t = 1; t + 1 < loopChain.size(); ++t) {
        if (loopChain[t] == fresh.pTerminal || loopChain[t] == fresh.qTerminal) {
            return nullptr;
        }
    }
    const Bead& pt = neck.beads[static_cast<size_t>(fresh.pTerminal)];
    const Bead& qt = neck.beads[static_cast<size_t>(fresh.qTerminal)];
    if ((pt.neighbors[0] != fresh.vertexBead && pt.neighbors[1] != fresh.vertexBead) ||
        (qt.neighbors[0] != fresh.vertexBead && qt.neighbors[1] != fresh.vertexBead)) {
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
    for (size_t t = 1; t + 1 < loopChain.size(); ++t) {
        Bead& b = neck.beads[static_cast<size_t>(loopChain[t])];
        b.removed = true;
        b.neighbors[0] = -1;
        b.neighbors[1] = -1;
    }
    Bead& x = neck.beads[static_cast<size_t>(fresh.vertexBead)];
    x.isVertex = false;
    x.crossing = -1;
    x.sigmaCycle.clear();
    for (int s = 0; s < 4; ++s) {
        x.neighborOfHalfEdge[s] = -1;
    }
    x.edge = std::abs(fresh.p);
    x.neighbors[0] = fresh.pTerminal;
    x.neighbors[1] = fresh.qTerminal;

    if (!ApplyRewriteToNecklace(neck, rw, *newLoop)) {
        neck = backup;
        return nullptr;
    }
    RemapStrandsThroughRewrite(oldLoop, rw, newLoop.get());

    MoveBookkeepingInput in;
    in.movedComponent = moved;
    in.exteriorArcEdges = {std::abs(fresh.p), std::abs(fresh.q)};
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
