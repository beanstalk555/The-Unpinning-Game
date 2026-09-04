#include "r3_move.h"

#include <algorithm>
#include <cstdlib>
#include <exception>
#include <set>
#include <unordered_map>
#include <vector>

namespace
{
    int Mod3(int i)
    {
        return ((i % 3) + 3) % 3;
    }

    // Ordered beads along one edge, walked from a vertex bead in the
    // direction of the given half-edge: [vertex, edge beads..., vertex].
    std::vector<int> WalkSideChain(const Necklace& neck, int vertexBead, int halfEdge)
    {
        std::vector<int> chain;
        chain.push_back(vertexBead);
        int prev = vertexBead;
        int cur = neck.NeighborToward(vertexBead, halfEdge);
        int guard = 0;
        const int limit = static_cast<int>(neck.beads.size()) + 8;
        while (cur >= 0 && !neck.beads[static_cast<size_t>(cur)].isVertex) {
            chain.push_back(cur);
            const Bead& b = neck.beads[static_cast<size_t>(cur)];
            const int nxt = (b.neighbors[0] == prev) ? b.neighbors[1] : b.neighbors[0];
            prev = cur;
            cur = nxt;
            if (++guard > limit) {
                break;
            }
        }
        if (cur >= 0) {
            chain.push_back(cur);
        }
        return chain;
    }

    int OverlapCount(const std::vector<int>& a, const std::vector<int>& b)
    {
        int count = 0;
        for (const int x : a) {
            if (std::find(b.begin(), b.end(), x) != b.end()) {
                ++count;
            }
        }
        return count;
    }
}

bool AnalyzeR3Face(const Multiloop& loop, const Necklace* neck, int face, R3FaceInfo* out)
{
    if (out == nullptr || face < 0 || face >= static_cast<int>(loop.faces.size())) {
        return false;
    }
    const std::vector<int>& cycle = loop.faces[static_cast<size_t>(face)];
    if (cycle.size() != 3) {
        return false;
    }

    R3FaceInfo info;
    info.face = face;
    for (int i = 0; i < 3; ++i) {
        info.h[i] = cycle[static_cast<size_t>(i)];
    }

    const std::vector<std::vector<int>> crossings = loop.sigma.Cycles();
    std::unordered_map<int, int> crossingOf;
    for (int c = 0; c < static_cast<int>(crossings.size()); ++c) {
        for (const int h : crossings[static_cast<size_t>(c)]) {
            crossingOf[h] = c;
        }
    }

    for (int i = 0; i < 3; ++i) {
        const auto it = crossingOf.find(info.h[i]);
        if (it == crossingOf.end()) {
            return false;
        }
        info.crossing[i] = it->second;
    }
    if (info.crossing[0] == info.crossing[1] || info.crossing[1] == info.crossing[2] ||
        info.crossing[0] == info.crossing[2]) {
        return false;
    }

    std::set<int> labels;
    for (int i = 0; i < 3; ++i) {
        const int prev = info.h[Mod3(i - 1)];
        if (loop.sigma.Apply(info.h[i]) != -prev) {
            return false; // face is not traversed in phi order
        }
        info.a[i] = loop.sigma.Apply(-prev);
        info.b[i] = loop.sigma.Apply(info.a[i]);
        if (loop.sigma.Apply(info.b[i]) != info.h[i]) {
            return false;
        }
        labels.insert(info.h[i]);
        labels.insert(-prev);
        labels.insert(info.a[i]);
        labels.insert(info.b[i]);
    }
    if (labels.size() != 12) {
        return false;
    }

    if (neck != nullptr) {
        std::set<int> terminals;
        for (int i = 0; i < 3; ++i) {
            if (info.crossing[i] >= static_cast<int>(neck->crossingBead.size())) {
                return false;
            }
            const int vb = neck->crossingBead[static_cast<size_t>(info.crossing[i])];
            info.vertexBead[i] = vb;
            info.aTerminal[i] = neck->NeighborToward(vb, info.a[i]);
            info.bTerminal[i] = neck->NeighborToward(vb, info.b[i]);
            for (const int t : {info.aTerminal[i], info.bTerminal[i]}) {
                if (t < 0 || neck->beads[static_cast<size_t>(t)].isVertex ||
                    neck->beads[static_cast<size_t>(t)].removed) {
                    return false;
                }
                terminals.insert(t);
            }
        }
        if (terminals.size() != 6) {
            return false;
        }
    }

    *out = info;
    return true;
}

bool ComputeR3Sigma(const Multiloop& loop, const R3FaceInfo& info, Permutation* out)
{
    if (out == nullptr) {
        return false;
    }
    std::vector<std::vector<int>> cycles = loop.sigma.Cycles();
    for (int i = 0; i < 3; ++i) {
        const int j = Mod3(i + 1);
        std::vector<int>& c = cycles[static_cast<size_t>(info.crossing[j])];
        if (info.ccw) {
            c = {info.h[j], -info.h[i], info.b[i], info.a[j]};
        }
        else {
            c = {info.h[j], info.a[j], info.b[i], -info.h[i]};
        }
    }
    try {
        Permutation sigma = Permutation::FromCycles(cycles);
        if (!ValidateMultiloopSigma(sigma).empty()) {
            return false;
        }
        *out = sigma;
    }
    catch (const std::exception&) {
        return false;
    }
    return true;
}

std::unique_ptr<Multiloop> FinalizeR3(const Multiloop& oldLoop, Necklace& neck,
                                      const R3FaceInfo& info, const int meeting[3],
                                      const std::vector<int>& oldExteriorFaces,
                                      std::vector<int>* newExteriorFaces)
{
    // ---- Phase 1: validate and compute; nothing is mutated until every
    // check has passed, so a failure leaves the necklace untouched. ----
    Permutation newSigma;
    if (!ComputeR3Sigma(oldLoop, info, &newSigma)) {
        return nullptr;
    }

    std::vector<int> side[3];
    int meetIdx[3] = {-1, -1, -1};
    int towardPrev[3];
    int towardNext[3];
    for (int i = 0; i < 3; ++i) {
        side[i] = WalkSideChain(neck, info.vertexBead[i], info.h[i]);
        if (side[i].size() < 3 || side[i].back() != info.vertexBead[Mod3(i + 1)]) {
            return nullptr;
        }
        for (size_t t = 1; t + 1 < side[i].size(); ++t) {
            if (side[i][t] == meeting[i]) {
                meetIdx[i] = static_cast<int>(t);
                break;
            }
        }
        if (meetIdx[i] < 0) {
            return nullptr;
        }
        towardPrev[i] = neck.NeighborToward(info.vertexBead[i], -info.h[Mod3(i - 1)]);
        towardNext[i] = neck.NeighborToward(info.vertexBead[i], info.h[i]);
        if (towardPrev[i] < 0 || towardNext[i] < 0) {
            return nullptr;
        }
    }

    // The terminals must still carry their animation cuts (one -1 neighbor).
    int aCut[3];
    int bCut[3];
    for (int i = 0; i < 3; ++i) {
        const Bead& at = neck.beads[static_cast<size_t>(info.aTerminal[i])];
        const Bead& bt = neck.beads[static_cast<size_t>(info.bTerminal[i])];
        aCut[i] = (at.neighbors[0] < 0) ? 0 : ((at.neighbors[1] < 0) ? 1 : -1);
        bCut[i] = (bt.neighbors[0] < 0) ? 0 : ((bt.neighbors[1] < 0) ? 1 : -1);
        if (aCut[i] < 0 || bCut[i] < 0) {
            return nullptr;
        }
    }

    std::unique_ptr<Multiloop> newLoop;
    try {
        newLoop.reset(new Multiloop(newSigma));
    }
    catch (const std::exception&) {
        return nullptr;
    }
    if (newLoop->numEdges != oldLoop.numEdges ||
        newLoop->faces.size() != oldLoop.faces.size()) {
        return nullptr;
    }

    // R3 redistributes a few edges among strands (the triangle sides trade
    // components) and Multiloop assigns strand IDs by first-seen root, so the
    // default numbering often permutes colors. Remap each new strand onto the
    // old strand it overlaps most so colors stay stable for the bulk of each
    // component.
    {
        struct StrandCand
        {
            int overlap;
            int newStrand;
            int oldStrand;
        };
        std::vector<StrandCand> strandCands;
        for (int ns = 0; ns < newLoop->numStrands; ++ns) {
            for (int os = 0; os < oldLoop.numStrands; ++os) {
                int ov = 0;
                for (int k = 1; k <= oldLoop.numEdges; ++k) {
                    if (oldLoop.strandOfEdge[static_cast<size_t>(k)] == os &&
                        newLoop->strandOfEdge[static_cast<size_t>(k)] == ns) {
                        ++ov;
                    }
                }
                if (ov > 0) {
                    strandCands.push_back({ov, ns, os});
                }
            }
        }
        std::sort(strandCands.begin(), strandCands.end(),
                  [](const StrandCand& a, const StrandCand& b) {
                      if (a.overlap != b.overlap) {
                          return a.overlap > b.overlap;
                      }
                      if (a.newStrand != b.newStrand) {
                          return a.newStrand < b.newStrand;
                      }
                      return a.oldStrand < b.oldStrand;
                  });
        std::vector<int> newToOld(static_cast<size_t>(newLoop->numStrands), -1);
        std::vector<bool> oldTaken(static_cast<size_t>(oldLoop.numStrands), false);
        for (const StrandCand& c : strandCands) {
            if (newToOld[static_cast<size_t>(c.newStrand)] >= 0 ||
                oldTaken[static_cast<size_t>(c.oldStrand)]) {
                continue;
            }
            newToOld[static_cast<size_t>(c.newStrand)] = c.oldStrand;
            oldTaken[static_cast<size_t>(c.oldStrand)] = true;
        }
        // Any unmatched new strand (shouldn't happen when counts match) gets
        // the first free old id.
        for (int ns = 0; ns < newLoop->numStrands; ++ns) {
            if (newToOld[static_cast<size_t>(ns)] >= 0) {
                continue;
            }
            for (int os = 0; os < oldLoop.numStrands; ++os) {
                if (!oldTaken[static_cast<size_t>(os)]) {
                    newToOld[static_cast<size_t>(ns)] = os;
                    oldTaken[static_cast<size_t>(os)] = true;
                    break;
                }
            }
        }
        for (int k = 1; k <= newLoop->numEdges; ++k) {
            const int ns = newLoop->strandOfEdge[static_cast<size_t>(k)];
            if (ns >= 0 && ns < newLoop->numStrands &&
                newToOld[static_cast<size_t>(ns)] >= 0) {
                newLoop->strandOfEdge[static_cast<size_t>(k)] =
                    newToOld[static_cast<size_t>(ns)];
            }
        }
    }

    // Match old faces (face beads plus every component exterior) to new faces
    // by half-edge overlap. Faces away from the triangle match exactly; the
    // affected ones keep enough labels for a best-overlap assignment.
    std::vector<int> oldFacesNeeded;
    std::vector<int> faceBeadAssignIndex(neck.faceBeads.size(), -1);
    oldFacesNeeded.reserve(neck.faceBeads.size() + oldExteriorFaces.size());
    for (size_t fi = 0; fi < neck.faceBeads.size(); ++fi) {
        const FaceBead& fb = neck.faceBeads[fi];
        if (fb.face < 0) {
            continue;
        }
        faceBeadAssignIndex[fi] = static_cast<int>(oldFacesNeeded.size());
        oldFacesNeeded.push_back(fb.face);
    }
    const size_t exteriorStart = oldFacesNeeded.size();
    for (const int of : oldExteriorFaces) {
        if (of < 0) {
            continue;
        }
        oldFacesNeeded.push_back(of);
    }

    struct Candidate
    {
        int overlap;
        int oldIdx;
        int newFace;
    };
    std::vector<Candidate> candidates;
    for (int oi = 0; oi < static_cast<int>(oldFacesNeeded.size()); ++oi) {
        const int of = oldFacesNeeded[static_cast<size_t>(oi)];
        if (of < 0 || of >= static_cast<int>(oldLoop.faces.size())) {
            return nullptr;
        }
        const std::vector<int>& oldCycle = oldLoop.faces[static_cast<size_t>(of)];
        for (int nf = 0; nf < static_cast<int>(newLoop->faces.size()); ++nf) {
            const int ov = OverlapCount(oldCycle, newLoop->faces[static_cast<size_t>(nf)]);
            if (ov > 0) {
                candidates.push_back({ov, oi, nf});
            }
        }
    }
    std::sort(candidates.begin(), candidates.end(), [](const Candidate& x, const Candidate& y) {
        if (x.overlap != y.overlap) {
            return x.overlap > y.overlap;
        }
        if (x.oldIdx != y.oldIdx) {
            return x.oldIdx < y.oldIdx;
        }
        return x.newFace < y.newFace;
    });
    std::vector<int> assignment(oldFacesNeeded.size(), -1);
    std::vector<bool> newTaken(newLoop->faces.size(), false);
    for (const Candidate& c : candidates) {
        if (assignment[static_cast<size_t>(c.oldIdx)] >= 0 || newTaken[static_cast<size_t>(c.newFace)]) {
            continue;
        }
        assignment[static_cast<size_t>(c.oldIdx)] = c.newFace;
        newTaken[static_cast<size_t>(c.newFace)] = true;
    }
    for (const int a : assignment) {
        if (a < 0) {
            return nullptr;
        }
    }

    const std::vector<std::vector<int>> newCycles = newSigma.Cycles();
    std::unordered_map<int, int> newCrossingOf;
    for (int c = 0; c < static_cast<int>(newCycles.size()); ++c) {
        for (const int h : newCycles[static_cast<size_t>(c)]) {
            newCrossingOf[h] = c;
        }
    }

    // ---- Phase 2: mutate the necklace. ----
    // Old triangle vertices become plain edge beads on the new side through
    // them (new edge |h[i]| runs meeting[i-1] -> X_i -> meeting[i]).
    for (int i = 0; i < 3; ++i) {
        Bead& v = neck.beads[static_cast<size_t>(info.vertexBead[i])];
        v.isVertex = false;
        v.crossing = -1;
        v.sigmaCycle.clear();
        for (int s = 0; s < 4; ++s) {
            v.neighborOfHalfEdge[s] = -1;
        }
        v.edge = std::abs(info.h[i]);
        if (info.h[i] > 0) {
            v.neighbors[0] = towardPrev[i];
            v.neighbors[1] = towardNext[i];
        }
        else {
            v.neighbors[0] = towardNext[i];
            v.neighbors[1] = towardPrev[i];
        }
    }

    // Beads of old side i-1 strictly between its meeting bead and X_i migrate
    // to edge |h[i]|. (The old side-i beads between X_i and meeting[i] already
    // carry label |h[i]| with the correct neighbor orientation.)
    for (int i = 0; i < 3; ++i) {
        const int im1 = Mod3(i - 1);
        const std::vector<int>& ch = side[im1];
        for (size_t t = static_cast<size_t>(meetIdx[im1]) + 1; t + 1 < ch.size(); ++t) {
            Bead& b = neck.beads[static_cast<size_t>(ch[t])];
            b.edge = std::abs(info.h[i]);
            const int p = ch[t - 1];
            const int s = ch[t + 1];
            if (info.h[i] > 0) {
                b.neighbors[0] = p;
                b.neighbors[1] = s;
            }
            else {
                b.neighbors[0] = s;
                b.neighbors[1] = p;
            }
        }
    }

    // Meeting beads become the new crossings; reconnect the cut terminals.
    for (int i = 0; i < 3; ++i) {
        const int j = Mod3(i + 1);
        Bead& m = neck.beads[static_cast<size_t>(meeting[i])];
        const int prevBead = side[i][static_cast<size_t>(meetIdx[i]) - 1];
        const int nextBead = side[i][static_cast<size_t>(meetIdx[i]) + 1];
        m.isVertex = true;
        m.removed = false;
        m.edge = -1;
        m.neighbors[0] = -1;
        m.neighbors[1] = -1;
        if (info.ccw) {
            m.sigmaCycle = {info.h[j], -info.h[i], info.b[i], info.a[j]};
            m.neighborOfHalfEdge[0] = nextBead;
            m.neighborOfHalfEdge[1] = prevBead;
            m.neighborOfHalfEdge[2] = info.bTerminal[i];
            m.neighborOfHalfEdge[3] = info.aTerminal[j];
        }
        else {
            m.sigmaCycle = {info.h[j], info.a[j], info.b[i], -info.h[i]};
            m.neighborOfHalfEdge[0] = nextBead;
            m.neighborOfHalfEdge[1] = info.aTerminal[j];
            m.neighborOfHalfEdge[2] = info.bTerminal[i];
            m.neighborOfHalfEdge[3] = prevBead;
        }
        neck.beads[static_cast<size_t>(info.bTerminal[i])].neighbors[bCut[i]] = meeting[i];
        neck.beads[static_cast<size_t>(info.aTerminal[j])].neighbors[aCut[j]] = meeting[i];
    }

    // Crossing indices under the new sigma.
    neck.crossingBead.assign(newCycles.size(), -1);
    for (size_t bi = 0; bi < neck.beads.size(); ++bi) {
        Bead& b = neck.beads[bi];
        if (!b.isVertex || b.sigmaCycle.empty()) {
            continue;
        }
        const auto it = newCrossingOf.find(b.sigmaCycle[0]);
        if (it == newCrossingOf.end()) {
            continue; // cannot happen: all vertex cycles come from newSigma
        }
        b.crossing = it->second;
        neck.crossingBead[static_cast<size_t>(it->second)] = static_cast<int>(bi);
    }

    // Face beads onto the new faces. Circle beads (face < 0) stay as-is.
    for (size_t fi = 0; fi < neck.faceBeads.size(); ++fi) {
        FaceBead& fb = neck.faceBeads[fi];
        if (fb.face < 0 || faceBeadAssignIndex[fi] < 0) {
            continue;
        }
        const int nf = assignment[static_cast<size_t>(faceBeadAssignIndex[fi])];
        fb.face = nf;
        fb.phiCycle = newLoop->faces[static_cast<size_t>(nf)];
        fb.vertexBeads.clear();
        for (const int h : fb.phiCycle) {
            fb.vertexBeads.push_back(neck.crossingBead[static_cast<size_t>(newCrossingOf[h])]);
        }
        fb.boundaryBeads.clear();
        fb.boundaryDots.clear();
        fb.cachedExcessArea = 0.0;
    }

    if (newExteriorFaces != nullptr) {
        newExteriorFaces->clear();
        size_t k = 0;
        for (const int of : oldExteriorFaces) {
            if (of < 0) {
                newExteriorFaces->push_back(-1);
                continue;
            }
            newExteriorFaces->push_back(assignment[exteriorStart + k]);
            ++k;
        }
    }
    return newLoop;
}

std::unique_ptr<Multiloop> FinalizeR3(const Multiloop& oldLoop, Necklace& neck,
                                      const R3FaceInfo& info, const int meeting[3],
                                      int oldExteriorFace, int* newExteriorFace)
{
    std::vector<int> olds;
    if (oldExteriorFace >= 0) {
        olds.push_back(oldExteriorFace);
    }
    std::vector<int> news;
    std::unique_ptr<Multiloop> result =
        FinalizeR3(oldLoop, neck, info, meeting, olds, &news);
    if (newExteriorFace != nullptr) {
        *newExteriorFace = news.empty() ? -1 : news.front();
    }
    return result;
}
