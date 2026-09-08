#include "physics.h"

#include "drawing.h"
#include "geometry.h"
#include "move_rewrite.h"
#include "permutation.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <unordered_map>
#include <unordered_set>

namespace
{
    constexpr float kPiF = 3.14159265358979323846f;
    constexpr int kMinStrandEdgeBeads = 6;

    class ScopedTimer
    {
    public:
        explicit ScopedTimer(double* acc)
            : acc_(acc), t0_(std::chrono::high_resolution_clock::now())
        {
        }
        ~ScopedTimer()
        {
            if (acc_ == nullptr) {
                return;
            }
            const auto t1 = std::chrono::high_resolution_clock::now();
            *acc_ += std::chrono::duration<double, std::milli>(t1 - t0_).count();
        }

        private:
        double* acc_;
        std::chrono::high_resolution_clock::time_point t0_;
    };
}

Physics::Physics(const Multiloop& loop, Necklace& necklace, PhysicsParams& params,
                 int exteriorFace)
    : Physics(loop, necklace, params, std::unordered_set<int>{exteriorFace})
{
}

Physics::Physics(const Multiloop& loop, Necklace& necklace, PhysicsParams& params,
                 std::unordered_set<int> exteriorFaces)
    : loop_(&loop),
      neck_(&necklace),
      params_(&params),
      exteriorFace_(exteriorFaces.empty() ? -1 : *exteriorFaces.begin()),
      exteriorFaces_(std::move(exteriorFaces)),
      radius_(necklace.radius),
      diameter_(2.0f * necklace.radius)
{
    vel_.assign(neck_->beads.size(), Vector2{0.0f, 0.0f});
    fvel_.assign(neck_->faceBeads.size(), Vector2{0.0f, 0.0f});

    const std::vector<std::vector<int>> cycles = loop.sigma.Cycles();
    for (size_t c = 0; c < cycles.size(); ++c) {
        for (const int h : cycles[c]) {
            crossingOf_[h] = static_cast<int>(c);
        }
    }
    for (int f = 0; f < static_cast<int>(loop.faces.size()); ++f) {
        for (const int h : loop.faces[f]) {
            faceOfHalfEdge_[h] = f;
        }
    }
    structureDirty_ = true;
}

Physics::Physics(Necklace& necklace, PhysicsParams& params)
    : loop_(nullptr),
      neck_(&necklace),
      params_(&params),
      exteriorFace_(1),
      exteriorFaces_{1},
      radius_(necklace.radius),
      diameter_(2.0f * necklace.radius),
      unknot_(true)
{
    dummyLoop_.reset(new Multiloop(Permutation()));
    dummyLoop_->numEdges = 1;
    dummyLoop_->numStrands = 1;
    dummyLoop_->strandOfEdge = {-1, 0};
    dummyLoop_->faces = {{1}, {-1}};
    loop_ = dummyLoop_.get();
    faceOfHalfEdge_[1] = 0;
    faceOfHalfEdge_[-1] = 1;
    vel_.assign(neck_->beads.size(), Vector2{0.0f, 0.0f});
    fvel_.assign(neck_->faceBeads.size(), Vector2{0.0f, 0.0f});
    structureDirty_ = true;
}

void Physics::ResetDebugTimings()
{
    debugTimings_ = {};
    debugTiming_ = true;
}

void Physics::SetDebugFrameTotalMs(double ms)
{
    debugTimings_.frameTotalMs = ms;
    debugTiming_ = false;
}

Vector2 Physics::PosOf(int body) const
{
    if (body < BeadCount()) {
        const BeadPos& p = neck_->beads[body].pos;
        return Vector2{p.x, p.y};
    }
    const BeadPos& p = neck_->faceBeads[body - BeadCount()].pos;
    return Vector2{p.x, p.y};
}

void Physics::SetPosOf(int body, Vector2 p)
{
    if (body < BeadCount()) {
        neck_->beads[body].pos = {p.x, p.y};
    }
    else {
        neck_->faceBeads[body - BeadCount()].pos = {p.x, p.y};
    }
}

Vector2& Physics::VelOf(int body)
{
    if (body < BeadCount()) {
        return vel_[body];
    }
    return fvel_[body - BeadCount()];
}

std::vector<int> Physics::WalkEdge(int edge) const
{
    std::vector<int> chain;
    const auto itStart = crossingOf_.find(edge);
    if (itStart == crossingOf_.end()) {
        return chain;
    }
    const int crossing = itStart->second;
    if (crossing < 0 || crossing >= static_cast<int>(neck_->crossingBead.size())) {
        return chain;
    }
    const int vStart = neck_->crossingBead[static_cast<size_t>(crossing)];
    if (vStart < 0 || vStart >= BeadCount()) {
        return chain;
    }
    chain.push_back(vStart);

    int prev = vStart;
    int cur = neck_->NeighborToward(vStart, edge);
    int guard = 0;
    const int limit = BeadCount() + 8;
    while (cur >= 0 && !neck_->beads[static_cast<size_t>(cur)].isVertex) {
        if (neck_->beads[static_cast<size_t>(cur)].removed) {
            break; // broken link into a removed bead
        }
        chain.push_back(cur);
        const Bead& b = neck_->beads[static_cast<size_t>(cur)];
        const int nxt = (b.neighbors[0] == prev) ? b.neighbors[1] : b.neighbors[0];
        if (nxt == prev) {
            break; // both neighbors equal (degenerate loop bead)
        }
        prev = cur;
        cur = nxt;
        if (++guard > limit) {
            break; // safety against a malformed link
        }
    }
    if (cur >= 0 && neck_->beads[static_cast<size_t>(cur)].isVertex) {
        chain.push_back(cur); // tail vertex bead
    }
    return chain;
}

bool Physics::EdgeIsLoop(int edge) const
{
    const auto itPos = crossingOf_.find(edge);
    const auto itNeg = crossingOf_.find(-edge);
    return itPos != crossingOf_.end() && itNeg != crossingOf_.end() &&
           itPos->second == itNeg->second;
}

bool Physics::NodeReferences(int node, int target) const
{
    if (node < 0 || node >= BeadCount() || target < 0) {
        return false;
    }
    const Bead& b = neck_->beads[static_cast<size_t>(node)];
    if (b.removed) {
        return false;
    }
    if (b.isVertex) {
        for (int i = 0; i < 4; ++i) {
            if (b.neighborOfHalfEdge[i] == target) {
                return true;
            }
        }
        return false;
    }
    return b.neighbors[0] == target || b.neighbors[1] == target;
}

int Physics::CountReferences(int node, int target) const
{
    if (node < 0 || node >= BeadCount()) {
        return 0;
    }
    const Bead& b = neck_->beads[static_cast<size_t>(node)];
    int count = 0;
    if (b.isVertex) {
        for (int i = 0; i < 4; ++i) {
            if (b.neighborOfHalfEdge[i] == target) {
                ++count;
            }
        }
    }
    else {
        for (int i = 0; i < 2; ++i) {
            if (b.neighbors[i] == target) {
                ++count;
            }
        }
    }
    return count;
}

void Physics::RebuildSpringsFromChains()
{
    const ScopedTimer timer(debugTiming_ ? &debugTimings_.rebuildChainsMs : nullptr);
    const int n = loop_->numEdges;
    edgeSprings_.assign(n + 1, {});
    for (int k = 1; k <= n; ++k) {
        const std::vector<int>& chain = chains_[k];
        std::vector<std::pair<int, int>>& es = edgeSprings_[static_cast<size_t>(k)];
        for (size_t i = 0; i + 1 < chain.size(); ++i) {
            const int a = chain[i];
            const int b = chain[i + 1];
            if (a == b) {
                continue;
            }
            if (a < BeadCount() && neck_->beads[static_cast<size_t>(a)].removed) {
                continue;
            }
            if (b < BeadCount() && neck_->beads[static_cast<size_t>(b)].removed) {
                continue;
            }
            es.push_back({a, b});
        }
    }
    FlattenSpringsFromEdges();
}

void Physics::RefreshEdgeSprings(int edge)
{
    if (edge < 1 || edge > loop_->numEdges) {
        return;
    }
    // Grow without wiping existing per-edge spring lists (assign would clear them).
    if (static_cast<int>(edgeSprings_.size()) < loop_->numEdges + 1) {
        edgeSprings_.resize(static_cast<size_t>(loop_->numEdges + 1));
    }
    const std::vector<int>& chain = chains_[edge];
    std::vector<std::pair<int, int>>& es = edgeSprings_[static_cast<size_t>(edge)];
    es.clear();
    for (size_t i = 0; i + 1 < chain.size(); ++i) {
        const int a = chain[i];
        const int b = chain[i + 1];
        if (a == b) {
            continue;
        }
        if (a < BeadCount() && neck_->beads[static_cast<size_t>(a)].removed) {
            continue;
        }
        if (b < BeadCount() && neck_->beads[static_cast<size_t>(b)].removed) {
            continue;
        }
        es.push_back({a, b});
    }
    FlattenSpringsFromEdges();
}

void Physics::FlattenSpringsFromEdges()
{
    springs_.clear();
    for (int k = 1; k <= loop_->numEdges; ++k) {
        if (static_cast<size_t>(k) >= edgeSprings_.size()) {
            continue;
        }
        const std::vector<std::pair<int, int>>& es = edgeSprings_[static_cast<size_t>(k)];
        springs_.insert(springs_.end(), es.begin(), es.end());
    }
    if (!unknot_) {
        AppendCircleSprings();
    }
}

void Physics::AppendCircleSprings()
{
    const auto appendRing = [&](const std::vector<int>& b) {
        if (b.size() < 2) {
            return;
        }
        for (size_t i = 0; i < b.size(); ++i) {
            const int a = b[i];
            const int c = b[(i + 1) % b.size()];
            if (a < 0 || c < 0 || a >= BeadCount() || c >= BeadCount()) {
                continue;
            }
            if (neck_->beads[static_cast<size_t>(a)].removed ||
                neck_->beads[static_cast<size_t>(c)].removed) {
                continue;
            }
            springs_.push_back({a, c});
        }
    };
    for (FaceBead& fb : neck_->faceBeads) {
        if (fb.face >= 0 || fb.boundaryBeads.size() < 3) {
            continue;
        }
        appendRing(fb.boundaryBeads);
        fb.RefreshBoundaryDots(*neck_);
        fb.cachedExcessArea = fb.ExcessArea(*neck_, radius_);
    }
    if (r0Active_) {
        appendRing(r0Pad_.boundaryBeads);
    }
}

bool Physics::IsExterior(int face) const
{
    return exteriorFaces_.count(face) != 0;
}

int Physics::ExteriorSlot(int face) const
{
    for (size_t i = 0; i < exteriorFaceOrder_.size(); ++i) {
        if (exteriorFaceOrder_[i] == face) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

void Physics::RebuildFaceGeometry()
{
    const ScopedTimer timer(debugTiming_ ? &debugTimings_.rebuildFaceGeometryMs : nullptr);
    sharpDotThreshold_ = std::cos(kPiF / 3.0f + params_->bendAngleEpsilon);

    neck_->sharpCornerEdges.clear();
    std::unordered_set<int> sharpSet;
    exteriorBoundaries_.clear();
    exteriorFaceOrder_.clear();
    for (int f = 0; f < static_cast<int>(loop_->faces.size()); ++f) {
        if (!IsExterior(f)) {
            continue;
        }
        std::vector<float> exteriorDots;
        std::vector<int> boundary;
        WalkFaceBoundaryFromPhiCycle(loop_->faces[f], *neck_, crossingOf_, boundary,
                                     exteriorDots, &sharpSet, sharpDotThreshold_);
        exteriorFaceOrder_.push_back(f);
        exteriorBoundaries_.push_back(std::move(boundary));
    }
    // Every face bead walks its own boundary (a merged face may briefly hold
    // two beads).
    for (FaceBead& fb : neck_->faceBeads) {
        if (fb.face < 0 || IsExterior(fb.face) ||
            fb.face >= static_cast<int>(loop_->faces.size())) {
            continue;
        }
        WalkFaceBoundaryFromPhiCycle(loop_->faces[static_cast<size_t>(fb.face)], *neck_,
                                     crossingOf_, fb.boundaryBeads, fb.boundaryDots, &sharpSet,
                                     sharpDotThreshold_);
    }
    ResolveHoles();
    for (FaceBead& fb : neck_->faceBeads) {
        if (fb.face >= 0) {
            if (IsExterior(fb.face)) {
                continue;
            }
            fb.cachedExcessArea = fb.ExcessArea(*neck_, radius_) - HoleArea(fb);
            continue;
        }
        if (fb.boundaryBeads.size() < 3) {
            continue;
        }
        fb.RefreshBoundaryDots(*neck_);
        fb.cachedExcessArea = fb.ExcessArea(*neck_, radius_) - HoleArea(fb);
    }
    neck_->sharpCornerEdges.assign(sharpSet.begin(), sharpSet.end());
}

double Physics::HoleArea(const FaceBead& fb) const
{
    double area = 0.0;
    for (const std::vector<int>& hole : fb.holeBeads) {
        std::vector<Vec2d> poly;
        poly.reserve(hole.size());
        for (const int bead : hole) {
            if (bead < 0 || bead >= BeadCount() ||
                neck_->beads[static_cast<size_t>(bead)].removed) {
                continue;
            }
            const Vector2 p = PosOf(bead);
            poly.push_back(Vec2d{static_cast<double>(p.x), static_cast<double>(p.y)});
        }
        area += ShoelaceArea(poly);
    }
    return area;
}

int Physics::PadStrand(const FaceBead& fb) const
{
    for (const int bead : fb.boundaryBeads) {
        if (bead >= 0 && bead < BeadCount() && !neck_->beads[static_cast<size_t>(bead)].removed) {
            return neck_->beads[static_cast<size_t>(bead)].strand;
        }
    }
    return -1;
}

std::vector<int> Physics::RingBeadsOfCircleStrand(int strand) const
{
    for (const FaceBead& fb : neck_->faceBeads) {
        if (IsCircleFaceBead(fb) && PadStrand(fb) == strand) {
            return fb.boundaryBeads;
        }
    }
    for (int i = 0; i < BeadCount(); ++i) {
        const Bead& b = neck_->beads[static_cast<size_t>(i)];
        if (!b.isVertex && !b.removed && b.edge == 0 && b.strand == strand) {
            return CircleRingThrough(*neck_, i);
        }
    }
    return {};
}

void Physics::ResolveHoles()
{
    for (FaceBead& fb : neck_->faceBeads) {
        fb.holeBeads.clear();
    }
    if (nesting_.empty()) {
        return;
    }
    for (const NestingHole& hole : nesting_) {
        std::vector<int> ring;
        if (hole.holeFace >= 0) {
            const int slot = ExteriorSlot(hole.holeFace);
            if (slot < 0) {
                continue;
            }
            ring = exteriorBoundaries_[static_cast<size_t>(slot)];
        }
        else if (hole.holeCircle >= 0) {
            ring = RingBeadsOfCircleStrand(hole.holeCircle);
        }
        ring.erase(std::remove_if(ring.begin(), ring.end(),
                                  [&](int bead) {
                                      return bead < 0 || bead >= BeadCount() ||
                                             neck_->beads[static_cast<size_t>(bead)].removed;
                                  }),
                   ring.end());
        if (ring.size() < 3) {
            continue;
        }
        for (FaceBead& fb : neck_->faceBeads) {
            const bool byFace = hole.parentFace >= 0 && fb.face == hole.parentFace;
            const bool byCircle = hole.parentCircle >= 0 && IsCircleFaceBead(fb) &&
                                  PadStrand(fb) == hole.parentCircle;
            if (byFace || byCircle) {
                fb.holeBeads.push_back(ring);
            }
        }
    }
}

void Physics::SetNesting(const std::vector<MultiloopComponent>& components)
{
    nesting_.clear();
    for (const MultiloopComponent& c : components) {
        if (!c.Nested()) {
            continue;
        }
        NestingHole hole;
        hole.parentFace = c.parentFace;
        hole.parentCircle = c.parentCircle;
        if (c.kind == MultiloopPart::Crossings) {
            hole.holeFace = c.exteriorFace;
        }
        else {
            hole.holeCircle = c.strand;
        }
        if (hole.holeFace < 0 && hole.holeCircle < 0) {
            continue;
        }
        nesting_.push_back(hole);
    }
    if (nesting_.empty()) {
        for (FaceBead& fb : neck_->faceBeads) {
            fb.holeBeads.clear();
        }
    }
    structureDirty_ = true;
}

bool Physics::FaceHasHoles(int face) const
{
    if (face < 0) {
        return false;
    }
    for (const NestingHole& hole : nesting_) {
        if (hole.parentFace == face) {
            return true;
        }
    }
    return false;
}

bool Physics::CircleHasHoles(int strand) const
{
    if (strand < 0) {
        return false;
    }
    for (const NestingHole& hole : nesting_) {
        if (hole.parentCircle == strand) {
            return true;
        }
    }
    return false;
}

void Physics::RebuildFaceBodiesFromCache()
{
    const ScopedTimer timer(debugTiming_ ? &debugTimings_.rebuildFaceBodiesMs : nullptr);
    faceBodies_.clear();
    faceBodies_.reserve(neck_->faceBeads.size() + 1);
    for (size_t fi = 0; fi < neck_->faceBeads.size(); ++fi) {
        const FaceBead& fb = neck_->faceBeads[fi];
        std::vector<int> bodies;
        bodies.push_back(BeadCount() + static_cast<int>(fi));
        bodies.insert(bodies.end(), fb.boundaryBeads.begin(), fb.boundaryBeads.end());
        // Nested components collide with this region's bead and boundary.
        for (const std::vector<int>& hole : fb.holeBeads) {
            bodies.insert(bodies.end(), hole.begin(), hole.end());
        }
        faceBodies_.push_back(std::move(bodies));
    }
    for (const std::vector<int>& boundary : exteriorBoundaries_) {
        if (!boundary.empty()) {
            faceBodies_.push_back(boundary);
        }
    }
    if (r0Active_ && !r0Pad_.boundaryBeads.empty()) {
        faceBodies_.push_back(r0Pad_.boundaryBeads);
    }
    PruneRemovedFromFaceBodies();
    RebuildHoneycombSpokes();
}

void Physics::PruneRemovedFromFaceBodies()
{
    for (std::vector<int>& group : faceBodies_) {
        group.erase(std::remove_if(group.begin(), group.end(),
                                   [&](int body) {
                                       return body < BeadCount() &&
                                              neck_->beads[static_cast<size_t>(body)].removed;
                                   }),
                    group.end());
    }
}

void Physics::RebuildHoneycombSpokes()
{
    honeycombFaceFi_.clear();
    honeycombSpokes_.clear();
    const int n = BeadCount();
    for (size_t fi = 0; fi < neck_->faceBeads.size(); ++fi) {
        const FaceBead& fb = neck_->faceBeads[fi];
        if (fb.phiCycle.size() > 3 || fb.boundaryBeads.size() != 6 || !fb.holeBeads.empty()) {
            continue;
        }
        if (mergeKeepLabel_ >= 0 && (fb.label == mergeKeepLabel_ || fb.label == mergeDropLabel_)) {
            continue;
        }
        int spokes[6];
        bool ok = true;
        for (int i = 0; i < 6; ++i) {
            const int b = fb.boundaryBeads[static_cast<size_t>(i)];
            if (b < 0 || b >= n || neck_->beads[static_cast<size_t>(b)].removed) {
                ok = false;
                break;
            }
            spokes[i] = b;
        }
        if (!ok) {
            continue;
        }
        honeycombFaceFi_.push_back(static_cast<int>(fi));
        honeycombSpokes_.insert(honeycombSpokes_.end(), spokes, spokes + 6);
    }
}

bool Physics::FaceHasHoneycombSpokes(int fi) const
{
    for (const int x : honeycombFaceFi_) {
        if (x == fi) {
            return true;
        }
    }
    return false;
}

int Physics::HoneycombSlotForFace(int face) const
{
    for (size_t h = 0; h < honeycombFaceFi_.size(); ++h) {
        const int fi = honeycombFaceFi_[h];
        if (fi >= 0 && fi < static_cast<int>(neck_->faceBeads.size()) &&
            neck_->faceBeads[static_cast<size_t>(fi)].face == face) {
            return static_cast<int>(h);
        }
    }
    return -1;
}

void Physics::AppendIncidentFaces(int edge, std::unordered_set<int>& out) const
{
    const auto itPos = faceOfHalfEdge_.find(edge);
    if (itPos != faceOfHalfEdge_.end() && !IsExterior(itPos->second)) {
        out.insert(itPos->second);
    }
    const auto itNeg = faceOfHalfEdge_.find(-edge);
    if (itNeg != faceOfHalfEdge_.end() && !IsExterior(itNeg->second)) {
        out.insert(itNeg->second);
    }
}

void Physics::WalkFaceBoundaryNoDots(int face)
{
    if (IsExterior(face)) {
        const int slot = ExteriorSlot(face);
        if (slot < 0) {
            return;
        }
        std::vector<float> unusedDots;
        WalkFaceBoundaryFromPhiCycle(loop_->faces[face], *neck_, crossingOf_,
                                     exteriorBoundaries_[static_cast<size_t>(slot)],
                                     unusedDots, nullptr, 0.0f, false);
        return;
    }
    for (FaceBead& fb : neck_->faceBeads) {
        if (fb.face != face) {
            continue;
        }
        WalkFaceBoundaryFromPhiCycle(loop_->faces[face], *neck_, crossingOf_, fb.boundaryBeads,
                                     fb.boundaryDots, nullptr, 0.0f, false);
    }
}

void Physics::RefreshFaceBody(int face)
{
    if (!nesting_.empty()) {
        // Hole rings may reference this face's exterior walk; rebuild all.
        ResolveHoles();
        RebuildFaceBodiesFromCache();
        return;
    }
    if (IsExterior(face)) {
        const int slot = ExteriorSlot(face);
        if (slot < 0) {
            return;
        }
        std::vector<int> bodies = exteriorBoundaries_[static_cast<size_t>(slot)];
        bodies.erase(std::remove_if(bodies.begin(), bodies.end(),
                                    [&](int body) {
                                        return body < BeadCount() &&
                                               neck_->beads[static_cast<size_t>(body)].removed;
                                    }),
                     bodies.end());
        size_t index = neck_->faceBeads.size();
        for (int i = 0; i < slot; ++i) {
            if (i < static_cast<int>(exteriorBoundaries_.size()) &&
                !exteriorBoundaries_[static_cast<size_t>(i)].empty()) {
                ++index;
            }
        }
        if (!bodies.empty()) {
            if (faceBodies_.size() > index) {
                faceBodies_[index] = std::move(bodies);
            }
            else {
                faceBodies_.push_back(std::move(bodies));
            }
        }
        RebuildHoneycombSpokes();
        return;
    }
    for (size_t fi = 0; fi < neck_->faceBeads.size(); ++fi) {
        if (neck_->faceBeads[fi].face != face) {
            continue;
        }
        const FaceBead& fb = neck_->faceBeads[fi];
        std::vector<int> bodies;
        bodies.push_back(BeadCount() + static_cast<int>(fi));
        bodies.insert(bodies.end(), fb.boundaryBeads.begin(), fb.boundaryBeads.end());
        bodies.erase(std::remove_if(bodies.begin(), bodies.end(),
                                    [&](int body) {
                                        return body < BeadCount() &&
                                               neck_->beads[static_cast<size_t>(body)].removed;
                                    }),
                     bodies.end());
        if (fi < faceBodies_.size()) {
            faceBodies_[fi] = std::move(bodies);
        }
        else {
            faceBodies_.push_back(std::move(bodies));
        }
    }
    RebuildHoneycombSpokes();
}

void Physics::RefreshFaceBodiesForEdge(int edge)
{
    std::unordered_set<int> faces;
    AppendIncidentFaces(edge, faces);
    // Inserts/removes on an exterior-bounding edge must re-walk that
    // exterior cycle. AppendIncidentFaces skips exterior faces (honeycomb),
    // so without this the new beads never join the exterior collision group.
    const auto itPos = faceOfHalfEdge_.find(edge);
    if (itPos != faceOfHalfEdge_.end() && IsExterior(itPos->second)) {
        faces.insert(itPos->second);
    }
    const auto itNeg = faceOfHalfEdge_.find(-edge);
    if (itNeg != faceOfHalfEdge_.end() && IsExterior(itNeg->second)) {
        faces.insert(itNeg->second);
    }
    for (const int f : faces) {
        WalkFaceBoundaryNoDots(f);
        RefreshFaceBody(f);
    }
}

void Physics::RefreshEdgeChain(int edge)
{
    if (edge < 1 || edge > loop_->numEdges) {
        return;
    }
    if (static_cast<int>(chains_.size()) < loop_->numEdges + 1) {
        chains_.resize(static_cast<size_t>(loop_->numEdges + 1));
    }
    chains_[edge] = WalkEdge(edge);
}

void Physics::RebuildStructures()
{
    const ScopedTimer timer(debugTiming_ ? &debugTimings_.rebuildStructuresMs : nullptr);
    if (unknot_) {
        RebuildUnknotStructures();
        return;
    }
    const int n = loop_->numEdges;
    chains_.assign(n + 1, {});
    for (int k = 1; k <= n; ++k) {
        chains_[k] = WalkEdge(k);
    }
    RebuildSpringsFromChains();
    RebuildFaceGeometry();
    RebuildFaceBodiesFromCache();
    structureDirty_ = false;
}

void Physics::RebuildUnknotStructures()
{
    std::vector<char> seen(static_cast<size_t>(std::max(0, BeadCount())), 0);
    std::vector<std::vector<int>> rings;
    for (int i = 0; i < BeadCount(); ++i) {
        const Bead& startBead = neck_->beads[static_cast<size_t>(i)];
        if (startBead.removed || startBead.isVertex || seen[static_cast<size_t>(i)]) {
            continue;
        }
        std::vector<int> ordered;
        int prev = -1;
        int cur = i;
        const int limit = BeadCount() + 2;
        for (int guard = 0; guard < limit; ++guard) {
            if (cur < 0 || cur >= BeadCount() || seen[static_cast<size_t>(cur)]) {
                break;
            }
            seen[static_cast<size_t>(cur)] = 1;
            ordered.push_back(cur);
            const Bead& b = neck_->beads[static_cast<size_t>(cur)];
            int nxt = -1;
            if (prev < 0) {
                nxt = (b.neighbors[1] >= 0) ? b.neighbors[1] : b.neighbors[0];
            }
            else {
                nxt = (b.neighbors[0] == prev) ? b.neighbors[1] : b.neighbors[0];
            }
            prev = cur;
            cur = nxt;
            if (cur == i) {
                break;
            }
        }
        if (!ordered.empty()) {
            ordered.push_back(ordered.front());
            rings.push_back(std::move(ordered));
        }
    }

    const int nRings = std::max(1, static_cast<int>(rings.size()));
    if (dummyLoop_ != nullptr) {
        dummyLoop_->numEdges = nRings;
        dummyLoop_->strandOfEdge.assign(static_cast<size_t>(nRings + 1), 0);
        dummyLoop_->strandOfEdge[0] = -1;
    }
    chains_.assign(static_cast<size_t>(nRings + 1), {});
    for (int r = 0; r < static_cast<int>(rings.size()); ++r) {
        chains_[static_cast<size_t>(r + 1)] = rings[static_cast<size_t>(r)];
    }
    RebuildSpringsFromChains();
    std::vector<char> ringUsed(rings.size(), 0);
    const auto ringStrand = [&](const std::vector<int>& ring) {
        for (const int bead : ring) {
            if (bead < 0 || bead >= BeadCount()) {
                continue;
            }
            const Bead& b = neck_->beads[static_cast<size_t>(bead)];
            if (!b.removed && !b.isVertex) {
                return b.strand;
            }
        }
        return -1;
    };
    for (size_t fi = 0; fi < neck_->faceBeads.size(); ++fi) {
        FaceBead& fb = neck_->faceBeads[fi];
        const int want = PadStrand(fb);
        int best = -1;
        int bestScore = -1;
        for (size_t r = 0; r < rings.size(); ++r) {
            if (ringUsed[r] != 0 || rings[r].size() < 2) {
                continue;
            }
            int overlap = 0;
            for (const int bead : fb.boundaryBeads) {
                if (std::find(rings[r].begin(), rings[r].end(), bead) != rings[r].end()) {
                    ++overlap;
                }
            }
            int score = overlap;
            const int rs = ringStrand(rings[r]);
            if (want >= 0 && rs == want) {
                score += 1000;
            }
            if (score > bestScore) {
                bestScore = score;
                best = static_cast<int>(r);
            }
        }
        if (best < 0) {
            for (size_t r = 0; r < rings.size(); ++r) {
                if (ringUsed[r] == 0 && rings[r].size() >= 2) {
                    best = static_cast<int>(r);
                    break;
                }
            }
        }
        if (best >= 0) {
            ringUsed[static_cast<size_t>(best)] = 1;
            const std::vector<int>& ring = rings[static_cast<size_t>(best)];
            fb.boundaryBeads.assign(ring.begin(), ring.end() - 1);
        }
        fb.RefreshBoundaryDots(*neck_);
        fb.cachedExcessArea = fb.ExcessArea(*neck_, radius_);
    }
    exteriorBoundaries_.clear();
    exteriorFaceOrder_.clear();
    ResolveHoles();
    RebuildFaceBodiesFromCache();
    structureDirty_ = false;
}

void Physics::ApplyFaceCollisions(std::vector<Vector2>& force, int total,
                                 const std::vector<char>& skipFaceBody)
{
    const ScopedTimer timer(debugTiming_ ? &debugTimings_.collisionMs : nullptr);
    collisionGrid_.SetCellSize(params_->collisionBoundsScale * diameter_);
    std::unordered_set<long long> done;
    const long long stride = total;
    const int beads = BeadCount();

    const auto applyPair = [&](int a, int b) {
        const long long lo = a;
        const long long hi = b;
        if (!done.insert(lo * stride + hi).second) {
            return;
        }
        const Vector2 pa = PosOf(a);
        const Vector2 pb = PosOf(b);
        const Vector2 d = {pb.x - pa.x, pb.y - pa.y};
        const float len = std::sqrt(d.x * d.x + d.y * d.y);
        if (len < 1e-6f || len >= diameter_) {
            return;
        }
        const float mag = params_->collisionStiffness * (diameter_ - len);
        const Vector2 dir = {d.x / len, d.y / len};
        force[a].x -= mag * dir.x;
        force[a].y -= mag * dir.y;
        force[b].x += mag * dir.x;
        force[b].y += mag * dir.y;
    };

    if (globalCollisions_) {
        for (int a = 0; a < beads; ++a) {
            if (neck_->beads[static_cast<size_t>(a)].removed) {
                continue;
            }
            for (int b = a + 1; b < beads; ++b) {
                if (neck_->beads[static_cast<size_t>(b)].removed) {
                    continue;
                }
                applyPair(a, b);
            }
        }
    }

    for (const std::vector<int>& group : faceBodies_) {
        collisionGrid_.Clear();
        collisionGrid_.Reserve(group.size());
        for (const int body : group) {
            if (body >= beads) {
                const int fi = body - beads;
                if (fi >= 0 && fi < static_cast<int>(skipFaceBody.size()) &&
                    skipFaceBody[static_cast<size_t>(fi)]) {
                    continue;
                }
            }
            const Vector2 p = PosOf(body);
            collisionGrid_.Insert(body, p.x, p.y);
        }
        collisionGrid_.ForEachCandidatePair([&](int a, int b) {
            applyPair(a, b);
        });
    }

    // Cross-component: one pass over every component's exterior boundary.
    // Interior faces stay in their per-face groups above.
    std::unordered_set<int> exteriorBodies;
    for (const std::vector<int>& boundary : exteriorBoundaries_) {
        for (const int bead : boundary) {
            if (bead >= 0 && bead < beads && !neck_->beads[static_cast<size_t>(bead)].removed) {
                exteriorBodies.insert(bead);
            }
        }
    }
    for (const FaceBead& fb : neck_->faceBeads) {
        if (!IsCircleFaceBead(fb)) {
            continue;
        }
        for (const int bead : fb.boundaryBeads) {
            if (bead >= 0 && bead < beads && !neck_->beads[static_cast<size_t>(bead)].removed) {
                exteriorBodies.insert(bead);
            }
        }
    }
    if (exteriorBodies.size() >= 2) {
        collisionGrid_.Clear();
        collisionGrid_.Reserve(exteriorBodies.size());
        for (const int bead : exteriorBodies) {
            const Vector2 p = PosOf(bead);
            collisionGrid_.Insert(bead, p.x, p.y);
        }
        collisionGrid_.ForEachCandidatePair([&](int a, int b) {
            applyPair(a, b);
        });
    }

    // R2 strip/pad: the bigon pad is gone, so opposing side interiors are not
    // in one face group. Collide them with each other (vertices stay out).
    if (r2Active_ && shortenSub_ != ShortenSub::Spoke) {
        RebuildMoveRepelSides();
        std::vector<int> sides;
        for (const std::vector<int>& side : moveRepelSides_) {
            sides.insert(sides.end(), side.begin(), side.end());
        }
        if (sides.size() >= 2) {
            collisionGrid_.Clear();
            collisionGrid_.Reserve(sides.size());
            for (const int bead : sides) {
                const Vector2 p = PosOf(bead);
                collisionGrid_.Insert(bead, p.x, p.y);
            }
            collisionGrid_.ForEachCandidatePair([&](int a, int b) {
                applyPair(a, b);
            });
        }
    }
}

void Physics::Integrate(float dt)
{
    const ScopedTimer timer(debugTiming_ ? &debugTimings_.integrateMs : nullptr);
    const int beads = BeadCount();
    const int faces = static_cast<int>(neck_->faceBeads.size());
    const int total = beads + faces;
    std::vector<Vector2> force(total, Vector2{0.0f, 0.0f});

    // Springs (rest length = diameter).
    for (size_t s = 0; s < springs_.size(); ++s) {
        const int a = springs_[s].first;
        const int b = springs_[s].second;
        const Vector2 pa = PosOf(a);
        const Vector2 pb = PosOf(b);
        const Vector2 d = {pb.x - pa.x, pb.y - pa.y};
        const float len = std::sqrt(d.x * d.x + d.y * d.y);
        if (len < 1e-6f) {
            continue;
        }
        const float mag = params_->springStiffness * (len - diameter_);
        const Vector2 dir = {d.x / len, d.y / len};
        force[a].x += mag * dir.x;
        force[a].y += mag * dir.y;
        force[b].x -= mag * dir.x;
        force[b].y -= mag * dir.y;
    }

    // Temporary Reidemeister 3 springs: each sliding terminal is pulled
    // toward its current target bead along the face boundary.
    if (r3Active_) {
        for (const R3Slider& s : r3Sliders_) {
            if (s.terminal < 0 || s.path.empty()) {
                continue;
            }
            const int target = s.path[static_cast<size_t>(s.targetIdx)];
            const Vector2 pa = PosOf(s.terminal);
            const Vector2 pb = PosOf(target);
            const Vector2 d = {pb.x - pa.x, pb.y - pa.y};
            const float len = std::sqrt(d.x * d.x + d.y * d.y);
            if (len < 1e-6f) {
                continue;
            }
            const float mag = params_->springStiffness * (len - diameter_);
            const Vector2 dir = {d.x / len, d.y / len};
            force[s.terminal].x += mag * dir.x;
            force[s.terminal].y += mag * dir.y;
            force[target].x -= mag * dir.x;
            force[target].y -= mag * dir.y;
        }
    }

    // Inverse-square repulsion between interiors on distinct sides of the
    // moved face (R3 triangle, R2 bigon) so it stays open while stripped.
    const bool repelSides = (r3Active_ && r3Phase_ == R3Phase::PreHoneycomb) ||
                            (r2Active_ && shortenSub_ != ShortenSub::Spoke);
    if (repelSides && params_->r3CycleRepelStiffness > 0.0f) {
        RebuildMoveRepelSides();
        const float cutoff = 4.0f * diameter_;
        const float rFloor = 0.25f * diameter_;
        const float k = params_->r3CycleRepelStiffness;
        const int sides = static_cast<int>(moveRepelSides_.size());
        for (int i = 0; i < sides; ++i) {
            for (int j = i + 1; j < sides; ++j) {
                for (const int a : moveRepelSides_[static_cast<size_t>(i)]) {
                    if (a < 0 || a >= beads || neck_->beads[static_cast<size_t>(a)].removed) {
                        continue;
                    }
                    for (const int b : moveRepelSides_[static_cast<size_t>(j)]) {
                        if (b < 0 || b >= beads ||
                            neck_->beads[static_cast<size_t>(b)].removed) {
                            continue;
                        }
                        const Vector2 pa = PosOf(a);
                        const Vector2 pb = PosOf(b);
                        const Vector2 d = {pb.x - pa.x, pb.y - pa.y};
                        const float len = std::sqrt(d.x * d.x + d.y * d.y);
                        if (len < 1e-6f || len >= cutoff) {
                            continue;
                        }
                        const float r = std::max(len, rFloor);
                        const float scale = diameter_ / r;
                        const float mag = k * scale * scale;
                        const Vector2 dir = {d.x / len, d.y / len};
                        force[a].x -= mag * dir.x;
                        force[a].y -= mag * dir.y;
                        force[b].x += mag * dir.x;
                        force[b].y += mag * dir.y;
                    }
                }
            }
        }
    }

    std::vector<char> skipFaceCollision(static_cast<size_t>(faces), 0);
    std::vector<char> faceArrived(static_cast<size_t>(faces), 0);
    std::vector<char> haveIncenter(static_cast<size_t>(faces), 0);
    std::vector<Vector2> incenterTarget(static_cast<size_t>(faces), Vector2{0.0f, 0.0f});
    std::unordered_set<int> r3Faces;
    if (MoveActive()) {
        CollectMoveInvolvedFaces(&r3Faces);
    }
    // Crossing the polygon is not arrival: the disk still overlaps edge
    // beads and chatters. Collisions and honeycomb spokes stay off until the
    // center is inside and within 0.2*diameter of the inscribed-circle center.
    const float incenterArriveTol = 0.2f * diameter_;
    const double incenterPrecision = std::max(1e-4, 0.05 * static_cast<double>(radius_));
    for (int fi = 0; fi < faces; ++fi) {
        const FaceBead& fb = neck_->faceBeads[static_cast<size_t>(fi)];
        std::vector<Vec2d> poly;
        std::vector<std::vector<Vec2d>> holes;
        const bool havePoly = FaceBeadBoundaryPolygon(fb, &poly, &holes);
        const bool inside =
            !havePoly || PointInPolygonWithHoles(Vec2d{static_cast<double>(fb.pos.x),
                                                       static_cast<double>(fb.pos.y)},
                                                 poly, holes);
        if (!havePoly) {
            faceArrived[static_cast<size_t>(fi)] = inside ? 1 : 0;
            if (!inside) {
                skipFaceCollision[static_cast<size_t>(fi)] = 1;
            }
            continue;
        }
        const InscribedCircle circle = holes.empty()
                                           ? LargestInscribedCircle(poly, incenterPrecision)
                                           : LargestInscribedCircle(poly, holes, incenterPrecision);
        Vec2d target = circle.center;
        if (!circle.found) {
            target = {0.0, 0.0};
            for (const Vec2d& v : poly) {
                target.x += v.x;
                target.y += v.y;
            }
            target.x /= static_cast<double>(poly.size());
            target.y /= static_cast<double>(poly.size());
        }
        const Vector2 p = PosOf(beads + fi);
        const float dx = static_cast<float>(target.x) - p.x;
        const float dy = static_cast<float>(target.y) - p.y;
        const float len = std::sqrt(dx * dx + dy * dy);
        incenterTarget[static_cast<size_t>(fi)] = {static_cast<float>(target.x),
                                                   static_cast<float>(target.y)};
        haveIncenter[static_cast<size_t>(fi)] = 1;
        // No inscribed circle: centroid can sit on a pinch, so never "arrive".
        const bool arrived = inside && circle.found && len <= incenterArriveTol;
        faceArrived[static_cast<size_t>(fi)] = arrived ? 1 : 0;
        if (!arrived) {
            skipFaceCollision[static_cast<size_t>(fi)] = 1;
        }
    }

    // Hexagonal honeycomb spokes: each 6-cycle 1-/2-/3-gon face bead is
    // pulled to its boundary beads (rest = diameter). Off until arrival so
    // the truss cannot pin a bead against the necklace on the way in.
    for (size_t h = 0; h < honeycombFaceFi_.size(); ++h) {
        if (honeycombSpokes_.size() < (h + 1) * 6) {
            break;
        }
        const int fi = honeycombFaceFi_[h];
        const int faceBody = beads + fi;
        if (fi < 0 || fi >= faces || faceBody < 0 || faceBody >= total ||
            !faceArrived[static_cast<size_t>(fi)]) {
            continue;
        }
        for (int i = 0; i < 6; ++i) {
            const int spoke = honeycombSpokes_[h * 6 + static_cast<size_t>(i)];
            if (spoke < 0 || spoke >= beads ||
                neck_->beads[static_cast<size_t>(spoke)].removed) {
                continue;
            }
            const Vector2 pa = PosOf(faceBody);
            const Vector2 pb = PosOf(spoke);
            const Vector2 d = {pb.x - pa.x, pb.y - pa.y};
            const float len = std::sqrt(d.x * d.x + d.y * d.y);
            if (len < 1e-6f) {
                continue;
            }
            const float mag = params_->honeycombFortifyStiffness * (len - diameter_);
            const Vector2 dir = {d.x / len, d.y / len};
            force[faceBody].x += mag * dir.x;
            force[faceBody].y += mag * dir.y;
            force[spoke].x -= mag * dir.x;
            force[spoke].y -= mag * dir.y;
        }
    }

    for (int fi = 0; fi < faces; ++fi) {
        const FaceBead& fb = neck_->faceBeads[static_cast<size_t>(fi)];
        const bool arrived = faceArrived[static_cast<size_t>(fi)] != 0;
        bool pullToIncenter = !arrived;
        if (MoveActive() && r3Faces.count(fb.face) != 0) {
            pullToIncenter = true;
        }
        // Both beads of a pending merge chase the merged face's incenter.
        if (mergeKeepLabel_ >= 0 && (fb.label == mergeKeepLabel_ || fb.label == mergeDropLabel_)) {
            pullToIncenter = true;
        }
        if (FaceHasHoneycombSpokes(fi) && arrived) {
            pullToIncenter = false;
        }
        if (!pullToIncenter || !haveIncenter[static_cast<size_t>(fi)]) {
            continue;
        }
        const int body = beads + fi;
        const Vector2 p = PosOf(body);
        const Vector2 target = incenterTarget[static_cast<size_t>(fi)];
        const Vector2 d = {target.x - p.x, target.y - p.y};
        const float len = std::sqrt(d.x * d.x + d.y * d.y);
        if (len > 1e-6f) {
            // Rest length 0: force = stiffness * displacement toward the incenter.
            const float mag = params_->springStiffness * len;
            const Vector2 dir = {d.x / len, d.y / len};
            force[body].x += mag * dir.x;
            force[body].y += mag * dir.y;
        }
    }

    // R2 spoke: pull the bigon's two vertex beads together (rest = diameter)
    // so they become tangent before the crossings are removed.
    if (r2Active_ && shortenSub_ == ShortenSub::Spoke) {
        const int x = r2Info_.vertexBead[0];
        const int y = r2Info_.vertexBead[1];
        if (x >= 0 && x < beads && y >= 0 && y < beads && x != y &&
            !neck_->beads[static_cast<size_t>(x)].removed &&
            !neck_->beads[static_cast<size_t>(y)].removed) {
            const Vector2 pa = PosOf(x);
            const Vector2 pb = PosOf(y);
            const Vector2 d = {pb.x - pa.x, pb.y - pa.y};
            const float len = std::sqrt(d.x * d.x + d.y * d.y);
            if (len > 1e-6f) {
                const float mag = params_->honeycombFortifyStiffness * (len - diameter_);
                const Vector2 dir = {d.x / len, d.y / len};
                force[x].x += mag * dir.x;
                force[x].y += mag * dir.y;
                force[y].x -= mag * dir.x;
                force[y].y -= mag * dir.y;
            }
        }
    }

    // Hard-body collisions within each face (spatial broad phase).
    ApplyFaceCollisions(force, total, skipFaceCollision);

    // Drag toward the cursor.
    if (dragBody_ >= 0 && dragBody_ < total) {
        const Vector2 p = PosOf(dragBody_);
        force[dragBody_].x += params_->dragStiffness * (dragTarget_.x - p.x);
        force[dragBody_].y += params_->dragStiffness * (dragTarget_.y - p.y);
    }

    // Damping + semi-implicit Euler. Edge beads have mass 1; during R3, face
    // beads use r3FaceBeadMass so they stay put while the necklace slides.
    for (int body = 0; body < total; ++body) {
        if (body < beads && neck_->beads[body].removed) {
            continue;
        }
        float mass = 1.0f;
        if (MoveActive() && body >= beads) {
            mass = std::max(params_->r3FaceBeadMass, 1e-3f);
        }
        Vector2& v = VelOf(body);
        Vector2 f = force[body];
        f.x -= params_->damping * v.x;
        f.y -= params_->damping * v.y;
        v.x += (f.x / mass) * dt;
        v.y += (f.y / mass) * dt;
        const float speed = std::sqrt(v.x * v.x + v.y * v.y);
        if (speed > params_->maxSpeed) {
            v.x *= params_->maxSpeed / speed;
            v.y *= params_->maxSpeed / speed;
        }
        Vector2 p = PosOf(body);
        p.x += v.x * dt;
        p.y += v.y * dt;
        SetPosOf(body, p);
    }
}

void Physics::ComputeTension()
{
    stretch_ = 0.0;
    compression_ = 0.0;
    int beadCount = 0;
    for (const Bead& bead : neck_->beads) {
        if (!bead.removed) {
            ++beadCount;
        }
    }
    for (int k = 1; k <= loop_->numEdges; ++k) {
        const std::vector<int>& chain = chains_[k];
        for (size_t i = 0; i + 1 < chain.size(); ++i) {
            const Vector2 a = PosOf(chain[i]);
            const Vector2 b = PosOf(chain[i + 1]);
            const float len = std::sqrt((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y));
            const float delta = len - diameter_;
            const double tension = params_->springStiffness * std::fabs(delta);
            if (delta > 0.0f) {
                stretch_ += tension;
            }
            else {
                compression_ += tension;
            }
        }
    }
    if (beadCount > 0) {
        stretch_ /= beadCount;
        compression_ /= beadCount;
    }
}

void Physics::RefreshTangle()
{
    tangled_ = DetectTangle();
}

bool Physics::DetectTangle() const
{
    const int n = BeadCount();
    for (int i = 0; i < n; ++i) {
        const Bead& b = neck_->beads[static_cast<size_t>(i)];
        if (!b.isVertex || b.removed) {
            continue;
        }
        Vec2d p[4];
        bool ok = true;
        for (int k = 0; k < 4; ++k) {
            const int nb = b.neighborOfHalfEdge[k];
            if (nb < 0 || nb >= n || neck_->beads[static_cast<size_t>(nb)].removed) {
                ok = false;
                break;
            }
            const BeadPos& q = neck_->beads[static_cast<size_t>(nb)].pos;
            p[k] = Vec2d{static_cast<double>(q.x), static_cast<double>(q.y)};
        }
        if (!ok) {
            continue;
        }
        const Vec2d c{static_cast<double>(b.pos.x), static_cast<double>(b.pos.y)};
        if (!OppositeSlotsLinked(c, p[0], p[1], p[2], p[3])) {
            return true;
        }
    }

    struct Seg
    {
        int a = -1;
        int b = -1;
        Vector2 pa = {0.0f, 0.0f};
        Vector2 pb = {0.0f, 0.0f};
    };
    std::vector<Seg> segs;
    segs.reserve(static_cast<size_t>(n) * 2);
    for (int i = 0; i < n; ++i) {
        const Bead& b = neck_->beads[static_cast<size_t>(i)];
        if (b.removed) {
            continue;
        }
        int neigh[4];
        int nc = 0;
        if (b.isVertex) {
            for (int k = 0; k < 4; ++k) {
                neigh[nc++] = b.neighborOfHalfEdge[k];
            }
        }
        else {
            neigh[nc++] = b.neighbors[0];
            neigh[nc++] = b.neighbors[1];
        }
        for (int k = 0; k < nc; ++k) {
            const int j = neigh[k];
            if (j <= i || j >= n || neck_->beads[static_cast<size_t>(j)].removed) {
                continue;
            }
            Seg s;
            s.a = i;
            s.b = j;
            s.pa = Vector2{b.pos.x, b.pos.y};
            const BeadPos& q = neck_->beads[static_cast<size_t>(j)].pos;
            s.pb = Vector2{q.x, q.y};
            segs.push_back(s);
        }
    }

    float maxLen = diameter_;
    for (const Seg& s : segs) {
        const float dx = s.pb.x - s.pa.x;
        const float dy = s.pb.y - s.pa.y;
        const float len = std::sqrt(dx * dx + dy * dy);
        if (len > maxLen) {
            maxLen = len;
        }
    }

    SpatialDiskGrid grid;
    grid.SetCellSize(std::max(3.0f * diameter_, maxLen));
    grid.Reserve(segs.size());
    for (int i = 0; i < static_cast<int>(segs.size()); ++i) {
        const float mx = 0.5f * (segs[static_cast<size_t>(i)].pa.x +
                                 segs[static_cast<size_t>(i)].pb.x);
        const float my = 0.5f * (segs[static_cast<size_t>(i)].pa.y +
                                 segs[static_cast<size_t>(i)].pb.y);
        grid.Insert(i, mx, my);
    }

    bool hit = false;
    grid.ForEachCandidatePair([&](int ia, int ib) {
        if (hit) {
            return;
        }
        const Seg& a = segs[static_cast<size_t>(ia)];
        const Seg& b = segs[static_cast<size_t>(ib)];
        if (a.a == b.a || a.a == b.b || a.b == b.a || a.b == b.b) {
            return;
        }
        if (SegmentsProperlyCross(a.pa, a.pb, b.pa, b.pb, nullptr)) {
            hit = true;
        }
    });
    return hit;
}

double Physics::EdgeTension(int edge) const
{
    if (edge < 1 || edge > loop_->numEdges) {
        return 0.0;
    }
    const std::vector<int>& chain = chains_[edge];
    if (chain.size() < 2) {
        return 0.0;
    }
    double stretch = 0.0;
    double compression = 0.0;
    for (size_t i = 0; i + 1 < chain.size(); ++i) {
        const Vector2 a = PosOf(chain[i]);
        const Vector2 b = PosOf(chain[i + 1]);
        const float len = std::sqrt((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y));
        const float delta = len - diameter_;
        const double tension = params_->springStiffness * std::fabs(delta);
        if (delta > 0.0f) {
            stretch += tension;
        }
        else {
            compression += tension;
        }
    }
    const int beadCount = static_cast<int>(chain.size());
    stretch /= beadCount;
    compression /= beadCount;
    return stretch + compression;
}

int Physics::SubstepsForStep() const
{
    if (!autoActive_ || pendingEdges_.empty()) {
        return params_->substeps;
    }
    bool anyUnsettled = false;
    bool allRelaxing = true;
    for (const PendingEdge& pending : pendingEdges_) {
        if (pending.settled) {
            continue;
        }
        anyUnsettled = true;
        const double threshold =
            params_->tensionThreshold * std::sqrt(static_cast<double>(pending.length));
        if (PendingTension(pending) > threshold * 0.5) {
            allRelaxing = false;
            break;
        }
    }
    if (!anyUnsettled || !allRelaxing) {
        return params_->substeps;
    }
    return std::max(1, params_->substeps / 4);
}

void Physics::Step()
{
    const ScopedTimer timer(debugTiming_ ? &debugTimings_.stepMs : nullptr);
    if (debugTiming_) {
        ++debugTimings_.stepCount;
    }
    if (dragBody_ >= 0 && dragBody_ < BeadCount()) {
        MaintainBeadSpacingDuringDrag();
    }
    if (structureDirty_) {
        if (r3Active_ && r3Phase_ == R3Phase::Slide) {
            RebuildStructuresR3();
        }
        else {
            RebuildStructures();
        }
    }
    const int substeps = SubstepsForStep();
    if (debugTiming_) {
        debugTimings_.substepsTotal += substeps;
    }
    for (int s = 0; s < substeps; ++s) {
        Integrate(params_->dt);
    }
    ComputeTension();
}

float Physics::Dot(int edgeBead) const
{
    const Bead& b = neck_->beads[edgeBead];
    const Vector2 center = {b.pos.x, b.pos.y};
    const Vector2 a = PosOf(b.neighbors[0]);
    const Vector2 c = PosOf(b.neighbors[1]);
    const Vector2 ba = {a.x - center.x, a.y - center.y};
    const Vector2 bc = {c.x - center.x, c.y - center.y};
    const float la = std::sqrt(ba.x * ba.x + ba.y * ba.y);
    const float lc = std::sqrt(bc.x * bc.x + bc.y * bc.y);
    if (la < 1e-6f || lc < 1e-6f) {
        return -1.0f;
    }
    return (ba.x * bc.x + ba.y * bc.y) / (la * lc);
}

int Physics::PickBody(Vector2 gridPos, float maxGridDist) const
{
    int best = -1;
    float bestDist = maxGridDist;
    for (int i = 0; i < BeadCount(); ++i) {
        if (neck_->beads[i].removed) {
            continue;
        }
        const Vector2 p = {neck_->beads[i].pos.x, neck_->beads[i].pos.y};
        const float d = std::sqrt((p.x - gridPos.x) * (p.x - gridPos.x) +
                                  (p.y - gridPos.y) * (p.y - gridPos.y));
        if (d < bestDist) {
            bestDist = d;
            best = i;
        }
    }
    for (size_t i = 0; i < neck_->faceBeads.size(); ++i) {
        const Vector2 p = {neck_->faceBeads[i].pos.x, neck_->faceBeads[i].pos.y};
        const float d = std::sqrt((p.x - gridPos.x) * (p.x - gridPos.x) +
                                  (p.y - gridPos.y) * (p.y - gridPos.y));
        if (d < bestDist) {
            bestDist = d;
            best = BeadCount() + static_cast<int>(i);
        }
    }
    return best;
}

void Physics::BeginDrag(int body, Vector2 gridPos)
{
    dragBody_ = body;
    dragTarget_ = gridPos;
}

void Physics::Drag(Vector2 gridPos)
{
    dragTarget_ = gridPos;
}

void Physics::EndDrag()
{
    dragBody_ = -1;
}

void Physics::ReplaceNeighbor(int node, int oldNeighbor, int newNeighbor)
{
    Bead& b = neck_->beads[node];
    if (b.isVertex) {
        for (int i = 0; i < 4; ++i) {
            if (b.neighborOfHalfEdge[i] == oldNeighbor) {
                b.neighborOfHalfEdge[i] = newNeighbor;
            }
        }
    }
    else {
        for (int i = 0; i < 2; ++i) {
            if (b.neighbors[i] == oldNeighbor) {
                b.neighbors[i] = newNeighbor;
            }
        }
    }
}

// Like ReplaceNeighbor, but only updates the first matching slot. Needed when
// reinserting into a kink (loop edge): both +k and -k at the shared vertex
// may temporarily point at the same bead, and blindly rewriting every match
// collapses the two ends of the loop onto one neighbor.
void Physics::ReplaceNeighborOnce(int node, int oldNeighbor, int newNeighbor)
{
    Bead& b = neck_->beads[node];
    if (b.isVertex) {
        for (int i = 0; i < 4; ++i) {
            if (b.neighborOfHalfEdge[i] == oldNeighbor) {
                b.neighborOfHalfEdge[i] = newNeighbor;
                return;
            }
        }
    }
    else {
        for (int i = 0; i < 2; ++i) {
            if (b.neighbors[i] == oldNeighbor) {
                b.neighbors[i] = newNeighbor;
                return;
            }
        }
    }
}

void Physics::ReinsertEdgeBead(int bead, int neighborA, int neighborC)
{
    if (bead < 0 || bead >= BeadCount() || neighborA < 0 || neighborC < 0 ||
        neighborA == neighborC || neighborA >= BeadCount() || neighborC >= BeadCount()) {
        return;
    }
    // A and C must still be linked to each other; otherwise splicing would
    // create a one-way orphan (the failure mode behind floating kink beads).
    if (!NodeReferences(neighborA, neighborC) || !NodeReferences(neighborC, neighborA)) {
        return;
    }
    if (CountReferences(neighborA, neighborC) != 1 ||
        CountReferences(neighborC, neighborA) != 1) {
        return;
    }
    Bead& b = neck_->beads[static_cast<size_t>(bead)];
    const Vector2 pa = PosOf(neighborA);
    const Vector2 pc = PosOf(neighborC);
    b.removed = false;
    b.pos = {(pa.x + pc.x) * 0.5f, (pa.y + pc.y) * 0.5f};
    b.neighbors[0] = neighborA;
    b.neighbors[1] = neighborC;
    ReplaceNeighborOnce(neighborA, neighborC, bead);
    ReplaceNeighborOnce(neighborC, neighborA, bead);
    if (bead < static_cast<int>(vel_.size())) {
        vel_[static_cast<size_t>(bead)] = Vector2{0.0f, 0.0f};
    }
    structureDirty_ = true;
}

void Physics::RemoveEdgeBead(int bead)
{
    if (bead < 0 || bead >= BeadCount()) {
        return;
    }
    Bead& b = neck_->beads[static_cast<size_t>(bead)];
    if (b.removed || b.isVertex) {
        return;
    }
    const int a = b.neighbors[0];
    const int c = b.neighbors[1];
    // Last bead on a loop edge has neighbors (v, v). Splicing that out
    // would write the vertex's own index into its neighbor slots.
    if (a == c || a < 0 || c < 0) {
        return;
    }
    // Each endpoint must reference this bead in exactly one slot. Two slots
    // (collapsed kink) or zero (already unlinked) means refuse.
    if (CountReferences(a, bead) != 1 || CountReferences(c, bead) != 1) {
        return;
    }
    ReplaceNeighborOnce(a, bead, c);
    ReplaceNeighborOnce(c, bead, a);
    b.removed = true;
    structureDirty_ = true;
}

int Physics::InsertEdgeBeadBetween(int neighborA, int neighborC, int edge)
{
    if (neighborA < 0 || neighborC < 0 || neighborA == neighborC ||
        !NodeReferences(neighborA, neighborC) || !NodeReferences(neighborC, neighborA) ||
        CountReferences(neighborA, neighborC) != 1 ||
        CountReferences(neighborC, neighborA) != 1) {
        return -1;
    }
    const int strandA = neck_->beads[static_cast<size_t>(neighborA)].strand;
    const int strandC = neck_->beads[static_cast<size_t>(neighborC)].strand;
    const int strand = (strandA >= 0) ? strandA : strandC;

    for (int i = 0; i < BeadCount(); ++i) {
        Bead& slot = neck_->beads[static_cast<size_t>(i)];
        if (!slot.removed || slot.isVertex || slot.edge != edge) {
            continue;
        }
        ReinsertEdgeBead(i, neighborA, neighborC);
        neck_->beads[static_cast<size_t>(i)].strand = strand;
        return i;
    }

    Bead bead;
    bead.isVertex = false;
    bead.removed = false;
    bead.edge = edge;
    bead.strand = strand;
    const Vector2 pa = PosOf(neighborA);
    const Vector2 pc = PosOf(neighborC);
    bead.pos = {(pa.x + pc.x) * 0.5f, (pa.y + pc.y) * 0.5f};
    bead.neighbors[0] = neighborA;
    bead.neighbors[1] = neighborC;
    const int idx = BeadCount();
    ReplaceNeighborOnce(neighborA, neighborC, idx);
    ReplaceNeighborOnce(neighborC, neighborA, idx);
    neck_->beads.push_back(std::move(bead));
    vel_.push_back(Vector2{0.0f, 0.0f});
    structureDirty_ = true;
    return idx;
}

void Physics::MaintainBeadSpacingDuringDrag()
{
    const Bead& dragged = neck_->beads[static_cast<size_t>(dragBody_)];
    if (dragged.removed) {
        return;
    }

    const float gapThreshold = params_->gapSize * diameter_;

    if (!dragged.isVertex && dragged.edge == 0) {
        FaceBead* pad = nullptr;
        for (FaceBead& fb : neck_->faceBeads) {
            if (!IsCircleFaceBead(fb)) {
                continue;
            }
            for (const int bead : fb.boundaryBeads) {
                if (bead == dragBody_) {
                    pad = &fb;
                    break;
                }
            }
            if (pad != nullptr) {
                break;
            }
        }

        std::vector<int> ring;
        if (pad != nullptr && pad->boundaryBeads.size() >= 2) {
            ring = pad->boundaryBeads;
        }
        else {
            int prev = -1;
            int cur = dragBody_;
            const int limit = BeadCount() + 2;
            for (int guard = 0; guard < limit; ++guard) {
                if (cur < 0 || cur >= BeadCount() ||
                    neck_->beads[static_cast<size_t>(cur)].removed) {
                    break;
                }
                ring.push_back(cur);
                const Bead& b = neck_->beads[static_cast<size_t>(cur)];
                int nxt = -1;
                if (prev < 0) {
                    nxt = (b.neighbors[1] >= 0) ? b.neighbors[1] : b.neighbors[0];
                }
                else {
                    nxt = (b.neighbors[0] == prev) ? b.neighbors[1] : b.neighbors[0];
                }
                prev = cur;
                cur = nxt;
                if (cur == dragBody_) {
                    break;
                }
            }
        }

        int gapA = -1;
        int gapB = -1;
        for (size_t i = 0; i < ring.size(); ++i) {
            const int a = ring[i];
            const int b = ring[(i + 1) % ring.size()];
            if (a < 0 || b < 0 || a >= BeadCount() || b >= BeadCount()) {
                continue;
            }
            if (neck_->beads[static_cast<size_t>(a)].removed ||
                neck_->beads[static_cast<size_t>(b)].removed) {
                continue;
            }
            const Vector2 pa = PosOf(a);
            const Vector2 pb = PosOf(b);
            const float dx = pb.x - pa.x;
            const float dy = pb.y - pa.y;
            const float len = std::sqrt(dx * dx + dy * dy);
            if (len > gapThreshold) {
                gapA = a;
                gapB = b;
                break;
            }
        }
        if (gapA < 0) {
            return;
        }
        if (InsertEdgeBeadBetween(gapA, gapB, 0) < 0) {
            return;
        }
        if (pad != nullptr) {
            RefreshCircleBoundary(*pad);
        }
        else {
            for (FaceBead& fb : neck_->faceBeads) {
                if (IsCircleFaceBead(fb)) {
                    RefreshCircleBoundary(fb);
                }
            }
        }
        RefreshCircleStructures();
        return;
    }

    std::vector<int> edgesToCheck;
    if (dragged.isVertex) {
        edgesToCheck.reserve(dragged.sigmaCycle.size());
        for (const int h : dragged.sigmaCycle) {
            edgesToCheck.push_back(std::abs(h));
        }
    }
    else {
        if (dragged.edge < 1 || dragged.edge > loop_->numEdges) {
            return;
        }
        edgesToCheck.push_back(dragged.edge);
    }

    if (structureDirty_) {
        if (r3Active_ && r3Phase_ == R3Phase::Slide) {
            RebuildStructuresR3();
        }
        else {
            RebuildStructures();
        }
    }

    int gapEdge = -1;
    int gapA = -1;
    int gapB = -1;
    for (const int edge : edgesToCheck) {
        if (edge < 1 || edge > loop_->numEdges) {
            continue;
        }
        const std::vector<int>& chain = chains_[edge];
        bool onChain = false;
        for (const int bead : chain) {
            if (bead == dragBody_) {
                onChain = true;
                break;
            }
        }
        if (!onChain) {
            continue;
        }

        for (size_t i = 0; i + 1 < chain.size(); ++i) {
            const int a = chain[i];
            const int b = chain[i + 1];
            const Vector2 pa = PosOf(a);
            const Vector2 pb = PosOf(b);
            const float dx = pb.x - pa.x;
            const float dy = pb.y - pa.y;
            const float len = std::sqrt(dx * dx + dy * dy);
            if (len > gapThreshold) {
                gapEdge = edge;
                gapA = a;
                gapB = b;
                break;
            }
        }
        if (gapEdge >= 0) {
            break;
        }
    }

    if (gapEdge < 0) {
        return;
    }

    InsertEdgeBeadBetween(gapA, gapB, gapEdge);
}

Vector2 Physics::FaceBeadPosition(int face) const
{
    for (size_t fi = 0; fi < neck_->faceBeads.size(); ++fi) {
        if (neck_->faceBeads[fi].face == face) {
            return PosOf(BeadCount() + static_cast<int>(fi));
        }
    }
    return Vector2{0.0f, 0.0f};
}

FaceBead* Physics::FaceBeadFor(int face)
{
    FaceBead* found = nullptr;
    for (FaceBead& fb : neck_->faceBeads) {
        if (fb.face != face) {
            continue;
        }
        // While a merge is pending, prefer the bead that survives.
        if (mergeKeepLabel_ >= 0 && fb.label == mergeKeepLabel_) {
            return &fb;
        }
        if (found == nullptr) {
            found = &fb;
        }
    }
    return found;
}

const FaceBead* Physics::FaceBeadFor(int face) const
{
    const FaceBead* found = nullptr;
    for (const FaceBead& fb : neck_->faceBeads) {
        if (fb.face != face) {
            continue;
        }
        if (mergeKeepLabel_ >= 0 && fb.label == mergeKeepLabel_) {
            return &fb;
        }
        if (found == nullptr) {
            found = &fb;
        }
    }
    return found;
}

bool Physics::FaceBoundaryPolygon(int face, std::vector<Vec2d>* out) const
{
    if (out == nullptr || face < 0 || face >= static_cast<int>(loop_->faces.size()) ||
        IsExterior(face)) {
        return false;
    }
    out->clear();
    std::vector<int> boundary;
    std::vector<float> unusedDots;
    WalkFaceBoundaryFromPhiCycle(loop_->faces[static_cast<size_t>(face)], *neck_,
                                 crossingOf_, boundary, unusedDots, nullptr, 0.0f, false);
    if (boundary.size() < 3) {
        const FaceBead* fb = FaceBeadFor(face);
        if (fb != nullptr) {
            boundary = fb->boundaryBeads;
        }
    }
    out->reserve(boundary.size());
    for (const int bead : boundary) {
        if (bead < 0 || bead >= BeadCount() || neck_->beads[static_cast<size_t>(bead)].removed) {
            continue;
        }
        const Vector2 p = PosOf(bead);
        out->push_back(Vec2d{static_cast<double>(p.x), static_cast<double>(p.y)});
    }
    return out->size() >= 3;
}

void Physics::CollectMoveInvolvedFaces(std::unordered_set<int>* out) const
{
    if (out == nullptr) {
        return;
    }
    out->clear();
    const auto addHalf = [&](int h) {
        const auto it = faceOfHalfEdge_.find(h);
        if (it != faceOfHalfEdge_.end() && !IsExterior(it->second)) {
            out->insert(it->second);
        }
    };
    const auto addEdge = [&](int h) {
        addHalf(h);
        addHalf(-h);
    };
    if (r3Active_) {
        if (r3Info_.face >= 0 && !IsExterior(r3Info_.face)) {
            out->insert(r3Info_.face);
        }
        for (int i = 0; i < 3; ++i) {
            addEdge(r3Info_.h[i]);
            addEdge(r3Info_.a[i]);
            addEdge(r3Info_.b[i]);
        }
    }
    if (r1Active_) {
        addEdge(r1Info_.e);
        addEdge(r1Info_.p);
        addEdge(r1Info_.q);
    }
    if (r2Active_) {
        addEdge(r2Info_.h0);
        addEdge(r2Info_.h1);
        addEdge(r2Info_.a);
        addEdge(r2Info_.b);
        addEdge(r2Info_.ap);
        addEdge(r2Info_.bp);
    }
}

bool Physics::FaceBeadBoundaryPolygon(const FaceBead& fb, std::vector<Vec2d>* out) const
{
    if (fb.face >= 0) {
        return FaceBoundaryPolygon(fb.face, out);
    }
    if (out == nullptr) {
        return false;
    }
    out->clear();
    out->reserve(fb.boundaryBeads.size());
    for (const int bead : fb.boundaryBeads) {
        if (bead < 0 || bead >= BeadCount() || neck_->beads[static_cast<size_t>(bead)].removed) {
            continue;
        }
        const Vector2 p = PosOf(bead);
        out->push_back(Vec2d{static_cast<double>(p.x), static_cast<double>(p.y)});
    }
    return out->size() >= 3;
}

bool Physics::FaceBeadBoundaryPolygon(const FaceBead& fb, std::vector<Vec2d>* out,
                                      std::vector<std::vector<Vec2d>>* holes) const
{
    if (holes != nullptr) {
        holes->clear();
    }
    if (!FaceBeadBoundaryPolygon(fb, out)) {
        return false;
    }
    if (holes == nullptr) {
        return true;
    }
    for (const std::vector<int>& ring : fb.holeBeads) {
        std::vector<Vec2d> poly;
        poly.reserve(ring.size());
        for (const int bead : ring) {
            if (bead < 0 || bead >= BeadCount() ||
                neck_->beads[static_cast<size_t>(bead)].removed) {
                continue;
            }
            const Vector2 p = PosOf(bead);
            poly.push_back(Vec2d{static_cast<double>(p.x), static_cast<double>(p.y)});
        }
        if (poly.size() >= 3) {
            holes->push_back(std::move(poly));
        }
    }
    return true;
}

bool Physics::FaceBeadInsideFace(const FaceBead& fb) const
{
    std::vector<Vec2d> poly;
    std::vector<std::vector<Vec2d>> holes;
    if (!FaceBeadBoundaryPolygon(fb, &poly, &holes)) {
        return true;
    }
    const Vec2d center{static_cast<double>(fb.pos.x), static_cast<double>(fb.pos.y)};
    return PointInPolygonWithHoles(center, poly, holes);
}

int Physics::PickRegionAt(Vector2 gridPos) const
{
    const Vec2d p{static_cast<double>(gridPos.x), static_cast<double>(gridPos.y)};
    int best = -1;
    double bestArea = 0.0;
    for (size_t fi = 0; fi < neck_->faceBeads.size(); ++fi) {
        const FaceBead& fb = neck_->faceBeads[fi];
        std::vector<Vec2d> poly;
        std::vector<std::vector<Vec2d>> holes;
        if (!FaceBeadBoundaryPolygon(fb, &poly, &holes)) {
            continue;
        }
        if (!PointInPolygonWithHoles(p, poly, holes)) {
            continue;
        }
        const double area = ShoelaceArea(poly);
        if (best < 0 || area < bestArea) {
            best = static_cast<int>(fi);
            bestArea = area;
        }
    }
    return best;
}

bool Physics::IsCircleFaceBead(const FaceBead& fb) const
{
    return fb.face < 0 || unknot_;
}

int Physics::LiveBoundaryCount(const FaceBead& fb) const
{
    int live = 0;
    for (const int bead : fb.boundaryBeads) {
        if (bead >= 0 && bead < BeadCount() && !neck_->beads[static_cast<size_t>(bead)].removed) {
            ++live;
        }
    }
    return live;
}

int Physics::FarthestLiveBoundaryBead(const FaceBead& fb) const
{
    const Vector2 center = {fb.pos.x, fb.pos.y};
    int best = -1;
    float bestDist = -1.0f;
    for (const int bead : fb.boundaryBeads) {
        if (bead < 0 || bead >= BeadCount()) {
            continue;
        }
        const Bead& b = neck_->beads[static_cast<size_t>(bead)];
        if (b.removed || b.isVertex) {
            continue;
        }
        if (b.neighbors[0] == b.neighbors[1] || b.neighbors[0] < 0 || b.neighbors[1] < 0) {
            continue;
        }
        const Vector2 p = PosOf(bead);
        const float dx = p.x - center.x;
        const float dy = p.y - center.y;
        const float dist = dx * dx + dy * dy;
        if (dist > bestDist) {
            bestDist = dist;
            best = bead;
        }
    }
    return best;
}

void Physics::RefreshCircleBoundary(FaceBead& fb)
{
    int start = -1;
    for (const int bead : fb.boundaryBeads) {
        if (bead >= 0 && bead < BeadCount() && !neck_->beads[static_cast<size_t>(bead)].removed) {
            start = bead;
            break;
        }
    }
    fb.boundaryBeads.clear();
    if (start < 0) {
        return;
    }
    int prev = -1;
    int cur = start;
    const int limit = BeadCount() + 2;
    for (int guard = 0; guard < limit; ++guard) {
        if (cur < 0 || cur >= BeadCount() || neck_->beads[static_cast<size_t>(cur)].removed) {
            break;
        }
        fb.boundaryBeads.push_back(cur);
        const Bead& b = neck_->beads[static_cast<size_t>(cur)];
        int nxt = -1;
        if (prev < 0) {
            nxt = (b.neighbors[1] >= 0) ? b.neighbors[1] : b.neighbors[0];
        }
        else {
            nxt = (b.neighbors[0] == prev) ? b.neighbors[1] : b.neighbors[0];
        }
        prev = cur;
        cur = nxt;
        if (cur == start) {
            break;
        }
    }
    fb.RefreshBoundaryDots(*neck_);
    fb.cachedExcessArea = fb.ExcessArea(*neck_, radius_);
}

void Physics::RefreshCircleStructures()
{
    if (unknot_) {
        RebuildUnknotStructures();
        return;
    }
    RebuildSpringsFromChains();
    RebuildFaceBodiesFromCache();
    structureDirty_ = false;
}

double Physics::CircleRingTension(const FaceBead& fb) const
{
    const std::vector<int>& ring = fb.boundaryBeads;
    if (ring.size() < 2) {
        return 0.0;
    }
    double stretch = 0.0;
    double compression = 0.0;
    int segments = 0;
    for (size_t i = 0; i < ring.size(); ++i) {
        const int a = ring[i];
        const int c = ring[(i + 1) % ring.size()];
        if (a < 0 || c < 0 || a >= BeadCount() || c >= BeadCount()) {
            continue;
        }
        if (neck_->beads[static_cast<size_t>(a)].removed ||
            neck_->beads[static_cast<size_t>(c)].removed) {
            continue;
        }
        const Vector2 pa = PosOf(a);
        const Vector2 pc = PosOf(c);
        const float len = std::sqrt((pc.x - pa.x) * (pc.x - pa.x) + (pc.y - pa.y) * (pc.y - pa.y));
        const float delta = len - diameter_;
        const double tension = params_->springStiffness * std::fabs(delta);
        if (delta > 0.0f) {
            stretch += tension;
        }
        else {
            compression += tension;
        }
        ++segments;
    }
    if (segments <= 0) {
        return 0.0;
    }
    const double beadCount = static_cast<double>(segments);
    return (stretch + compression) / beadCount;
}

float Physics::CircleRingLength(const FaceBead& fb) const
{
    const std::vector<int>& ring = fb.boundaryBeads;
    float length = 0.0f;
    for (size_t i = 0; i < ring.size(); ++i) {
        const int a = ring[i];
        const int c = ring[(i + 1) % ring.size()];
        if (a < 0 || c < 0 || a >= BeadCount() || c >= BeadCount()) {
            continue;
        }
        if (neck_->beads[static_cast<size_t>(a)].removed ||
            neck_->beads[static_cast<size_t>(c)].removed) {
            continue;
        }
        const Vector2 pa = PosOf(a);
        const Vector2 pc = PosOf(c);
        length += std::sqrt((pc.x - pa.x) * (pc.x - pa.x) + (pc.y - pa.y) * (pc.y - pa.y));
    }
    return length;
}

bool Physics::RemoveCircleExcessIfAllowed(bool oneBead)
{
    if (unknot_) {
        if (structureDirty_) {
            RebuildUnknotStructures();
        }
    }
    else if (structureDirty_) {
        RebuildStructures();
    }

    bool any = false;
    for (FaceBead& fb : neck_->faceBeads) {
        if (!IsCircleFaceBead(fb)) {
            continue;
        }
        if (LiveBoundaryCount(fb) <= std::max(6, HoleFloor(fb))) {
            continue;
        }
        const int bead = FarthestLiveBoundaryBead(fb);
        if (bead < 0) {
            continue;
        }
        RemoveEdgeBead(bead);
        if (!neck_->beads[static_cast<size_t>(bead)].removed) {
            continue;
        }
        RefreshCircleBoundary(fb);
        any = true;
        if (oneBead) {
            break;
        }
    }
    if (!any) {
        return false;
    }
    RefreshCircleStructures();
    return true;
}

void Physics::QueueCircleSettlePendings()
{
    for (size_t fi = 0; fi < neck_->faceBeads.size(); ++fi) {
        const FaceBead& fb = neck_->faceBeads[fi];
        if (!IsCircleFaceBead(fb) || fb.boundaryBeads.size() < 2) {
            continue;
        }
        float length = 0.0f;
        const std::vector<int>& ring = fb.boundaryBeads;
        for (size_t i = 0; i < ring.size(); ++i) {
            const int a = ring[i];
            const int c = ring[(i + 1) % ring.size()];
            if (a < 0 || c < 0 || a >= BeadCount() || c >= BeadCount()) {
                continue;
            }
            if (neck_->beads[static_cast<size_t>(a)].removed ||
                neck_->beads[static_cast<size_t>(c)].removed) {
                continue;
            }
            const Vector2 pa = PosOf(a);
            const Vector2 pc = PosOf(c);
            length += std::sqrt((pc.x - pa.x) * (pc.x - pa.x) + (pc.y - pa.y) * (pc.y - pa.y));
        }
        PendingEdge pending;
        pending.edge = 0;
        pending.faceBead = static_cast<int>(fi);
        pending.length = length;
        pendingEdges_.push_back(pending);
    }
}

double Physics::PendingTension(const PendingEdge& pending) const
{
    if (pending.faceBead >= 0 &&
        pending.faceBead < static_cast<int>(neck_->faceBeads.size())) {
        return CircleRingTension(neck_->faceBeads[static_cast<size_t>(pending.faceBead)]);
    }
    return EdgeTension(pending.edge);
}

double Physics::FaceExcessArea(const FaceBead& fb) const
{
    if (fb.boundaryBeads.size() < 3) {
        return 0.0;
    }
    fb.RefreshBoundaryDots(*neck_);
    return fb.ExcessArea(*neck_, radius_);
}

double Physics::ExcessArea(int face) const
{
    if (IsExterior(face)) {
        return std::numeric_limits<double>::infinity();
    }
    const FaceBead* fb = FaceBeadFor(face);
    if (fb == nullptr) {
        return 0.0;
    }
    return FaceExcessArea(*fb);
}

double Physics::IsoperimetricRatio(int face) const
{
    if (IsExterior(face)) {
        return std::numeric_limits<double>::infinity();
    }
    const FaceBead* fb = FaceBeadFor(face);
    if (fb == nullptr) {
        return 0.0;
    }
    if (fb->boundaryBeads.size() < 3) {
        return 0.0;
    }
    fb->RefreshBoundaryDots(*neck_);
    return fb->IsoperimetricRatio(*neck_);
}

void Physics::StrikeHoneycombEdges(bool boundariesCurrent,
                                   const std::unordered_set<int>& refreshFaces)
{
    const ScopedTimer timer(debugTiming_ ? &debugTimings_.strikeHoneycombMs : nullptr);
    if (struckEdges_.size() != static_cast<size_t>(loop_->numEdges + 1)) {
        struckEdges_.assign(loop_->numEdges + 1, false);
    }
    for (int f = 0; f < static_cast<int>(loop_->faces.size()); ++f) {
        if (IsExterior(f)) {
            continue;
        }
        FaceBead* fb = FaceBeadFor(f);
        if (fb == nullptr) {
            continue;
        }
        if (!boundariesCurrent && refreshFaces.count(f) != 0) {
            WalkFaceBoundaryNoDots(f);
        }
        // Honeycomb faces keep their 6 beads; faces with holes keep enough
        // boundary to stay clear of the nested rings.
        const int live = static_cast<int>(fb->boundaryBeads.size());
        const bool honeycomb = live <= 6 && fb->holeBeads.empty();
        const bool holeFloor = !fb->holeBeads.empty() && live <= HoleFloor(*fb);
        if (!honeycomb && !holeFloor) {
            continue;
        }
        for (const int h : loop_->faces[f]) {
            const int k = std::abs(h);
            if (k >= 1 && k < static_cast<int>(struckEdges_.size())) {
                struckEdges_[k] = true;
            }
        }
    }
}

namespace
{
    void PrintIntCycle(const char* label, const std::vector<int>& cycle)
    {
        std::printf("%s(", label);
        for (size_t i = 0; i < cycle.size(); ++i) {
            if (i > 0) {
                std::printf(",");
            }
            std::printf("%d", cycle[i]);
        }
        std::printf(")");
    }
}

void Physics::PrintShoelaceDebug(int face) const
{
    if (IsExterior(face)) {
        std::printf("[shoelace] face %d is exterior (infinite area)\n", face);
        return;
    }

    const FaceBead* fb = FaceBeadFor(face);
    if (fb == nullptr) {
        return;
    }
    if (fb->boundaryBeads.empty()) {
        fb->WalkFaceBoundary(*loop_, *neck_, crossingOf_);
    }
    else {
        fb->RefreshBoundaryDots(*neck_);
    }
    std::printf("[shoelace] face %d polygon (%zu vertices):\n", face,
                fb->boundaryBeads.size());

    for (size_t i = 0; i < fb->boundaryBeads.size(); ++i) {
        const int bead = fb->boundaryBeads[i];
        const Bead& b = neck_->beads[static_cast<size_t>(bead)];
        const Vector2 p = PosOf(bead);

        std::printf("  [%zu] bead %d ", i, bead);
        if (b.isVertex) {
            std::printf("vertex (crossing %d) ", b.crossing);
            PrintIntCycle("sigma", b.sigmaCycle);
        }
        else {
            std::printf("edge (edge %d)", b.edge);
            if (b.removed) {
                std::printf(" [removed]");
            }
        }
        std::printf(" pos=(%.6f, %.6f)", p.x, p.y);
        if (i < fb->boundaryDots.size()) {
            std::printf(" dot=%.6f", fb->boundaryDots[i]);
        }
        std::printf("\n");
    }

    std::printf("[shoelace] polyArea=%.6f excessArea=%.6f ratio=%.6f\n",
                fb->PolyArea(*neck_), fb->ExcessArea(*neck_, radius_),
                fb->IsoperimetricRatio(*neck_));
}

double Physics::EdgeExcessScore(int edge) const
{
    std::vector<double> excesses;
    const auto itPos = faceOfHalfEdge_.find(edge);
    const auto itNeg = faceOfHalfEdge_.find(-edge);
    if (itPos != faceOfHalfEdge_.end() && !IsExterior(itPos->second)) {
        const FaceBead* fb = FaceBeadFor(itPos->second);
        if (fb != nullptr) {
            excesses.push_back(fb->cachedExcessArea);
        }
    }
    if (itNeg != faceOfHalfEdge_.end() && !IsExterior(itNeg->second)) {
        const FaceBead* fb = FaceBeadFor(itNeg->second);
        if (fb != nullptr) {
            excesses.push_back(fb->cachedExcessArea);
        }
    }
    if (excesses.empty()) {
        return 0.0;
    }
    if (excesses.size() == 1) {
        return excesses[0];
    }
    switch (params_->ratioCombineMode) {
    case RatioCombineMode::Max:
        return std::max(excesses[0], excesses[1]);
    case RatioCombineMode::Avg:
        return 0.5 * (excesses[0] + excesses[1]);
    case RatioCombineMode::Min:
    default:
        return std::min(excesses[0], excesses[1]);
    }
}

int Physics::RemoveBeadFromEdge(int edge, bool bendPassOnly, int* neighborA, int* neighborC,
                                bool skipExcessCheck)
{
    if (edge < 1 || edge > loop_->numEdges) {
        return -1;
    }
    const std::vector<int>& chain = chains_[edge];
    const int edgeBeads = static_cast<int>(chain.size()) - 2;
    // Loop edges (both ends the same crossing) must keep at least two edge
    // beads: reducing to one merges the +k/-k neighbor slots at that vertex
    // and subsequent reinserts corrupt the chain into a cycle.
    const bool loopEdge = EdgeIsLoop(edge);
    const int minBeads = loopEdge ? 3 : 2;
    if (edgeBeads < minBeads) {
        return -1;
    }
    const int strand = StrandOfEdgeLabel(edge);
    if (strand >= 0 && StrandLiveEdgeBeadCount(strand) <= kMinStrandEdgeBeads) {
        return -1;
    }

    const float maxDot = sharpDotThreshold_;

    if (bendPassOnly) {
        for (size_t i = 1; i + 1 < chain.size(); ++i) {
            if (Dot(chain[i]) >= maxDot) {
                const int bead = chain[i];
                const int a = neck_->beads[bead].neighbors[0];
                const int c = neck_->beads[bead].neighbors[1];
                if (a == c) {
                    continue;
                }
                RemoveEdgeBead(bead);
                if (!neck_->beads[static_cast<size_t>(bead)].removed) {
                    continue;
                }
                if (neighborA != nullptr) {
                    *neighborA = a;
                }
                if (neighborC != nullptr) {
                    *neighborC = c;
                }
                return bead;
            }
        }
        return -1;
    }

    const auto facePos = faceOfHalfEdge_.find(edge);
    const auto faceNeg = faceOfHalfEdge_.find(-edge);
    if (facePos == faceOfHalfEdge_.end() || faceNeg == faceOfHalfEdge_.end()) {
        return -1;
    }
    if (!skipExcessCheck && EdgeExcessScore(edge) <= params_->minExcessArea) {
        return -1;
    }

    const Vector2 facePosCenter = FaceBeadPosition(facePos->second);
    const Vector2 faceNegCenter = FaceBeadPosition(faceNeg->second);

    int bestBead = -1;
    float bestScore = -1.0f;
    for (size_t i = 1; i + 1 < chain.size(); ++i) {
        const int bead = chain[i];
        const Vector2 p = PosOf(bead);
        const float dPos = std::sqrt((p.x - facePosCenter.x) * (p.x - facePosCenter.x) +
                                     (p.y - facePosCenter.y) * (p.y - facePosCenter.y));
        const float dNeg = std::sqrt((p.x - faceNegCenter.x) * (p.x - faceNegCenter.x) +
                                     (p.y - faceNegCenter.y) * (p.y - faceNegCenter.y));
        const float score = std::min(dPos, dNeg);
        if (score > bestScore) {
            bestScore = score;
            bestBead = bead;
        }
    }
    if (bestBead < 0) {
        return -1;
    }
    const int a = neck_->beads[bestBead].neighbors[0];
    const int c = neck_->beads[bestBead].neighbors[1];
    if (a == c) {
        return -1;
    }
    RemoveEdgeBead(bestBead);
    if (!neck_->beads[static_cast<size_t>(bestBead)].removed) {
        return -1;
    }
    if (neighborA != nullptr) {
        *neighborA = a;
    }
    if (neighborC != nullptr) {
        *neighborC = c;
    }
    return bestBead;
}

bool Physics::BeginRemovalRound()
{
    const ScopedTimer timer(debugTiming_ ? &debugTimings_.beginRemovalRoundMs : nullptr);
    if (unknot_) {
        pendingEdges_.clear();
        if (structureDirty_) {
            RebuildUnknotStructures();
        }
        if (!RemoveCircleExcessIfAllowed()) {
            return false;
        }
        QueueCircleSettlePendings();
        framesWaited_ = 0;
        return !pendingEdges_.empty();
    }
    const bool evaluateGeometry = needsGeometryEval_;

    std::unordered_set<int> honeycombFaces;
    for (const int k : removedEdgesLastRound_) {
        AppendIncidentFaces(k, honeycombFaces);
    }
    removedEdgesLastRound_.clear();

    if (evaluateGeometry) {
        RebuildStructures();
    }
    else if (structureDirty_) {
        const int n = loop_->numEdges;
        chains_.assign(n + 1, {});
        for (int k = 1; k <= n; ++k) {
            chains_[k] = WalkEdge(k);
        }
        RebuildSpringsFromChains();
        structureDirty_ = false;
    }

    StrikeHoneycombEdges(evaluateGeometry, honeycombFaces);

    if (!evaluateGeometry) {
        for (const int f : honeycombFaces) {
            RefreshFaceBody(f);
        }
    }

    const auto edgeLength = [&](const std::vector<int>& chain) {
        float length = 0.0f;
        for (size_t i = 0; i + 1 < chain.size(); ++i) {
            const Vector2 a = PosOf(chain[i]);
            const Vector2 b = PosOf(chain[i + 1]);
            length += std::sqrt((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y));
        }
        return length;
    };

    std::vector<int> excessEdges;
    std::vector<bool> hasSharpBead(static_cast<size_t>(loop_->numEdges + 1), false);
    for (const int bead : neck_->sharpCornerEdges) {
        if (bead < 0 || bead >= BeadCount() || neck_->beads[bead].removed) {
            continue;
        }
        const int k = neck_->beads[bead].edge;
        if (k >= 1 && k <= loop_->numEdges) {
            hasSharpBead[static_cast<size_t>(k)] = true;
        }
    }
    for (int k = 1; k <= loop_->numEdges; ++k) {
        if (k < static_cast<int>(struckEdges_.size()) && struckEdges_[k]) {
            continue;
        }
        const std::vector<int>& chain = chains_[k];
        if (static_cast<int>(chain.size()) - 2 < 2) {
            continue;
        }
        if (!hasSharpBead[static_cast<size_t>(k)]) {
            excessEdges.push_back(k);
        }
    }

    pendingEdges_.clear();
    bool anyRemoval = false;
    std::unordered_set<int> edgesRemovedThisRound;

    std::vector<int> edgeBeadsRemaining(static_cast<size_t>(loop_->numEdges + 1), 0);
    for (int k = 1; k <= loop_->numEdges; ++k) {
        edgeBeadsRemaining[static_cast<size_t>(k)] =
            static_cast<int>(chains_[k].size()) - 2;
    }

    const std::vector<int> sharpRemovals = neck_->sharpCornerEdges;
    for (const int bead : sharpRemovals) {
        if (bead < 0 || bead >= BeadCount() || neck_->beads[bead].removed ||
            neck_->beads[bead].isVertex) {
            continue;
        }
        const int k = neck_->beads[bead].edge;
        if (k < 1 || k > loop_->numEdges ||
            (k < static_cast<int>(struckEdges_.size()) && struckEdges_[k])) {
            continue;
        }
        const bool loopEdge = EdgeIsLoop(k);
        const int minBeads = loopEdge ? 3 : 2;
        if (edgeBeadsRemaining[static_cast<size_t>(k)] < minBeads) {
            continue;
        }
        // Also refuse beads whose two neighbors are the same vertex.
        if (neck_->beads[bead].neighbors[0] == neck_->beads[bead].neighbors[1]) {
            continue;
        }
        const float length = edgeLength(chains_[k]);
        const int a = neck_->beads[bead].neighbors[0];
        const int c = neck_->beads[bead].neighbors[1];
        RemoveEdgeBead(bead);
        if (neck_->beads[bead].removed == false) {
            continue; // RemoveEdgeBead refused (e.g. a == c)
        }
        PendingEdge pending;
        pending.edge = k;
        pending.length = length;
        pending.removedBead = bead;
        pending.neighborA = a;
        pending.neighborC = c;
        pendingEdges_.push_back(pending);
        edgeBeadsRemaining[static_cast<size_t>(k)]--;
        RefreshEdgeChain(k);
        edgesRemovedThisRound.insert(k);
        anyRemoval = true;
    }

    const bool skipExcessCheck = !evaluateGeometry;
    for (const int k : excessEdges) {
        const bool loopEdge = EdgeIsLoop(k);
        const int minBeads = loopEdge ? 3 : 2;
        if (edgeBeadsRemaining[static_cast<size_t>(k)] < minBeads) {
            continue;
        }
        const float length = edgeLength(chains_[k]);
        int neighborA = -1;
        int neighborC = -1;
        const int removed =
            RemoveBeadFromEdge(k, false, &neighborA, &neighborC, skipExcessCheck);
        if (removed < 0) {
            continue;
        }
        PendingEdge pending;
        pending.edge = k;
        pending.length = length;
        pending.removedBead = removed;
        pending.neighborA = neighborA;
        pending.neighborC = neighborC;
        pendingEdges_.push_back(pending);
        edgeBeadsRemaining[static_cast<size_t>(k)]--;
        RefreshEdgeChain(k);
        edgesRemovedThisRound.insert(k);
        anyRemoval = true;
    }

    if (anyRemoval) {
        std::unordered_set<int> affectedFaces;
        for (const int k : edgesRemovedThisRound) {
            AppendIncidentFaces(k, affectedFaces);
        }
        // Full spring rebuild: safer than per-edge refresh after kink edits.
        RebuildSpringsFromChains();
        for (const int f : affectedFaces) {
            WalkFaceBoundaryNoDots(f);
            RefreshFaceBody(f);
        }
        removedEdgesLastRound_ = std::move(edgesRemovedThisRound);
        structureDirty_ = false;
    }

    if (evaluateGeometry) {
        needsGeometryEval_ = false;
    }

    const bool circleRemoved = RemoveCircleExcessIfAllowed();
    if (circleRemoved) {
        QueueCircleSettlePendings();
    }

    if (pendingEdges_.empty()) {
        return false;
    }

    framesWaited_ = 0;
    return true;
}

void Physics::StartAutoTighten()
{
    if (autoActive_ || MoveActive()) {
        return;
    }
    autoActive_ = true;
    needsGeometryEval_ = !unknot_;
    removedEdgesLastRound_.clear();
    if (!unknot_) {
        struckEdges_.assign(loop_->numEdges + 1, false);
    }
    pendingEdges_.clear();
    framesWaited_ = 0;
    if (!BeginRemovalRound()) {
        autoActive_ = false;
    }
}

void Physics::StopAutoTighten()
{
    // Restore any beads that were removed but not yet settled, so cancelling
    // mid-round cannot leave the necklace with dangling cuts.
    for (const PendingEdge& pending : pendingEdges_) {
        if (!pending.settled && pending.removedBead >= 0 &&
            pending.removedBead < BeadCount() &&
            neck_->beads[static_cast<size_t>(pending.removedBead)].removed) {
            ReinsertEdgeBead(pending.removedBead, pending.neighborA, pending.neighborC);
        }
    }
    autoActive_ = false;
    needsGeometryEval_ = false;
    pendingEdges_.clear();
    structureDirty_ = true;
}

void Physics::AutoStep()
{
    if (!autoActive_) {
        return;
    }
    if (pendingEdges_.empty()) {
        autoActive_ = false;
        return;
    }

    {
        const ScopedTimer settleTimer(debugTiming_ ? &debugTimings_.autoStepOverheadMs : nullptr);
        for (int i = 0; i < params_->autoStepsPerFrame; ++i) {
            Step();
        }
        framesWaited_ += params_->autoStepsPerFrame;

    std::vector<int> failedEdges;
    for (PendingEdge& pending : pendingEdges_) {
        if (pending.settled) {
            continue;
        }
        const double threshold =
            params_->tensionThreshold * std::sqrt(static_cast<double>(pending.length));
        if (PendingTension(pending) < threshold) {
            pending.settled = true;
        }
        else if (framesWaited_ >= params_->relaxStepsMax) {
            failedEdges.push_back(pending.edge);
            pending.settled = true;
        }
    }

        if (!failedEdges.empty()) {
            std::unordered_set<int> reinsertedEdges;
            for (const int k : failedEdges) {
                for (PendingEdge& pending : pendingEdges_) {
                    if (pending.edge != k || pending.removedBead < 0) {
                        continue;
                    }
                    // Circle strips never go back above 6 beads.
                    if (pending.faceBead >= 0 || pending.edge == 0) {
                        pending.settled = true;
                        continue;
                    }
                    ReinsertEdgeBead(pending.removedBead, pending.neighborA, pending.neighborC);
                    pending.settled = true;
                    reinsertedEdges.insert(k);
                }
                if (k >= 1 && k < static_cast<int>(struckEdges_.size())) {
                    struckEdges_[k] = true;
                }
            }
            for (const int k : reinsertedEdges) {
                RefreshEdgeChain(k);
                RefreshFaceBodiesForEdge(k);
                removedEdgesLastRound_.insert(k);
            }
            RebuildSpringsFromChains();
            structureDirty_ = false;
            needsGeometryEval_ = true;
        }

        bool allSettledNow = true;
        for (const PendingEdge& pending : pendingEdges_) {
            if (!pending.settled) {
                allSettledNow = false;
                break;
            }
        }
        if (allSettledNow) {
            pendingEdges_.clear();
            if (!BeginRemovalRound()) {
                autoActive_ = false;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Reidemeister 3 (honeycomb-tighten, then slide)
// ---------------------------------------------------------------------------

void Physics::ClearR3State()
{
    for (R3Slider& s : r3Sliders_) {
        s = R3Slider{};
    }
    r3CutEdges_.clear();
    r3ExteriorEdges_.clear();
    moveRepelSides_.clear();
    r3Meeting_[0] = r3Meeting_[1] = r3Meeting_[2] = -1;
    r3StallFrames_ = 0;
    r3SettleEdge_ = -1;
    r3SettleLength_ = 0.0f;
    r3SettleFrames_ = 0;
    r3Phase_ = R3Phase::PreHoneycomb;
    r3PreSub_ = R3PreSub::NeedRemove;
}

std::vector<int> Physics::WalkSideChain(int vertexBead, int halfEdge) const
{
    std::vector<int> chain;
    chain.push_back(vertexBead);
    int prev = vertexBead;
    int cur = neck_->NeighborToward(vertexBead, halfEdge);
    int guard = 0;
    const int limit = BeadCount() + 8;
    while (cur >= 0 && !neck_->beads[static_cast<size_t>(cur)].isVertex) {
        if (neck_->beads[static_cast<size_t>(cur)].removed) {
            break;
        }
        chain.push_back(cur);
        const Bead& b = neck_->beads[static_cast<size_t>(cur)];
        const int nxt = (b.neighbors[0] == prev) ? b.neighbors[1] : b.neighbors[0];
        if (nxt == prev) {
            break;
        }
        prev = cur;
        cur = nxt;
        if (++guard > limit) {
            break;
        }
    }
    if (cur >= 0 && neck_->beads[static_cast<size_t>(cur)].isVertex) {
        chain.push_back(cur);
    }
    return chain;
}

void Physics::RebuildMoveRepelSides()
{
    moveRepelSides_.clear();
    std::vector<std::pair<int, int>> sides; // (vertex bead, half-edge leaving it)
    if (r3Active_) {
        for (int i = 0; i < 3; ++i) {
            sides.emplace_back(r3Info_.vertexBead[i], r3Info_.h[i]);
        }
    }
    else if (r2Active_) {
        sides.emplace_back(r2Info_.vertexBead[0], r2Info_.h0);
        sides.emplace_back(r2Info_.vertexBead[1], r2Info_.h1);
    }
    for (const std::pair<int, int>& side : sides) {
        const int vb = side.first;
        const int h = side.second;
        std::vector<int> interiors;
        if (vb >= 0) {
            const std::vector<int> chain = WalkSideChain(vb, h);
            for (size_t t = 1; t + 1 < chain.size(); ++t) {
                const int b = chain[t];
                if (b < 0 || b >= BeadCount() || neck_->beads[static_cast<size_t>(b)].removed) {
                    continue;
                }
                interiors.push_back(b);
            }
        }
        moveRepelSides_.push_back(std::move(interiors));
    }
}

int Physics::EdgeInteriorCount(int edge) const
{
    if (edge < 1 || edge >= static_cast<int>(chains_.size())) {
        return 0;
    }
    const std::vector<int>& chain = chains_[edge];
    if (chain.size() < 2) {
        return 0;
    }
    return static_cast<int>(chain.size()) - 2;
}

float Physics::EdgeChainLength(int edge) const
{
    if (edge < 1 || edge >= static_cast<int>(chains_.size())) {
        return 0.0f;
    }
    const std::vector<int>& chain = chains_[edge];
    float length = 0.0f;
    for (size_t i = 0; i + 1 < chain.size(); ++i) {
        const Vector2 a = PosOf(chain[i]);
        const Vector2 b = PosOf(chain[i + 1]);
        length += std::sqrt((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y));
    }
    return length;
}

int Physics::InsertBeadAtFarthestGap(int edge)
{
    if (structureDirty_) {
        if (r3Phase_ == R3Phase::Slide) {
            RebuildStructuresR3();
        }
        else {
            RebuildStructures();
        }
    }
    else {
        RefreshEdgeChain(edge);
    }
    if (edge < 1 || edge >= static_cast<int>(chains_.size())) {
        return -1;
    }
    const std::vector<int>& chain = chains_[edge];
    if (chain.size() < 2) {
        return -1;
    }
    size_t best = 0;
    float bestDist = -1.0f;
    for (size_t i = 0; i + 1 < chain.size(); ++i) {
        const Vector2 a = PosOf(chain[i]);
        const Vector2 b = PosOf(chain[i + 1]);
        const float d = std::sqrt((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y));
        if (d > bestDist) {
            bestDist = d;
            best = i;
        }
    }
    const int beadsBefore = BeadCount();
    const int inserted = InsertEdgeBeadBetween(chain[best], chain[best + 1], edge);
    if (inserted >= 0) {
        RefreshEdgeChain(edge);
        RebuildSpringsFromChains();
        if (r3Phase_ == R3Phase::Slide) {
            // Collision groups stay frozen during the slide; only adjust face
            // body ids if the necklace grew, then add the new bead locally.
            if (BeadCount() > beadsBefore) {
                const int delta = BeadCount() - beadsBefore;
                for (std::vector<int>& group : faceBodies_) {
                    for (int& body : group) {
                        if (body >= beadsBefore) {
                            body += delta;
                        }
                    }
                }
            }
            AddBeadToFrozenGroups(inserted, chain[best]);
        }
        else {
            // Phase A (and any non-slide use): refresh boundaries and rebuild
            // faceBodies_ so face-bead ids stay BeadCount()+fi.
            RefreshFaceBodiesForEdge(edge);
            RebuildFaceBodiesFromCache();
        }
        structureDirty_ = false;
    }
    return inserted;
}

int Physics::RemoveFarthestFromFaceCenter(int edge, int face)
{
    if (structureDirty_) {
        RebuildStructures();
    }
    return RemoveFarthestFromPoint(edge, FaceBeadPosition(face));
}

int Physics::RemoveFarthestFromPoint(int edge, Vector2 fc)
{
    if (structureDirty_) {
        RebuildStructures();
    }
    else {
        RefreshEdgeChain(edge);
    }
    if (edge < 1 || edge >= static_cast<int>(chains_.size())) {
        return -1;
    }
    const std::vector<int>& chain = chains_[edge];
    if (static_cast<int>(chain.size()) - 2 < 2) {
        return -1; // need more than one interior to remove
    }
    int best = -1;
    float bestDist = -1.0f;
    for (size_t i = 1; i + 1 < chain.size(); ++i) {
        const int bead = chain[i];
        if (neck_->beads[static_cast<size_t>(bead)].removed) {
            continue;
        }
        const Vector2 p = PosOf(bead);
        const float d = std::sqrt((p.x - fc.x) * (p.x - fc.x) + (p.y - fc.y) * (p.y - fc.y));
        if (d > bestDist) {
            bestDist = d;
            best = bead;
        }
    }
    if (best < 0) {
        return -1;
    }
    RemoveEdgeBead(best);
    if (!neck_->beads[static_cast<size_t>(best)].removed) {
        return -1;
    }
    RefreshEdgeChain(edge);
    RebuildSpringsFromChains();
    RefreshFaceBodiesForEdge(edge);
    RebuildFaceBodiesFromCache();
    structureDirty_ = false;
    return best;
}

bool Physics::R3SettleTick(int edge, float lengthAtStart)
{
    for (int i = 0; i < params_->autoStepsPerFrame; ++i) {
        Step();
    }
    r3SettleFrames_ += params_->autoStepsPerFrame;
    if (structureDirty_) {
        RebuildStructures();
    }
    else {
        RefreshEdgeChain(edge);
    }
    const double threshold =
        params_->tensionThreshold * std::sqrt(static_cast<double>(std::max(lengthAtStart, 1.0f)));
    if (EdgeTension(edge) < threshold) {
        return true;
    }
    if (r3SettleFrames_ >= params_->relaxStepsMax) {
        return true;
    }
    return false;
}

bool Physics::NeighborHoneycombBlocksRemoval(int edge)
{
    // Match StrikeHoneycombEdges: a face with <=6 boundary beads must not lose
    // more edge beads. Protect neighbors of the clicked triangle only.
    std::unordered_set<int> faces;
    AppendIncidentFaces(edge, faces);
    for (const int f : faces) {
        if (f == r3Info_.face || IsExterior(f)) {
            continue;
        }
        WalkFaceBoundaryNoDots(f);
        const FaceBead* fb = FaceBeadFor(f);
        if (fb != nullptr && fb->boundaryBeads.size() <= 6) {
            return true;
        }
    }
    return false;
}

bool Physics::HoneycombFloorBlocksRemoval(int edge)
{
    // Only the infinite face itself being a monogon (phi length 1) must keep
    // its long exterior arc. A finite monogon that merely touches the outer
    // face should still be allowed to tighten after R3.
    if (EdgeIsLoop(edge) && loop_ != nullptr) {
        const auto itPos = faceOfHalfEdge_.find(edge);
        const auto itNeg = faceOfHalfEdge_.find(-edge);
        const auto exteriorMonogon = [&](int face) {
            return IsExterior(face) && face >= 0 &&
                   face < static_cast<int>(loop_->faces.size()) &&
                   loop_->faces[static_cast<size_t>(face)].size() == 1;
        };
        if ((itPos != faceOfHalfEdge_.end() && exteriorMonogon(itPos->second)) ||
            (itNeg != faceOfHalfEdge_.end() && exteriorMonogon(itNeg->second))) {
            return true;
        }
    }
    std::unordered_set<int> faces;
    AppendIncidentFaces(edge, faces);
    for (const int f : faces) {
        if (IsExterior(f)) {
            continue;
        }
        WalkFaceBoundaryNoDots(f);
        const FaceBead* fb = FaceBeadFor(f);
        if (fb == nullptr) {
            continue;
        }
        const int live = static_cast<int>(fb->boundaryBeads.size());
        if (live <= 6 || (!fb->holeBeads.empty() && live <= HoleFloor(*fb))) {
            return true;
        }
    }
    return false;
}

int Physics::HoleFloor(const FaceBead& fb) const
{
    if (fb.holeBeads.empty()) {
        return 0;
    }
    // Enough perimeter for the nested rings plus an annulus one bead wide
    // (2*pi*2 diameters ~ 12 beads) so the region's own bead still fits.
    int live = 0;
    for (const std::vector<int>& ring : fb.holeBeads) {
        for (const int b : ring) {
            if (b >= 0 && b < BeadCount() && !neck_->beads[static_cast<size_t>(b)].removed) {
                ++live;
            }
        }
    }
    return live + 12;
}

bool Physics::InsertLargestGapIfNeeded()
{
    if (structureDirty_) {
        RebuildStructures();
    }
    const float gapThreshold = params_->r3PostRepairGapSize * diameter_;
    float bestDist = gapThreshold;
    int bestEdge = -1;
    int bestA = -1;
    int bestC = -1;
    for (int k = 1; k <= loop_->numEdges; ++k) {
        if (k >= static_cast<int>(chains_.size())) {
            break;
        }
        const std::vector<int>& chain = chains_[k];
        for (size_t i = 0; i + 1 < chain.size(); ++i) {
            const Vector2 a = PosOf(chain[i]);
            const Vector2 b = PosOf(chain[i + 1]);
            const float d =
                std::sqrt((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y));
            if (d > bestDist) {
                bestDist = d;
                bestEdge = k;
                bestA = chain[i];
                bestC = chain[i + 1];
            }
        }
    }
    if (bestEdge < 0) {
        return false;
    }
    const int inserted = InsertEdgeBeadBetween(bestA, bestC, bestEdge);
    if (inserted < 0) {
        return false;
    }
    RefreshEdgeChain(bestEdge);
    RebuildSpringsFromChains();
    RefreshFaceBodiesForEdge(bestEdge);
    RebuildFaceBodiesFromCache();
    structureDirty_ = false;
    return true;
}

bool Physics::RemoveSharpBendIfAllowed()
{
    if (structureDirty_) {
        RebuildStructures();
    }
    sharpDotThreshold_ = std::cos(kPiF / 3.0f + params_->bendAngleEpsilon);
    for (int k = 1; k <= loop_->numEdges; ++k) {
        if (HoneycombFloorBlocksRemoval(k)) {
            continue;
        }
        const int removed = RemoveBeadFromEdge(k, true);
        if (removed < 0) {
            continue;
        }
        RefreshEdgeChain(k);
        RebuildSpringsFromChains();
        RefreshFaceBodiesForEdge(k);
        RebuildFaceBodiesFromCache();
        structureDirty_ = false;
        return true;
    }
    return false;
}

bool Physics::RemoveMonogonExcessIfAllowed()
{
    if (structureDirty_) {
        RebuildStructures();
    }
    for (int k = 1; k <= loop_->numEdges; ++k) {
        if (!EdgeIsLoop(k) || HoneycombFloorBlocksRemoval(k)) {
            continue;
        }
        const int removed = RemoveBeadFromEdge(k, false, nullptr, nullptr, true);
        if (removed < 0) {
            continue;
        }
        RefreshEdgeChain(k);
        RebuildSpringsFromChains();
        RefreshFaceBodiesForEdge(k);
        RebuildFaceBodiesFromCache();
        structureDirty_ = false;
        return true;
    }
    return false;
}

void Physics::StartPostRepair()
{
    if (autoActive_ || MoveActive()) {
        return;
    }
    ClearR3State();
    r3Active_ = true;
    r3Complete_ = false;
    r3Phase_ = R3Phase::PostRepair;
    r3PreSub_ = R3PreSub::NeedInsert;
    structureDirty_ = true;
}

void Physics::R3PostRepairStep()
{
    TryFaceMerge();
    if (++r3StallFrames_ > 800) {
        r3Active_ = false;
        r3Complete_ = false;
        ClearR3State();
        structureDirty_ = true;
        return;
    }
    if (unknot_) {
        // No crossings left: ring excess is shed passively (IdleStep).
        Step();
        if (!FaceMergePending()) {
            r3Active_ = false;
            r3Complete_ = false;
            ClearR3State();
            structureDirty_ = true;
        }
        return;
    }
    if (r3PreSub_ == R3PreSub::NeedInsert) {
        int n = 0;
        bool any = false;
        while (n++ < 16 && InsertLargestGapIfNeeded()) {
            any = true;
        }
        Step();
        if (!any) {
            r3PreSub_ = R3PreSub::NeedRemove;
        }
        return;
    }

    // Monogon / circle excess is not stripped here: it is shed slowly in the
    // resting state (see IdleStep) so the move ends without a snap.
    int n = 0;
    bool any = false;
    while (n++ < 8 && RemoveSharpBendIfAllowed()) {
        any = true;
    }
    Step();
    if (!any && !FaceMergePending()) {
        r3Active_ = false;
        r3Complete_ = false;
        ClearR3State();
        structureDirty_ = true;
    }
}

// ---------------------------------------------------------------------------
// Resting state: passive tightening of monogons (to 5 beads) and circles
// (to 6 beads), one bead every params_->passiveTightenInterval frames.
// ---------------------------------------------------------------------------

void Physics::IdleStep()
{
    Step();
    if (params_->passiveTightenInterval <= 0 || MoveActive() || autoActive_ ||
        dragBody_ >= 0 || FaceMergePending()) {
        idleFrames_ = 0;
        return;
    }
    if (++idleFrames_ < params_->passiveTightenInterval) {
        return;
    }
    idleFrames_ = 0;
    PassiveTightenOnce();
}

bool Physics::PassiveTightenOnce()
{
    if (unknot_) {
        return RemoveCircleExcessIfAllowed(true);
    }
    if (loop_ == nullptr || loop_->numEdges <= 0) {
        return false;
    }
    if (RemoveMonogonExcessIfAllowed()) {
        return true;
    }
    return RemoveCircleExcessIfAllowed(true);
}

// ---------------------------------------------------------------------------
// Face-bead merge (connected R2): two beads share the merged face.
// ---------------------------------------------------------------------------

void Physics::SetFaceMergePair(int keepLabel, int dropLabel)
{
    mergeKeepLabel_ = -1;
    mergeDropLabel_ = -1;
    mergeFrames_ = 0;
    if (keepLabel < 0 || dropLabel < 0 || keepLabel == dropLabel) {
        return;
    }
    const FaceBead* keep = nullptr;
    const FaceBead* drop = nullptr;
    for (const FaceBead& fb : neck_->faceBeads) {
        if (fb.label == keepLabel) {
            keep = &fb;
        }
        else if (fb.label == dropLabel) {
            drop = &fb;
        }
    }
    if (keep == nullptr || drop == nullptr || keep->face < 0 || keep->face != drop->face) {
        return;
    }
    mergeKeepLabel_ = keepLabel;
    mergeDropLabel_ = dropLabel;
    structureDirty_ = true;
}

void Physics::TryFaceMerge()
{
    if (mergeKeepLabel_ < 0) {
        return;
    }
    int keepIdx = -1;
    int dropIdx = -1;
    for (size_t fi = 0; fi < neck_->faceBeads.size(); ++fi) {
        const int label = neck_->faceBeads[fi].label;
        if (label == mergeKeepLabel_) {
            keepIdx = static_cast<int>(fi);
        }
        else if (label == mergeDropLabel_) {
            dropIdx = static_cast<int>(fi);
        }
    }
    if (keepIdx < 0 || dropIdx < 0) {
        mergeKeepLabel_ = -1;
        mergeDropLabel_ = -1;
        mergeFrames_ = 0;
        structureDirty_ = true;
        return;
    }
    ++mergeFrames_;
    const Vector2 a = PosOf(BeadCount() + keepIdx);
    const Vector2 b = PosOf(BeadCount() + dropIdx);
    const float dx = a.x - b.x;
    const float dy = a.y - b.y;
    const float dist = std::sqrt(dx * dx + dy * dy);
    constexpr int kMergeFrameCap = 600;
    if (dist > 0.25f * diameter_ && mergeFrames_ < kMergeFrameCap) {
        return;
    }
    neck_->faceBeads[static_cast<size_t>(keepIdx)].highlighted =
        neck_->faceBeads[static_cast<size_t>(keepIdx)].highlighted ||
        neck_->faceBeads[static_cast<size_t>(dropIdx)].highlighted;
    neck_->faceBeads.erase(neck_->faceBeads.begin() + dropIdx);
    if (static_cast<size_t>(dropIdx) < fvel_.size()) {
        fvel_.erase(fvel_.begin() + dropIdx);
    }
    mergeKeepLabel_ = -1;
    mergeDropLabel_ = -1;
    mergeFrames_ = 0;
    structureDirty_ = true;
}

// ---------------------------------------------------------------------------
// Face-bead erase / restore for shortening moves
// ---------------------------------------------------------------------------

bool Physics::EraseFaceBeadAt(size_t fi)
{
    if (fi >= neck_->faceBeads.size()) {
        return false;
    }
    erasedFaceBead_ = neck_->faceBeads[fi];
    erasedFaceBeadIndex_ = fi;
    erasedFaceBeadValid_ = true;
    neck_->faceBeads.erase(neck_->faceBeads.begin() + static_cast<long>(fi));
    if (fi < fvel_.size()) {
        fvel_.erase(fvel_.begin() + static_cast<long>(fi));
    }
    structureDirty_ = true;
    return true;
}

bool Physics::EraseFaceBead(int face)
{
    for (size_t fi = 0; fi < neck_->faceBeads.size(); ++fi) {
        if (neck_->faceBeads[fi].face != face) {
            continue;
        }
        return EraseFaceBeadAt(fi);
    }
    return false;
}

void Physics::RestoreErasedFaceBead()
{
    if (!erasedFaceBeadValid_) {
        return;
    }
    const size_t at = std::min(erasedFaceBeadIndex_, neck_->faceBeads.size());
    neck_->faceBeads.insert(neck_->faceBeads.begin() + static_cast<long>(at), erasedFaceBead_);
    if (fvel_.size() + 1 == neck_->faceBeads.size()) {
        fvel_.insert(fvel_.begin() + static_cast<long>(std::min(at, fvel_.size())),
                     Vector2{0.0f, 0.0f});
    }
    else {
        fvel_.assign(neck_->faceBeads.size(), Vector2{0.0f, 0.0f});
    }
    erasedFaceBeadValid_ = false;
    structureDirty_ = true;
}

bool Physics::PadEdgeToward(int vertexBead, int halfEdge, int minInteriors)
{
    for (int guard = 0; guard < 8; ++guard) {
        const std::vector<int> chain = WalkSideChain(vertexBead, halfEdge);
        if (chain.size() < 2) {
            return false;
        }
        if (static_cast<int>(chain.size()) - 2 >= minInteriors) {
            return true;
        }
        if (InsertEdgeBeadBetween(chain[0], chain[1], std::abs(halfEdge)) < 0) {
            return false;
        }
        structureDirty_ = true;
    }
    return false;
}

int Physics::StrandOfLiveBead(const Bead& b) const
{
    if (b.removed || b.isVertex) {
        return -1;
    }
    if (b.edge == 0) {
        return b.strand;
    }
    return StrandOfEdgeLabel(b.edge);
}

int Physics::StrandOfEdgeLabel(int edge) const
{
    if (loop_ == nullptr || edge < 1 ||
        edge >= static_cast<int>(loop_->strandOfEdge.size())) {
        return -1;
    }
    return loop_->strandOfEdge[static_cast<size_t>(edge)];
}

int Physics::StrandLiveEdgeBeadCount(int strand) const
{
    if (strand < 0) {
        return 0;
    }
    int n = 0;
    for (int i = 0; i < BeadCount(); ++i) {
        if (StrandOfLiveBead(neck_->beads[static_cast<size_t>(i)]) == strand) {
            ++n;
        }
    }
    return n;
}

int Physics::PickR2ConserveEdge() const
{
    std::unordered_set<int> shortStrands;
    for (const int k : shortenPadEdges_) {
        const int s = StrandOfEdgeLabel(k);
        if (s >= 0 && StrandLiveEdgeBeadCount(s) < kMinStrandEdgeBeads) {
            shortStrands.insert(s);
        }
    }
    std::vector<int> candidates;
    for (const int k : shortenPadEdges_) {
        if (shortStrands.empty() ||
            shortStrands.count(StrandOfEdgeLabel(k)) != 0) {
            candidates.push_back(k);
        }
    }
    if (candidates.empty()) {
        candidates = shortenPadEdges_;
    }
    int bestEdge = -1;
    double bestTension = -1.0;
    for (const int k : candidates) {
        const double t = EdgeTension(k);
        if (t > bestTension) {
            bestTension = t;
            bestEdge = k;
        }
    }
    return bestEdge;
}

bool Physics::PadStrandsToMinBeads(int minBeads, const std::vector<int>& avoidEdges)
{
    const auto isAvoided = [&](int edge) {
        return std::find(avoidEdges.begin(), avoidEdges.end(), edge) != avoidEdges.end();
    };
    for (int guard = 0; guard < 64; ++guard) {
        if (structureDirty_) {
            RebuildStructures();
        }
        std::unordered_map<int, int> counts;
        for (int i = 0; i < BeadCount(); ++i) {
            const int s = StrandOfLiveBead(neck_->beads[static_cast<size_t>(i)]);
            if (s >= 0) {
                ++counts[s];
            }
        }
        int shortStrand = -1;
        for (const auto& entry : counts) {
            if (entry.second < minBeads) {
                shortStrand = entry.first;
                break;
            }
        }
        if (shortStrand < 0) {
            return true;
        }

        int bestAny = -1;
        int bestPrefer = -1;
        double tAny = -1.0;
        double tPrefer = -1.0;
        if (loop_ != nullptr) {
            for (int k = 1; k <= loop_->numEdges; ++k) {
                if (StrandOfEdgeLabel(k) != shortStrand) {
                    continue;
                }
                RefreshEdgeChain(k);
                const double t = EdgeTension(k);
                if (t > tAny) {
                    tAny = t;
                    bestAny = k;
                }
                if (!isAvoided(k) && t > tPrefer) {
                    tPrefer = t;
                    bestPrefer = k;
                }
            }
        }
        const int pick = (bestPrefer >= 0) ? bestPrefer : bestAny;
        if (pick >= 0) {
            if (InsertBeadAtFarthestGap(pick) < 0) {
                return false;
            }
            continue;
        }

        const std::vector<int> ring = RingBeadsOfCircleStrand(shortStrand);
        if (ring.size() < 2) {
            return false;
        }
        float bestDist = -1.0f;
        int bestA = -1;
        int bestC = -1;
        for (size_t i = 0; i < ring.size(); ++i) {
            const int a = ring[i];
            const int c = ring[(i + 1) % ring.size()];
            if (a < 0 || c < 0 || a >= BeadCount() || c >= BeadCount()) {
                continue;
            }
            if (neck_->beads[static_cast<size_t>(a)].removed ||
                neck_->beads[static_cast<size_t>(c)].removed) {
                continue;
            }
            const Vector2 pa = PosOf(a);
            const Vector2 pc = PosOf(c);
            const float d =
                std::sqrt((pc.x - pa.x) * (pc.x - pa.x) + (pc.y - pa.y) * (pc.y - pa.y));
            if (d > bestDist) {
                bestDist = d;
                bestA = a;
                bestC = c;
            }
        }
        if (bestA < 0 || InsertEdgeBeadBetween(bestA, bestC, 0) < 0) {
            return false;
        }
        RefreshCircleStructures();
        structureDirty_ = true;
    }
    return false;
}

bool Physics::StartR3Move(int face)
{
    if (unknot_ || MoveActive() || autoActive_ || FaceHasHoles(face) || FaceMergePending()) {
        return false;
    }
    ClearR3State();

    R3FaceInfo info;
    if (!AnalyzeR3Face(*loop_, nullptr, face, &info)) {
        return false;
    }
    for (int i = 0; i < 3; ++i) {
        if (info.crossing[i] < 0 ||
            info.crossing[i] >= static_cast<int>(neck_->crossingBead.size())) {
            return false;
        }
        info.vertexBead[i] = neck_->crossingBead[static_cast<size_t>(info.crossing[i])];
    }

    const auto padToward = [&](int vertexBead, int halfEdge, int minInteriors) -> bool {
        return PadEdgeToward(vertexBead, halfEdge, minInteriors);
    };

    // Pad any empty face-cycle side up to one interior bead.
    for (int i = 0; i < 3; ++i) {
        if (!padToward(info.vertexBead[i], info.h[i], 1)) {
            return false;
        }
    }

    // An exterior a/b edge that runs between two triangle vertices (a bigon, or
    // a loop at one vertex) needs two interior beads. With only one, both
    // vertices' NeighborToward results alias to the same terminal and R3 refuses.
    for (int i = 0; i < 3; ++i) {
        for (const int h : {info.a[i], info.b[i]}) {
            const std::vector<int> chain = WalkSideChain(info.vertexBead[i], h);
            int minInteriors = 1;
            if (chain.size() >= 2) {
                const int other = chain.back();
                for (int j = 0; j < 3; ++j) {
                    if (info.vertexBead[j] == other) {
                        minInteriors = 2;
                        break;
                    }
                }
            }
            if (!padToward(info.vertexBead[i], h, minInteriors)) {
                return false;
            }
        }
    }

    if (!AnalyzeR3Face(*loop_, neck_, face, &info)) {
        return false;
    }

    {
        std::vector<int> boundary;
        std::vector<float> dots;
        WalkFaceBoundaryFromPhiCycle(loop_->faces[static_cast<size_t>(face)], *neck_,
                                     crossingOf_, boundary, dots, nullptr, 0.0f, false);
        if (boundary.size() < 3) {
            return false;
        }
        std::vector<Vec2d> poly;
        poly.reserve(boundary.size());
        for (const int bead : boundary) {
            const Vector2 p = PosOf(bead);
            poly.push_back(Vec2d{static_cast<double>(p.x), static_cast<double>(p.y)});
        }
        info.ccw = SignedShoelaceArea(poly) > 0.0;
    }

    {
        Permutation newSigma;
        if (!ComputeR3Sigma(*loop_, info, &newSigma)) {
            return false;
        }
    }

    RebuildStructures();

    r3ExteriorEdges_.clear();
    for (int i = 0; i < 3; ++i) {
        for (const int h : {info.a[i], info.b[i]}) {
            const int k = std::abs(h);
            bool have = false;
            for (const int e : r3ExteriorEdges_) {
                if (e == k) {
                    have = true;
                    break;
                }
            }
            if (!have) {
                r3ExteriorEdges_.push_back(k);
            }
        }
    }

    r3Info_ = info;
    r3Active_ = true;
    r3Complete_ = false;
    r3Phase_ = R3Phase::PreHoneycomb;
    r3PreSub_ = R3PreSub::NeedRemove;
    structureDirty_ = true;
    return true;
}

void Physics::R3PreHoneycombStep()
{
    if (structureDirty_) {
        RebuildStructures();
    }

    if (r3PreSub_ == R3PreSub::NeedRemove) {
        bool allOne = true;
        std::vector<int> oversized;
        for (int i = 0; i < 3; ++i) {
            const int k = std::abs(r3Info_.h[i]);
            RefreshEdgeChain(k);
            const int n = EdgeInteriorCount(k);
            if (n < 1) {
                const std::vector<int>& chain = chains_[k];
                if (chain.size() >= 2) {
                    InsertEdgeBeadBetween(chain[0], chain[1], k);
                    RefreshEdgeChain(k);
                    RebuildSpringsFromChains();
                    RefreshFaceBodiesForEdge(k);
                    RebuildFaceBodiesFromCache();
                }
            }
            else if (n > 1) {
                allOne = false;
                oversized.push_back(k);
            }
        }
        if (allOne) {
            BeginFortifyOrSlide();
            return;
        }
        if (oversized.empty()) {
            CancelR3Move();
            return;
        }

        std::vector<int> candidates;
        for (const int k : oversized) {
            if (!NeighborHoneycombBlocksRemoval(k)) {
                candidates.push_back(k);
            }
        }
        if (candidates.empty()) {
            // Cannot strip further without breaking neighbor honeycombs.
            if (!BeginR3Slide()) {
                CancelR3Move();
            }
            return;
        }

        const auto phiIndex = [this](int edge) {
            for (int i = 0; i < 3; ++i) {
                if (std::abs(r3Info_.h[i]) == edge) {
                    return i;
                }
            }
            return 3;
        };
        int pick = candidates[0];
        float bestLen = EdgeChainLength(pick);
        int bestPhi = phiIndex(pick);
        for (size_t c = 1; c < candidates.size(); ++c) {
            const int k = candidates[c];
            const float len = EdgeChainLength(k);
            const int phi = phiIndex(k);
            if (len > bestLen || (len == bestLen && phi < bestPhi)) {
                pick = k;
                bestLen = len;
                bestPhi = phi;
            }
        }
        r3SettleLength_ = EdgeChainLength(pick);
        if (RemoveFarthestFromFaceCenter(pick, r3Info_.face) < 0) {
            CancelR3Move();
            return;
        }
        r3SettleEdge_ = pick;
        r3SettleFrames_ = 0;
        r3PreSub_ = R3PreSub::SettlingRemove;
        return;
    }

    if (r3PreSub_ == R3PreSub::SettlingRemove) {
        if (R3SettleTick(r3SettleEdge_, r3SettleLength_)) {
            r3PreSub_ = R3PreSub::NeedInsert;
        }
        return;
    }

    if (r3PreSub_ == R3PreSub::NeedInsert) {
        if (r3ExteriorEdges_.empty()) {
            r3PreSub_ = R3PreSub::NeedRemove;
            return;
        }
        for (const int k : r3ExteriorEdges_) {
            RefreshEdgeChain(k);
        }
        int bestEdge = -1;
        double bestTension = -1.0;
        for (const int k : r3ExteriorEdges_) {
            const double t = EdgeTension(k);
            if (t > bestTension) {
                bestTension = t;
                bestEdge = k;
            }
        }
        if (bestEdge < 0) {
            r3PreSub_ = R3PreSub::NeedRemove;
            return;
        }
        const float length = EdgeChainLength(bestEdge);
        const double threshold =
            params_->tensionThreshold * std::sqrt(static_cast<double>(std::max(length, 1.0f)));
        if (bestTension < threshold) {
            r3PreSub_ = R3PreSub::NeedRemove;
            return;
        }
        r3SettleLength_ = length;
        if (InsertBeadAtFarthestGap(bestEdge) < 0) {
            r3PreSub_ = R3PreSub::NeedRemove;
            return;
        }
        r3SettleEdge_ = bestEdge;
        r3SettleFrames_ = 0;
        r3PreSub_ = R3PreSub::SettlingInsert;
        return;
    }

    if (R3SettleTick(r3SettleEdge_, r3SettleLength_)) {
        r3PreSub_ = R3PreSub::NeedRemove;
    }
}

bool Physics::R3HoneycombSpokesActive() const
{
    return HoneycombSlotForFace(r3Info_.face) >= 0;
}

bool Physics::BindR3HoneycombSpokes()
{
    if (structureDirty_) {
        RebuildStructures();
    }
    return HoneycombSlotForFace(r3Info_.face) >= 0;
}

bool Physics::R3HoneycombSpokesSettled() const
{
    const int slot = HoneycombSlotForFace(r3Info_.face);
    if (slot < 0 || honeycombSpokes_.size() < static_cast<size_t>(slot + 1) * 6) {
        return false;
    }
    const int fi = honeycombFaceFi_[static_cast<size_t>(slot)];
    const FaceBead& fb = neck_->faceBeads[static_cast<size_t>(fi)];
    if (!FaceBeadInsideFace(fb)) {
        return false;
    }
    std::vector<Vec2d> poly;
    if (FaceBoundaryPolygon(fb.face, &poly)) {
        const double precision = std::max(1e-4, 0.05 * static_cast<double>(radius_));
        const InscribedCircle circle = LargestInscribedCircle(poly, precision);
        if (!circle.found) {
            return false;
        }
        const float dx = static_cast<float>(circle.center.x) - fb.pos.x;
        const float dy = static_cast<float>(circle.center.y) - fb.pos.y;
        if (std::sqrt(dx * dx + dy * dy) > 0.2f * diameter_) {
            return false;
        }
    }
    const int faceBody = BeadCount() + fi;
    const float band = 0.35f * diameter_;
    const Vector2 pa = PosOf(faceBody);
    for (int i = 0; i < 6; ++i) {
        const int spoke =
            honeycombSpokes_[static_cast<size_t>(slot) * 6 + static_cast<size_t>(i)];
        if (spoke < 0 || spoke >= BeadCount() ||
            neck_->beads[static_cast<size_t>(spoke)].removed) {
            return false;
        }
        const Vector2 pb = PosOf(spoke);
        const float dx = pb.x - pa.x;
        const float dy = pb.y - pa.y;
        const float len = std::sqrt(dx * dx + dy * dy);
        if (std::fabs(len - diameter_) > band) {
            return false;
        }
    }
    return true;
}

void Physics::BeginFortifyOrSlide()
{
    if (BindR3HoneycombSpokes()) {
        r3Phase_ = R3Phase::Fortify;
        r3StallFrames_ = 0;
        r3SettleFrames_ = 0;
        return;
    }
    if (!BeginR3Slide()) {
        CancelR3Move();
    }
}

void Physics::R3FortifyStep()
{
    for (int i = 0; i < params_->autoStepsPerFrame; ++i) {
        Step();
    }
    r3SettleFrames_ += params_->autoStepsPerFrame;
    if (R3HoneycombSpokesSettled() || r3SettleFrames_ >= params_->relaxStepsMax) {
        if (!BeginR3Slide()) {
            CancelR3Move();
        }
    }
}

bool Physics::BeginR3Slide()
{
    // Snapshot complete face collision groups from the intact necklace.
    // Cuts below break face walks, so we must not rebuild faceBodies_ after.
    RebuildStructures();

    for (int i = 0; i < 3; ++i) {
        r3Info_.aTerminal[i] =
            neck_->NeighborToward(r3Info_.vertexBead[i], r3Info_.a[i]);
        r3Info_.bTerminal[i] =
            neck_->NeighborToward(r3Info_.vertexBead[i], r3Info_.b[i]);
        if (r3Info_.aTerminal[i] < 0 || r3Info_.bTerminal[i] < 0 ||
            neck_->beads[static_cast<size_t>(r3Info_.aTerminal[i])].isVertex ||
            neck_->beads[static_cast<size_t>(r3Info_.bTerminal[i])].isVertex) {
            return false;
        }
    }

    for (int i = 0; i < 3; ++i) {
        const int j = (i + 1) % 3;
        const std::vector<int> chain =
            WalkSideChain(r3Info_.vertexBead[i], r3Info_.h[i]);
        if (static_cast<int>(chain.size()) < 3 ||
            chain.back() != r3Info_.vertexBead[j]) {
            return false;
        }
        const int interiors = static_cast<int>(chain.size()) - 2;
        const int mid = (interiors - 1) / 2;
        const int meetIdx = 1 + mid;
        if (neck_->beads[static_cast<size_t>(chain[static_cast<size_t>(meetIdx)])].isVertex) {
            return false;
        }
        r3Meeting_[i] = chain[static_cast<size_t>(meetIdx)];

        R3Slider& f = r3Sliders_[2 * i];
        f = R3Slider{};
        f.terminal = r3Info_.bTerminal[i];
        f.cutVertex = r3Info_.vertexBead[i];
        f.path.assign(chain.begin(), chain.begin() + meetIdx + 1);
        f.targetIdx = 0;

        R3Slider& g = r3Sliders_[2 * i + 1];
        g = R3Slider{};
        g.terminal = r3Info_.aTerminal[j];
        g.cutVertex = r3Info_.vertexBead[j];
        g.path.clear();
        for (int idx = static_cast<int>(chain.size()) - 1; idx >= meetIdx; --idx) {
            g.path.push_back(chain[static_cast<size_t>(idx)]);
        }
        g.targetIdx = 0;
    }

    const auto cutLink = [&](R3Slider& s, int halfEdge) {
        Bead& v = neck_->beads[static_cast<size_t>(s.cutVertex)];
        s.cutVertexSlot = -1;
        for (size_t k = 0; k < v.sigmaCycle.size(); ++k) {
            if (v.sigmaCycle[k] == halfEdge) {
                s.cutVertexSlot = static_cast<int>(k);
                break;
            }
        }
        if (s.cutVertexSlot < 0) {
            return false;
        }
        v.neighborOfHalfEdge[s.cutVertexSlot] = -1;
        Bead& t = neck_->beads[static_cast<size_t>(s.terminal)];
        s.cutNeighborSlot = (t.neighbors[0] == s.cutVertex) ? 0
                            : ((t.neighbors[1] == s.cutVertex) ? 1 : -1);
        if (s.cutNeighborSlot < 0) {
            return false;
        }
        t.neighbors[s.cutNeighborSlot] = -1;
        return true;
    };

    r3CutEdges_.clear();
    bool cutsOk = true;
    for (int i = 0; i < 3 && cutsOk; ++i) {
        const int j = (i + 1) % 3;
        cutsOk = cutLink(r3Sliders_[2 * i], r3Info_.b[i]) &&
                 cutLink(r3Sliders_[2 * i + 1], r3Info_.a[j]);
    }
    for (int i = 0; i < 3; ++i) {
        const int j = (i + 1) % 3;
        const int bEdge = std::abs(r3Info_.b[i]);
        const int aEdge = std::abs(r3Info_.a[j]);
        bool haveB = false;
        bool haveA = false;
        for (const auto& entry : r3CutEdges_) {
            haveB = haveB || entry.first == bEdge;
            haveA = haveA || entry.first == aEdge;
        }
        if (!haveB) {
            r3CutEdges_.push_back({bEdge, r3Info_.bTerminal[i]});
        }
        if (!haveA) {
            r3CutEdges_.push_back({aEdge, r3Info_.aTerminal[j]});
        }
    }

    if (!cutsOk) {
        // Restore any partial cuts.
        for (R3Slider& s : r3Sliders_) {
            if (s.terminal < 0) {
                continue;
            }
            if (s.cutNeighborSlot >= 0) {
                neck_->beads[static_cast<size_t>(s.terminal)].neighbors[s.cutNeighborSlot] =
                    s.cutVertex;
            }
            if (s.cutVertexSlot >= 0) {
                neck_->beads[static_cast<size_t>(s.cutVertex)]
                    .neighborOfHalfEdge[s.cutVertexSlot] = s.terminal;
            }
            s = R3Slider{};
        }
        r3CutEdges_.clear();
        return false;
    }

    // Springs/chains follow the cut terminals; collision groups stay pre-cut.
    RebuildStructuresR3();
    r3Phase_ = R3Phase::Slide;
    r3StallFrames_ = 0;
    return true;
}

void Physics::CancelR3Move()
{
    if (!r3Active_) {
        return;
    }
    if (r3Phase_ == R3Phase::Slide) {
        for (R3Slider& s : r3Sliders_) {
            if (s.terminal < 0) {
                continue;
            }
            if (s.cutNeighborSlot >= 0 && s.terminal < BeadCount() &&
                !neck_->beads[static_cast<size_t>(s.terminal)].removed) {
                neck_->beads[static_cast<size_t>(s.terminal)].neighbors[s.cutNeighborSlot] =
                    s.cutVertex;
            }
            if (s.cutVertexSlot >= 0 && s.cutVertex < BeadCount()) {
                neck_->beads[static_cast<size_t>(s.cutVertex)]
                    .neighborOfHalfEdge[s.cutVertexSlot] = s.terminal;
            }
        }
    }
    r3Active_ = false;
    r3Complete_ = false;
    ClearR3State();
    structureDirty_ = true;
}

void Physics::GetR3Result(R3FaceInfo* info, int meeting[3]) const
{
    if (info != nullptr) {
        *info = r3Info_;
    }
    if (meeting != nullptr) {
        for (int i = 0; i < 3; ++i) {
            meeting[i] = r3Meeting_[i];
        }
    }
}

std::vector<int> Physics::WalkFromTerminal(int terminal) const
{
    std::vector<int> chain;
    chain.push_back(terminal);
    const Bead& t = neck_->beads[static_cast<size_t>(terminal)];
    int prev = terminal;
    int cur = (t.neighbors[0] >= 0) ? t.neighbors[0] : t.neighbors[1];
    int guard = 0;
    const int limit = BeadCount() + 8;
    while (cur >= 0 && !neck_->beads[static_cast<size_t>(cur)].isVertex) {
        chain.push_back(cur);
        const Bead& b = neck_->beads[static_cast<size_t>(cur)];
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

void Physics::RebuildStructuresR3()
{
    // Chains and springs only. Face collision groups stay frozen from the
    // pre-cut snapshot; rebuilding them here would drop the cut a/b edges.
    const int n = loop_->numEdges;
    chains_.assign(n + 1, {});
    for (int k = 1; k <= n; ++k) {
        chains_[k] = WalkEdge(k);
    }
    for (const auto& entry : r3CutEdges_) {
        chains_[entry.first] = WalkFromTerminal(entry.second);
    }
    RebuildSpringsFromChains();
    structureDirty_ = false;
}

void Physics::AddBeadToFrozenGroups(int bead, int adjacent)
{
    for (std::vector<int>& group : faceBodies_) {
        if (std::find(group.begin(), group.end(), adjacent) != group.end()) {
            group.push_back(bead);
        }
    }
}

void Physics::MaintainR3BranchSpacing()
{
    if (structureDirty_) {
        RebuildStructuresR3();
    }
    const float gapThreshold = params_->gapSize * diameter_;
    for (const auto& entry : r3CutEdges_) {
        const std::vector<int> chain = chains_[entry.first];
        for (size_t i = 0; i + 1 < chain.size(); ++i) {
            const Vector2 pa = PosOf(chain[i]);
            const Vector2 pb = PosOf(chain[i + 1]);
            const float dx = pb.x - pa.x;
            const float dy = pb.y - pa.y;
            if (std::sqrt(dx * dx + dy * dy) <= gapThreshold) {
                continue;
            }
            const int beadsBefore = BeadCount();
            const int inserted = InsertEdgeBeadBetween(chain[i], chain[i + 1], entry.first);
            if (inserted >= 0) {
                if (BeadCount() > beadsBefore) {
                    const int delta = BeadCount() - beadsBefore;
                    for (std::vector<int>& group : faceBodies_) {
                        for (int& body : group) {
                            if (body >= beadsBefore) {
                                body += delta;
                            }
                        }
                    }
                }
                AddBeadToFrozenGroups(inserted, chain[i]);
            }
            break;
        }
    }
}

bool Physics::R3SliderSettled(const R3Slider& s) const
{
    if (s.path.empty() || s.targetIdx < 0 ||
        s.targetIdx >= static_cast<int>(s.path.size())) {
        return false;
    }
    const Vector2 pa = PosOf(s.terminal);
    const Vector2 pb = PosOf(s.path[static_cast<size_t>(s.targetIdx)]);
    const float dx = pb.x - pa.x;
    const float dy = pb.y - pa.y;
    return std::sqrt(dx * dx + dy * dy) < 1.35f * diameter_;
}

void Physics::AdvanceR3Slider(R3Slider& s)
{
    if (s.targetIdx + 1 >= static_cast<int>(s.path.size())) {
        return;
    }
    ++s.targetIdx;
    const Vector2 from = PosOf(s.path[static_cast<size_t>(s.targetIdx) - 1]);
    const Vector2 to = PosOf(s.path[static_cast<size_t>(s.targetIdx)]);
    const float dx = to.x - from.x;
    const float dy = to.y - from.y;
    const float len = std::sqrt(dx * dx + dy * dy);
    if (len < 1e-6f) {
        return;
    }
    const float boost = 0.5f * params_->maxSpeed;
    vel_[s.terminal].x += boost * dx / len;
    vel_[s.terminal].y += boost * dy / len;
}

void Physics::R3Step()
{
    if (!r3Active_ || r3Complete_) {
        return;
    }
    if (r3Phase_ == R3Phase::PreHoneycomb) {
        R3PreHoneycombStep();
        return;
    }
    if (r3Phase_ == R3Phase::Fortify) {
        R3FortifyStep();
        return;
    }
    if (r3Phase_ == R3Phase::PostRepair) {
        R3PostRepairStep();
        return;
    }

    for (int i = 0; i < params_->autoStepsPerFrame; ++i) {
        Step();
    }
    MaintainR3BranchSpacing();

    bool advanced = false;
    bool allConverged = true;
    for (int i = 0; i < 3; ++i) {
        R3Slider& f = r3Sliders_[2 * i];
        R3Slider& g = r3Sliders_[2 * i + 1];
        const int m = static_cast<int>(f.path.size()) - 1;
        // Precomputed meeting: both paths end at the same honeycomb bead.
        if (f.targetIdx + g.targetIdx >= m + 1) {
            continue;
        }
        allConverged = false;
        if (!R3SliderSettled(f) || !R3SliderSettled(g)) {
            continue;
        }
        const int sum = f.targetIdx + g.targetIdx;
        if (sum <= m - 1) {
            AdvanceR3Slider(f);
            AdvanceR3Slider(g);
        }
        else {
            AdvanceR3Slider(f.targetIdx <= g.targetIdx ? f : g);
        }
        advanced = true;
    }

    if (allConverged) {
        bool settled = true;
        for (const R3Slider& s : r3Sliders_) {
            if (!R3SliderSettled(s)) {
                settled = false;
                break;
            }
        }
        if (settled) {
            // Prefer the precomputed meetings from the honeycomb phase.
            for (int i = 0; i < 3; ++i) {
                if (r3Meeting_[i] < 0) {
                    const R3Slider& f = r3Sliders_[2 * i];
                    r3Meeting_[i] = f.path[static_cast<size_t>(f.targetIdx)];
                }
            }
            r3Complete_ = true;
            return;
        }
    }

    if (advanced) {
        r3StallFrames_ = 0;
    }
    else if (++r3StallFrames_ > 240) {
        CancelR3Move();
    }
}

void Physics::EnsureR3ExteriorHoneycombMins(const Multiloop& newLoop)
{
    // Exterior labels from the pre-move analysis (half-edge labels unchanged).
    std::vector<int> edges;
    for (int i = 0; i < 3; ++i) {
        edges.push_back(r3Info_.a[i]);
        edges.push_back(r3Info_.b[i]);
    }
    EnsureExteriorHoneycombMins(newLoop, edges);
}

void Physics::EnsureExteriorHoneycombMins(const Multiloop& newLoop, const std::vector<int>& edges)
{
    if (newLoop.numEdges <= 0 || newLoop.sigma.MaxAbsLabel() == 0) {
        return;
    }
    // Build half-edge → crossing and half-edge → face maps for the rewritten loop.
    std::unordered_map<int, int> crossingOf;
    const std::vector<std::vector<int>> cycles = newLoop.sigma.Cycles();
    for (int c = 0; c < static_cast<int>(cycles.size()); ++c) {
        for (const int h : cycles[static_cast<size_t>(c)]) {
            crossingOf[h] = c;
        }
    }
    std::unordered_map<int, int> faceOf;
    for (int f = 0; f < static_cast<int>(newLoop.faces.size()); ++f) {
        for (const int h : newLoop.faces[static_cast<size_t>(f)]) {
            faceOf[h] = f;
        }
    }

    const auto walkEdgeNew = [&](int edge) {
        std::vector<int> chain;
        const auto itStart = crossingOf.find(edge);
        if (itStart == crossingOf.end()) {
            return chain;
        }
        if (itStart->second < 0 ||
            itStart->second >= static_cast<int>(neck_->crossingBead.size())) {
            return chain;
        }
        const int vStart = neck_->crossingBead[static_cast<size_t>(itStart->second)];
        if (vStart < 0 || vStart >= BeadCount()) {
            return chain;
        }
        chain.push_back(vStart);
        int prev = vStart;
        int cur = neck_->NeighborToward(vStart, edge);
        int guard = 0;
        const int limit = BeadCount() + 8;
        while (cur >= 0 && !neck_->beads[static_cast<size_t>(cur)].isVertex) {
            if (neck_->beads[static_cast<size_t>(cur)].removed) {
                break;
            }
            chain.push_back(cur);
            const Bead& b = neck_->beads[static_cast<size_t>(cur)];
            const int nxt = (b.neighbors[0] == prev) ? b.neighbors[1] : b.neighbors[0];
            if (nxt == prev) {
                break;
            }
            prev = cur;
            cur = nxt;
            if (++guard > limit) {
                break;
            }
        }
        if (cur >= 0 && neck_->beads[static_cast<size_t>(cur)].isVertex) {
            chain.push_back(cur);
        }
        return chain;
    };

    const auto isLoop = [&](int edge) {
        const auto itPos = crossingOf.find(edge);
        const auto itNeg = crossingOf.find(-edge);
        return itPos != crossingOf.end() && itNeg != crossingOf.end() &&
               itPos->second == itNeg->second;
    };

    const auto faceSize = [&](int halfEdge) {
        const auto it = faceOf.find(halfEdge);
        if (it == faceOf.end()) {
            return 0;
        }
        return static_cast<int>(newLoop.faces[static_cast<size_t>(it->second)].size());
    };

    const auto insertFarthest = [&](int edge) {
        const std::vector<int> chain = walkEdgeNew(edge);
        if (chain.size() < 2) {
            return -1;
        }
        size_t best = 0;
        float bestDist = -1.0f;
        for (size_t i = 0; i + 1 < chain.size(); ++i) {
            const Vector2 a = PosOf(chain[i]);
            const Vector2 b = PosOf(chain[i + 1]);
            const float d =
                std::sqrt((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y));
            if (d > bestDist) {
                bestDist = d;
                best = i;
            }
        }
        return InsertEdgeBeadBetween(chain[best], chain[best + 1], edge);
    };

    std::vector<int> exteriors;
    for (const int h : edges) {
        const int k = std::abs(h);
        if (k < 1 || k > newLoop.numEdges) {
            continue;
        }
        if (std::find(exteriors.begin(), exteriors.end(), k) == exteriors.end()) {
            exteriors.push_back(k);
        }
    }

    for (const int k : exteriors) {
        int minBeads = 0;
        if (isLoop(k)) {
            minBeads = 5; // monogon
        }
        else if (faceSize(k) == 2 || faceSize(-k) == 2) {
            minBeads = 2; // bigon
        }
        if (minBeads <= 0) {
            continue;
        }
        for (int guard = 0; guard < 16; ++guard) {
            const std::vector<int> chain = walkEdgeNew(k);
            const int interiors =
                (chain.size() >= 2) ? static_cast<int>(chain.size()) - 2 : 0;
            if (interiors >= minBeads) {
                break;
            }
            if (insertFarthest(k) < 0) {
                break;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Shortening moves (R1 / R2)
// ---------------------------------------------------------------------------

void Physics::MoveStep()
{
    if (r0Active_) {
        R0Step();
        return;
    }
    if (r1Active_ || r2Active_) {
        ShortenStep();
        return;
    }
    if (r3Active_) {
        R3Step();
    }
}

void Physics::CancelMove()
{
    if (r0Active_) {
        CancelR0Move();
        return;
    }
    if (r1Active_ || r2Active_) {
        CancelShorten();
        return;
    }
    CancelR3Move();
}

void Physics::CancelShorten()
{
    r1Active_ = false;
    r1Complete_ = false;
    r2Active_ = false;
    r2Complete_ = false;
    shortenStripEdges_.clear();
    shortenStripTarget_.clear();
    shortenPadEdges_.clear();
    shortenStallFrames_ = 0;
    spokeFrames_ = 0;
    shortenSub_ = ShortenSub::NeedRemove;
    ClearR3State();
    RestoreErasedFaceBead();
    structureDirty_ = true;
}

void Physics::GetR1Result(R1FaceInfo* info) const
{
    if (info != nullptr) {
        *info = r1Info_;
    }
}

void Physics::GetR2Result(R2FaceInfo* info) const
{
    if (info != nullptr) {
        *info = r2Info_;
    }
}

void Physics::GetR0Result(R0MoveInfo* info) const
{
    if (info == nullptr) {
        return;
    }
    info->strand = r0Strand_;
    info->circleOrdinal = r0CircleOrdinal_;
    info->ringBeads = r0Pad_.boundaryBeads;
}

bool Physics::CircleSettleTick()
{
    for (int i = 0; i < params_->autoStepsPerFrame; ++i) {
        Step();
    }
    r3SettleFrames_ += params_->autoStepsPerFrame;
    if (structureDirty_) {
        RebuildStructures();
    }
    RefreshCircleBoundary(r0Pad_);
    const double threshold =
        params_->tensionThreshold *
        std::sqrt(static_cast<double>(std::max(r0SettleLength_, 1.0f)));
    if (CircleRingTension(r0Pad_) < threshold) {
        return true;
    }
    if (r3SettleFrames_ >= params_->relaxStepsMax) {
        return true;
    }
    return false;
}

void Physics::CancelR0Move()
{
    r0Active_ = false;
    r0Complete_ = false;
    r0Pad_ = FaceBead{};
    r0Strand_ = -1;
    r0CircleOrdinal_ = -1;
    r0Sub_ = ShortenSub::NeedRemove;
    r0StallFrames_ = 0;
    r0SettleLength_ = 0.0f;
    ClearR3State();
    RestoreErasedFaceBead();
    structureDirty_ = true;
}

bool Physics::StartR0Move(int faceBeadIndex)
{
    if (MoveActive() || autoActive_ || FaceMergePending()) {
        return false;
    }
    if (faceBeadIndex < 0 ||
        faceBeadIndex >= static_cast<int>(neck_->faceBeads.size())) {
        return false;
    }
    const FaceBead& fb = neck_->faceBeads[static_cast<size_t>(faceBeadIndex)];
    if (!IsCircleFaceBead(fb) || !fb.holeBeads.empty()) {
        return false;
    }
    const int strand = PadStrand(fb);
    if (CircleHasHoles(strand)) {
        return false;
    }
    for (const int bead : fb.boundaryBeads) {
        if (bead < 0 || bead >= BeadCount()) {
            continue;
        }
        if (neck_->beads[static_cast<size_t>(bead)].isVertex) {
            return false;
        }
    }
    if (LiveBoundaryCount(fb) < 3) {
        return false;
    }

    int ordinal = 0;
    for (int i = 0; i < faceBeadIndex; ++i) {
        if (IsCircleFaceBead(neck_->faceBeads[static_cast<size_t>(i)])) {
            ++ordinal;
        }
    }

    ClearR3State();
    r0Pad_ = fb;
    r0Strand_ = strand;
    r0CircleOrdinal_ = ordinal;
    r0Active_ = true;
    r0Complete_ = false;
    r0Sub_ = ShortenSub::NeedRemove;
    r0StallFrames_ = 0;
    r0SettleLength_ = 0.0f;
    if (!EraseFaceBeadAt(static_cast<size_t>(faceBeadIndex))) {
        r0Active_ = false;
        r0Pad_ = FaceBead{};
        return false;
    }
    RefreshCircleBoundary(r0Pad_);
    structureDirty_ = true;
    RebuildStructures();
    return true;
}

void Physics::R0Step()
{
    if (!r0Active_ || r0Complete_) {
        return;
    }
    if (++r0StallFrames_ > 4000) {
        CancelR0Move();
        return;
    }
    if (structureDirty_) {
        RebuildStructures();
    }

    if (r0Sub_ == ShortenSub::NeedRemove) {
        RefreshCircleBoundary(r0Pad_);
        if (LiveBoundaryCount(r0Pad_) <= 3) {
            r0Complete_ = true;
            return;
        }
        const int bead = FarthestLiveBoundaryBead(r0Pad_);
        if (bead < 0) {
            CancelR0Move();
            return;
        }
        r0SettleLength_ = CircleRingLength(r0Pad_);
        RemoveEdgeBead(bead);
        RefreshCircleBoundary(r0Pad_);
        r3SettleFrames_ = 0;
        r0Sub_ = ShortenSub::SettlingRemove;
        structureDirty_ = true;
        return;
    }

    if (r0Sub_ == ShortenSub::SettlingRemove) {
        if (CircleSettleTick()) {
            r0Sub_ = ShortenSub::NeedRemove;
        }
    }
}

Vector2 Physics::ShortenCenter() const
{
    if (r1Active_) {
        if (r1Info_.vertexBead >= 0 && r1Info_.vertexBead < BeadCount()) {
            return PosOf(r1Info_.vertexBead);
        }
    }
    else if (r2Active_) {
        const int x = r2Info_.vertexBead[0];
        const int y = r2Info_.vertexBead[1];
        if (x >= 0 && x < BeadCount() && y >= 0 && y < BeadCount()) {
            const Vector2 px = PosOf(x);
            const Vector2 py = PosOf(y);
            return Vector2{0.5f * (px.x + py.x), 0.5f * (px.y + py.y)};
        }
    }
    if (erasedFaceBeadValid_) {
        return Vector2{erasedFaceBead_.pos.x, erasedFaceBead_.pos.y};
    }
    return Vector2{0.0f, 0.0f};
}

bool Physics::StartR1Move(int face)
{
    if (unknot_ || MoveActive() || autoActive_ || FaceHasHoles(face) || FaceMergePending()) {
        return false;
    }
    R1FaceInfo info;
    if (!AnalyzeR1Face(*loop_, nullptr, face, &info)) {
        return false;
    }
    if (info.crossing < 0 || info.crossing >= static_cast<int>(neck_->crossingBead.size())) {
        return false;
    }
    const int vb = neck_->crossingBead[static_cast<size_t>(info.crossing)];
    if (vb < 0 || vb >= BeadCount()) {
        return false;
    }
    // The loop keeps two beads; p and q need distinct terminals (two beads
    // when they are the same edge, i.e. the one-crossing kink diagram).
    if (!PadEdgeToward(vb, info.e, 2)) {
        return false;
    }
    const int minPQ = (std::abs(info.p) == std::abs(info.q)) ? 2 : 1;
    if (!PadEdgeToward(vb, info.p, minPQ) || !PadEdgeToward(vb, info.q, minPQ)) {
        return false;
    }
    if (!AnalyzeR1Face(*loop_, neck_, face, &info)) {
        return false;
    }
    ClearR3State();
    r1Info_ = info;
    r1Active_ = true;
    r1Complete_ = false;
    shortenSub_ = ShortenSub::NeedRemove;
    shortenStripEdges_ = {std::abs(info.e)};
    shortenStripTarget_ = {2};
    shortenPadEdges_.clear();
    shortenStallFrames_ = 0;
    EraseFaceBead(face);
    structureDirty_ = true;
    RebuildStructures();
    return true;
}

bool Physics::StartR2Move(int face)
{
    if (unknot_ || MoveActive() || autoActive_ || FaceHasHoles(face) || FaceMergePending()) {
        return false;
    }
    R2FaceInfo info;
    if (!AnalyzeR2Face(*loop_, nullptr, face, &info)) {
        return false;
    }
    for (int i = 0; i < 2; ++i) {
        if (info.crossing[i] < 0 ||
            info.crossing[i] >= static_cast<int>(neck_->crossingBead.size())) {
            return false;
        }
        info.vertexBead[i] = neck_->crossingBead[static_cast<size_t>(info.crossing[i])];
        if (info.vertexBead[i] < 0 || info.vertexBead[i] >= BeadCount()) {
            return false;
        }
    }
    // Each bigon side keeps one bead.
    if (!PadEdgeToward(info.vertexBead[0], info.h0, 1) ||
        !PadEdgeToward(info.vertexBead[1], info.h1, 1)) {
        return false;
    }
    // Exterior arcs need distinct terminals: two beads when an arc runs
    // between X and Y (or loops back to its own vertex).
    const auto padExterior = [&](int vertexBead, int h) -> bool {
        const std::vector<int> chain = WalkSideChain(vertexBead, h);
        int minInteriors = 1;
        if (chain.size() >= 2) {
            const int other = chain.back();
            if (other == info.vertexBead[0] || other == info.vertexBead[1]) {
                minInteriors = 2;
            }
        }
        return PadEdgeToward(vertexBead, h, minInteriors);
    };
    if (!padExterior(info.vertexBead[0], info.a) || !padExterior(info.vertexBead[0], info.b) ||
        !padExterior(info.vertexBead[1], info.ap) || !padExterior(info.vertexBead[1], info.bp)) {
        return false;
    }
    if (info.disconnects && !IsExterior(info.faceL)) {
        int ext = -1;
        for (const int f : exteriorFaceOrder_) {
            bool hit = false;
            for (const int h : loop_->faces[static_cast<size_t>(f)]) {
                const int k = std::abs(h);
                if (k == std::abs(info.a) || k == std::abs(info.b) || k == std::abs(info.ap) ||
                    k == std::abs(info.bp)) {
                    hit = true;
                    break;
                }
            }
            if (hit) {
                ext = f;
                break;
            }
        }
        if (ext < 0 && !exteriorFaceOrder_.empty()) {
            ext = exteriorFaceOrder_[0];
        }
        const int opp = R2OpposingHostHalfEdge(*loop_, info, ext);
        if (opp != 0) {
            const int vb = (opp == info.a || opp == info.b) ? info.vertexBead[0]
                                                            : info.vertexBead[1];
            if (!PadEdgeToward(vb, opp, 3)) {
                return false;
            }
        }
    }
    {
        const std::vector<int> strip = {std::abs(info.h0), std::abs(info.h1)};
        PadStrandsToMinBeads(kMinStrandEdgeBeads, strip);
    }
    if (!AnalyzeR2Face(*loop_, neck_, face, &info)) {
        return false;
    }
    ClearR3State();
    r2Info_ = info;
    r2Active_ = true;
    r2Complete_ = false;
    shortenSub_ = ShortenSub::NeedRemove;
    shortenStripEdges_ = {std::abs(info.h0), std::abs(info.h1)};
    shortenStripTarget_ = {1, 1};
    shortenPadEdges_.clear();
    for (const int h : {info.a, info.b, info.ap, info.bp}) {
        const int k = std::abs(h);
        if (std::find(shortenStripEdges_.begin(), shortenStripEdges_.end(), k) !=
                shortenStripEdges_.end() ||
            std::find(shortenPadEdges_.begin(), shortenPadEdges_.end(), k) !=
                shortenPadEdges_.end()) {
            continue;
        }
        shortenPadEdges_.push_back(k);
    }
    shortenStallFrames_ = 0;
    spokeFrames_ = 0;
    EraseFaceBead(face);
    structureDirty_ = true;
    RebuildStructures();
    return true;
}

void Physics::ShortenStep()
{
    if ((!r1Active_ && !r2Active_) || r1Complete_ || r2Complete_) {
        return;
    }
    if (++shortenStallFrames_ > 4000) {
        CancelShorten();
        return;
    }
    if (structureDirty_) {
        RebuildStructures();
    }

    const auto interiorsOf = [&](int edge) {
        RefreshEdgeChain(edge);
        return EdgeInteriorCount(edge);
    };

    if (shortenSub_ == ShortenSub::NeedRemove) {
        int pick = -1;
        float bestLen = -1.0f;
        for (size_t i = 0; i < shortenStripEdges_.size(); ++i) {
            const int k = shortenStripEdges_[i];
            if (interiorsOf(k) <= shortenStripTarget_[i]) {
                continue;
            }
            const float len = EdgeChainLength(k);
            if (len > bestLen) {
                bestLen = len;
                pick = k;
            }
        }
        if (pick < 0) {
            if (r1Active_) {
                r1Complete_ = true;
                return;
            }
            shortenSub_ = ShortenSub::Spoke;
            spokeFrames_ = 0;
            return;
        }
        r3SettleLength_ = EdgeChainLength(pick);
        if (RemoveFarthestFromPoint(pick, ShortenCenter()) < 0) {
            CancelShorten();
            return;
        }
        r3SettleEdge_ = pick;
        r3SettleFrames_ = 0;
        shortenSub_ = ShortenSub::SettlingRemove;
        return;
    }

    if (shortenSub_ == ShortenSub::SettlingRemove) {
        if (R3SettleTick(r3SettleEdge_, r3SettleLength_)) {
            shortenSub_ = shortenPadEdges_.empty() ? ShortenSub::NeedRemove
                                                   : ShortenSub::NeedInsert;
        }
        return;
    }

    if (shortenSub_ == ShortenSub::NeedInsert) {
        for (const int k : shortenPadEdges_) {
            RefreshEdgeChain(k);
        }
        int bestEdge = -1;
        double bestTension = -1.0;
        if (r2Active_) {
            bestEdge = PickR2ConserveEdge();
            if (bestEdge >= 0) {
                bestTension = EdgeTension(bestEdge);
            }
        }
        else {
            for (const int k : shortenPadEdges_) {
                const double t = EdgeTension(k);
                if (t > bestTension) {
                    bestTension = t;
                    bestEdge = k;
                }
            }
        }
        if (bestEdge < 0) {
            shortenSub_ = ShortenSub::NeedRemove;
            return;
        }
        const float length = EdgeChainLength(bestEdge);
        if (!r2Active_) {
            const double threshold =
                params_->tensionThreshold * std::sqrt(static_cast<double>(std::max(length, 1.0f)));
            if (bestTension < threshold) {
                shortenSub_ = ShortenSub::NeedRemove;
                return;
            }
        }
        r3SettleLength_ = length;
        if (InsertBeadAtFarthestGap(bestEdge) < 0) {
            shortenSub_ = ShortenSub::NeedRemove;
            return;
        }
        r3SettleEdge_ = bestEdge;
        r3SettleFrames_ = 0;
        shortenSub_ = ShortenSub::SettlingInsert;
        return;
    }

    if (shortenSub_ == ShortenSub::SettlingInsert) {
        if (R3SettleTick(r3SettleEdge_, r3SettleLength_)) {
            shortenSub_ = ShortenSub::NeedRemove;
        }
        return;
    }

    // Spoke: the X-Y spring (see Integrate) draws the vertex beads together.
    // The move fires on a short timer rather than on a settling threshold,
    // which could stall (or drag the diagram away) on a tense necklace.
    for (int i = 0; i < params_->autoStepsPerFrame; ++i) {
        Step();
    }
    ++spokeFrames_;
    const int x = r2Info_.vertexBead[0];
    const int y = r2Info_.vertexBead[1];
    if (x < 0 || y < 0 || x >= BeadCount() || y >= BeadCount()) {
        CancelShorten();
        return;
    }
    if (spokeFrames_ >= params_->r2SpokeFrames) {
        // Timer, not a settle threshold: pull X and Y to near-tangency so
        // the splice does not start with a long chord between them.
        const Vector2 px = PosOf(x);
        const Vector2 py = PosOf(y);
        const float dx = py.x - px.x;
        const float dy = py.y - px.y;
        const float len = std::sqrt(dx * dx + dy * dy);
        const float maxLen = 1.05f * diameter_;
        if (len > maxLen && len > 1e-6f) {
            const float extra = 0.5f * (len - maxLen);
            const Vector2 dir = {dx / len, dy / len};
            SetPosOf(x, {px.x + extra * dir.x, px.y + extra * dir.y});
            SetPosOf(y, {py.x - extra * dir.x, py.y - extra * dir.y});
        }
        r2Complete_ = true;
    }
}
