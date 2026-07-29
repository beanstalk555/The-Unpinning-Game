#include "physics.h"

#include "geometry.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <unordered_set>

namespace
{
    constexpr float kPiF = 3.14159265358979323846f;

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
    : loop_(&loop),
      neck_(&necklace),
      params_(&params),
      exteriorFace_(exteriorFace),
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
    const int vStart = neck_->crossingBead[itStart->second];
    chain.push_back(vStart);

    int prev = vStart;
    int cur = neck_->NeighborToward(vStart, edge);
    int guard = 0;
    const int limit = BeadCount() + 8;
    while (cur >= 0 && !neck_->beads[cur].isVertex) {
        chain.push_back(cur);
        const Bead& b = neck_->beads[cur];
        const int nxt = (b.neighbors[0] == prev) ? b.neighbors[1] : b.neighbors[0];
        prev = cur;
        cur = nxt;
        if (++guard > limit) {
            break; // safety against a malformed link
        }
    }
    if (cur >= 0) {
        chain.push_back(cur); // tail vertex bead
    }
    return chain;
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
            es.push_back({chain[i], chain[i + 1]});
        }
    }
    FlattenSpringsFromEdges();
}

void Physics::RefreshEdgeSprings(int edge)
{
    if (edge < 1 || edge > loop_->numEdges) {
        return;
    }
    if (static_cast<int>(edgeSprings_.size()) < loop_->numEdges + 1) {
        edgeSprings_.assign(loop_->numEdges + 1, {});
    }
    const std::vector<int>& chain = chains_[edge];
    std::vector<std::pair<int, int>>& es = edgeSprings_[static_cast<size_t>(edge)];
    es.clear();
    for (size_t i = 0; i + 1 < chain.size(); ++i) {
        es.push_back({chain[i], chain[i + 1]});
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
}

void Physics::RebuildFaceGeometry()
{
    const ScopedTimer timer(debugTiming_ ? &debugTimings_.rebuildFaceGeometryMs : nullptr);
    sharpDotThreshold_ = std::cos(kPiF / 3.0f + params_->bendAngleEpsilon);

    neck_->sharpCornerEdges.clear();
    std::unordered_set<int> sharpSet;
    exteriorBoundary_.clear();
    for (int f = 0; f < static_cast<int>(loop_->faces.size()); ++f) {
        if (f == exteriorFace_) {
            std::vector<float> exteriorDots;
            WalkFaceBoundaryFromPhiCycle(loop_->faces[f], *neck_, crossingOf_,
                                         exteriorBoundary_, exteriorDots, &sharpSet,
                                         sharpDotThreshold_);
            continue;
        }
        FaceBead* fb = FaceBeadFor(f);
        if (fb == nullptr) {
            continue;
        }
        WalkFaceBoundaryFromPhiCycle(loop_->faces[f], *neck_, crossingOf_,
                                     fb->boundaryBeads, fb->boundaryDots, &sharpSet,
                                     sharpDotThreshold_);
        fb->cachedExcessArea = fb->ExcessArea(*neck_, radius_);
    }
    neck_->sharpCornerEdges.assign(sharpSet.begin(), sharpSet.end());
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
        faceBodies_.push_back(std::move(bodies));
    }
    if (!exteriorBoundary_.empty()) {
        faceBodies_.push_back(exteriorBoundary_);
    }
    PruneRemovedFromFaceBodies();
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

void Physics::AppendIncidentFaces(int edge, std::unordered_set<int>& out) const
{
    const auto itPos = faceOfHalfEdge_.find(edge);
    if (itPos != faceOfHalfEdge_.end() && itPos->second != exteriorFace_) {
        out.insert(itPos->second);
    }
    const auto itNeg = faceOfHalfEdge_.find(-edge);
    if (itNeg != faceOfHalfEdge_.end() && itNeg->second != exteriorFace_) {
        out.insert(itNeg->second);
    }
}

void Physics::WalkFaceBoundaryNoDots(int face)
{
    if (face == exteriorFace_) {
        std::vector<float> unusedDots;
        WalkFaceBoundaryFromPhiCycle(loop_->faces[face], *neck_, crossingOf_, exteriorBoundary_,
                                     unusedDots, nullptr, 0.0f, false);
        return;
    }
    FaceBead* fb = FaceBeadFor(face);
    if (fb == nullptr) {
        return;
    }
    WalkFaceBoundaryFromPhiCycle(loop_->faces[face], *neck_, crossingOf_, fb->boundaryBeads,
                                 fb->boundaryDots, nullptr, 0.0f, false);
}

void Physics::RefreshFaceBody(int face)
{
    if (face == exteriorFace_) {
        std::vector<int> bodies = exteriorBoundary_;
        bodies.erase(std::remove_if(bodies.begin(), bodies.end(),
                                    [&](int body) {
                                        return body < BeadCount() &&
                                               neck_->beads[static_cast<size_t>(body)].removed;
                                    }),
                     bodies.end());
        const size_t interiorCount = neck_->faceBeads.size();
        if (!bodies.empty()) {
            if (faceBodies_.size() > interiorCount) {
                faceBodies_[interiorCount] = std::move(bodies);
            }
            else {
                faceBodies_.push_back(std::move(bodies));
            }
        }
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
        return;
    }
}

void Physics::RefreshFaceBodiesForEdge(int edge)
{
    std::unordered_set<int> faces;
    AppendIncidentFaces(edge, faces);
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
        chains_.assign(loop_->numEdges + 1, {});
    }
    chains_[edge] = WalkEdge(edge);
}

void Physics::RebuildStructures()
{
    const ScopedTimer timer(debugTiming_ ? &debugTimings_.rebuildStructuresMs : nullptr);
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

void Physics::ApplyFaceCollisions(std::vector<Vector2>& force, int total)
{
    const ScopedTimer timer(debugTiming_ ? &debugTimings_.collisionMs : nullptr);
    collisionGrid_.SetCellSize(params_->collisionBoundsScale * diameter_);
    std::unordered_set<long long> done;
    const long long stride = total;

    for (const std::vector<int>& group : faceBodies_) {
        collisionGrid_.Clear();
        collisionGrid_.Reserve(group.size());
        for (const int body : group) {
            const Vector2 p = PosOf(body);
            collisionGrid_.Insert(body, p.x, p.y);
        }
        collisionGrid_.ForEachCandidatePair([&](int a, int b) {
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
        });
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

    // Hard-body collisions within each face (spatial broad phase).
    ApplyFaceCollisions(force, total);

    // Drag toward the cursor.
    if (dragBody_ >= 0 && dragBody_ < total) {
        const Vector2 p = PosOf(dragBody_);
        force[dragBody_].x += params_->dragStiffness * (dragTarget_.x - p.x);
        force[dragBody_].y += params_->dragStiffness * (dragTarget_.y - p.y);
    }

    // Damping + semi-implicit Euler integration (mass = 1, removed beads frozen).
    for (int body = 0; body < total; ++body) {
        if (body < beads && neck_->beads[body].removed) {
            continue;
        }
        Vector2& v = VelOf(body);
        Vector2 f = force[body];
        f.x -= params_->damping * v.x;
        f.y -= params_->damping * v.y;
        v.x += f.x * dt;
        v.y += f.y * dt;
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
        if (EdgeTension(pending.edge) > threshold * 0.5) {
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
        if (r3Active_) {
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

void Physics::ReinsertEdgeBead(int bead, int neighborA, int neighborC)
{
    Bead& b = neck_->beads[bead];
    const Vector2 pa = PosOf(neighborA);
    const Vector2 pc = PosOf(neighborC);
    b.removed = false;
    b.pos = {(pa.x + pc.x) * 0.5f, (pa.y + pc.y) * 0.5f};
    b.neighbors[0] = neighborA;
    b.neighbors[1] = neighborC;
    ReplaceNeighbor(neighborA, neighborC, bead);
    ReplaceNeighbor(neighborC, neighborA, bead);
    if (bead >= 0 && bead < static_cast<int>(vel_.size())) {
        vel_[bead] = Vector2{0.0f, 0.0f};
    }
    structureDirty_ = true;
}

void Physics::RemoveEdgeBead(int bead)
{
    Bead& b = neck_->beads[bead];
    const int a = b.neighbors[0];
    const int c = b.neighbors[1];
    ReplaceNeighbor(a, bead, c);
    ReplaceNeighbor(c, bead, a);
    b.removed = true;
    structureDirty_ = true;
}

int Physics::InsertEdgeBeadBetween(int neighborA, int neighborC, int edge)
{
    for (int i = 0; i < BeadCount(); ++i) {
        Bead& slot = neck_->beads[static_cast<size_t>(i)];
        if (!slot.removed || slot.isVertex || slot.edge != edge) {
            continue;
        }
        ReinsertEdgeBead(i, neighborA, neighborC);
        return i;
    }

    Bead bead;
    bead.isVertex = false;
    bead.removed = false;
    bead.edge = edge;
    const Vector2 pa = PosOf(neighborA);
    const Vector2 pc = PosOf(neighborC);
    bead.pos = {(pa.x + pc.x) * 0.5f, (pa.y + pc.y) * 0.5f};
    bead.neighbors[0] = neighborA;
    bead.neighbors[1] = neighborC;
    const int idx = BeadCount();
    ReplaceNeighbor(neighborA, neighborC, idx);
    ReplaceNeighbor(neighborC, neighborA, idx);
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
        if (r3Active_) {
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
    for (FaceBead& fb : neck_->faceBeads) {
        if (fb.face == face) {
            return &fb;
        }
    }
    return nullptr;
}

const FaceBead* Physics::FaceBeadFor(int face) const
{
    for (const FaceBead& fb : neck_->faceBeads) {
        if (fb.face == face) {
            return &fb;
        }
    }
    return nullptr;
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
    if (face == exteriorFace_) {
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
    if (face == exteriorFace_) {
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
        if (f == exteriorFace_) {
            continue;
        }
        FaceBead* fb = FaceBeadFor(f);
        if (fb == nullptr) {
            continue;
        }
        if (!boundariesCurrent && refreshFaces.count(f) != 0) {
            WalkFaceBoundaryNoDots(f);
        }
        if (fb->boundaryBeads.size() > 6) {
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
    if (face == exteriorFace_) {
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
    if (itPos != faceOfHalfEdge_.end() && itPos->second != exteriorFace_) {
        const FaceBead* fb = FaceBeadFor(itPos->second);
        if (fb != nullptr) {
            excesses.push_back(fb->cachedExcessArea);
        }
    }
    if (itNeg != faceOfHalfEdge_.end() && itNeg->second != exteriorFace_) {
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
    if (edgeBeads < 2) {
        return -1;
    }

    const float maxDot = sharpDotThreshold_;

    if (bendPassOnly) {
        for (size_t i = 1; i + 1 < chain.size(); ++i) {
            if (Dot(chain[i]) >= maxDot) {
                const int bead = chain[i];
                const int a = neck_->beads[bead].neighbors[0];
                const int c = neck_->beads[bead].neighbors[1];
                RemoveEdgeBead(bead);
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
    RemoveEdgeBead(bestBead);
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
        if (edgeBeadsRemaining[static_cast<size_t>(k)] < 2) {
            continue;
        }
        const float length = edgeLength(chains_[k]);
        const int a = neck_->beads[bead].neighbors[0];
        const int c = neck_->beads[bead].neighbors[1];
        RemoveEdgeBead(bead);
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
        if (edgeBeadsRemaining[static_cast<size_t>(k)] < 2) {
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
            RefreshEdgeSprings(k);
            AppendIncidentFaces(k, affectedFaces);
        }
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

    if (pendingEdges_.empty()) {
        return false;
    }

    framesWaited_ = 0;
    return true;
}

void Physics::StartAutoTighten()
{
    if (autoActive_ || r3Active_) {
        return;
    }
    autoActive_ = true;
    needsGeometryEval_ = true;
    removedEdgesLastRound_.clear();
    struckEdges_.assign(loop_->numEdges + 1, false);
    pendingEdges_.clear();
    framesWaited_ = 0;
    if (!BeginRemovalRound()) {
        autoActive_ = false;
    }
}

void Physics::StopAutoTighten()
{
    autoActive_ = false;
    needsGeometryEval_ = false;
    pendingEdges_.clear();
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
        if (EdgeTension(pending.edge) < threshold) {
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
                    if (pending.edge == k && pending.removedBead >= 0) {
                        ReinsertEdgeBead(pending.removedBead, pending.neighborA, pending.neighborC);
                        pending.settled = true;
                        reinsertedEdges.insert(k);
                    }
                }
                if (k >= 1 && k < static_cast<int>(struckEdges_.size())) {
                    struckEdges_[k] = true;
                }
            }
            for (const int k : reinsertedEdges) {
                RefreshEdgeChain(k);
                RefreshEdgeSprings(k);
                RefreshFaceBodiesForEdge(k);
                removedEdgesLastRound_.insert(k);
            }
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
// Reidemeister 3
// ---------------------------------------------------------------------------

bool Physics::StartR3Move(int face)
{
    if (r3Active_ || autoActive_) {
        return false;
    }
    R3FaceInfo info;
    if (!AnalyzeR3Face(*loop_, neck_, face, &info)) {
        return false;
    }

    // Orientation of the face walk in grid coordinates decides between the
    // two mirror-image sigma rewrites.
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

    // The rewritten sigma must be valid before we cut anything.
    {
        Permutation newSigma;
        if (!ComputeR3Sigma(*loop_, info, &newSigma)) {
            return false;
        }
    }

    // Ordered beads along one side, walked from a triangle vertex in the
    // direction of the given half-edge.
    const auto walkSide = [&](int vertexBead, int halfEdge) {
        std::vector<int> chain;
        chain.push_back(vertexBead);
        int prev = vertexBead;
        int cur = neck_->NeighborToward(vertexBead, halfEdge);
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
    };

    // Each side needs room for a meeting bead: at least three interior beads.
    for (int i = 0; i < 3; ++i) {
        for (int guard = 0; guard < 8; ++guard) {
            const std::vector<int> chain = walkSide(info.vertexBead[i], info.h[i]);
            if (chain.size() < 2) {
                return false;
            }
            if (static_cast<int>(chain.size()) - 2 >= 3) {
                break;
            }
            InsertEdgeBeadBetween(chain[0], chain[1], std::abs(info.h[i]));
        }
    }

    // Fresh structures pre-cut; the collision groups built here stay frozen
    // for the whole slide (each side's outer face already contains both
    // sliding branches and the side beads, so the groups cover the move).
    RebuildStructures();

    for (int i = 0; i < 3; ++i) {
        const int j = (i + 1) % 3;
        const std::vector<int> chain = walkSide(info.vertexBead[i], info.h[i]);
        if (static_cast<int>(chain.size()) < 5 || chain.back() != info.vertexBead[j]) {
            return false;
        }
        // bTerminal[i] slides forward along side i, from X_i toward X_{i+1}.
        R3Slider& f = r3Sliders_[2 * i];
        f = R3Slider{};
        f.terminal = info.bTerminal[i];
        f.cutVertex = info.vertexBead[i];
        f.path.assign(chain.begin(), chain.end() - 1);
        // aTerminal[i+1] slides backward along side i, from X_{i+1} toward X_i.
        R3Slider& g = r3Sliders_[2 * i + 1];
        g = R3Slider{};
        g.terminal = info.aTerminal[j];
        g.cutVertex = info.vertexBead[j];
        g.path.assign(chain.rbegin(), chain.rend() - 1);
    }

    // Cut the six terminal links (identified by half-edge, so kink loops at a
    // corner are cut on the correct slot).
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
        cutsOk = cutLink(r3Sliders_[2 * i], info.b[i]) &&
                 cutLink(r3Sliders_[2 * i + 1], info.a[j]);
    }
    for (int i = 0; i < 3; ++i) {
        const int j = (i + 1) % 3;
        const int bEdge = std::abs(info.b[i]);
        const int aEdge = std::abs(info.a[j]);
        bool haveB = false;
        bool haveA = false;
        for (const auto& entry : r3CutEdges_) {
            haveB = haveB || entry.first == bEdge;
            haveA = haveA || entry.first == aEdge;
        }
        if (!haveB) {
            r3CutEdges_.push_back({bEdge, info.bTerminal[i]});
        }
        if (!haveA) {
            r3CutEdges_.push_back({aEdge, info.aTerminal[j]});
        }
    }

    r3Info_ = info;
    r3Active_ = true;
    r3Complete_ = false;
    r3StallFrames_ = 0;
    r3Meeting_[0] = r3Meeting_[1] = r3Meeting_[2] = -1;
    structureDirty_ = true;

    if (!cutsOk) {
        CancelR3Move();
        return false;
    }
    return true;
}

void Physics::CancelR3Move()
{
    if (!r3Active_) {
        return;
    }
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
    r3Active_ = false;
    r3Complete_ = false;
    r3CutEdges_.clear();
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
    const int n = loop_->numEdges;
    chains_.assign(n + 1, {});
    for (int k = 1; k <= n; ++k) {
        chains_[k] = WalkEdge(k);
    }
    // Cut edges no longer reach both vertices; walk them from a terminal so
    // their whole bead chain keeps its springs.
    for (const auto& entry : r3CutEdges_) {
        chains_[entry.first] = WalkFromTerminal(entry.second);
    }
    RebuildSpringsFromChains();
    // faceBodies_ intentionally left frozen from the pre-cut rebuild.
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
            const int inserted = InsertEdgeBeadBetween(chain[i], chain[i + 1], entry.first);
            AddBeadToFrozenGroups(inserted, chain[i]);
            break; // one insert per branch per frame
        }
    }
}

bool Physics::R3SliderSettled(const R3Slider& s) const
{
    const Vector2 pa = PosOf(s.terminal);
    const Vector2 pb = PosOf(s.path[static_cast<size_t>(s.targetIdx)]);
    const float dx = pb.x - pa.x;
    const float dy = pb.y - pa.y;
    return std::sqrt(dx * dx + dy * dy) < 1.35f * diameter_;
}

void Physics::AdvanceR3Slider(R3Slider& s)
{
    ++s.targetIdx;
    const Vector2 from = PosOf(s.path[static_cast<size_t>(s.targetIdx) - 1]);
    const Vector2 to = PosOf(s.path[static_cast<size_t>(s.targetIdx)]);
    const float dx = to.x - from.x;
    const float dy = to.y - from.y;
    const float len = std::sqrt(dx * dx + dy * dy);
    if (len < 1e-6f) {
        return;
    }
    // A boost tangent to the boundary helps the terminal hop past collisions.
    const float boost = 0.5f * params_->maxSpeed;
    vel_[s.terminal].x += boost * dx / len;
    vel_[s.terminal].y += boost * dy / len;
}

void Physics::R3Step()
{
    if (!r3Active_ || r3Complete_) {
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
        if (f.targetIdx + g.targetIdx >= m + 1) {
            continue; // this pair already pulls toward a shared bead
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
            for (int i = 0; i < 3; ++i) {
                const R3Slider& f = r3Sliders_[2 * i];
                r3Meeting_[i] = f.path[static_cast<size_t>(f.targetIdx)];
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
