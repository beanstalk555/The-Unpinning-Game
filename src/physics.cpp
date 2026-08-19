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
        chains_.resize(static_cast<size_t>(loop_->numEdges + 1));
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

    // Hexagonal honeycomb spokes: face bead to the 6-cycle (rest = diameter).
    if (r3Active_ && R3HoneycombSpokesActive()) {
        const int faceBody = beads + r3HoneycombFaceFi_;
        if (faceBody >= 0 && faceBody < total) {
            for (int i = 0; i < 6; ++i) {
                const int spoke = r3HoneycombSpokes_[i];
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
                const float mag = params_->r3FortifyStiffness * (len - diameter_);
                const Vector2 dir = {d.x / len, d.y / len};
                force[faceBody].x += mag * dir.x;
                force[faceBody].y += mag * dir.y;
                force[spoke].x -= mag * dir.x;
                force[spoke].y -= mag * dir.y;
            }
        }
    }

    std::vector<char> skipFaceCollision(static_cast<size_t>(faces), 0);
    std::unordered_set<int> r3Faces;
    if (r3Active_) {
        CollectR3InvolvedFaces(&r3Faces);
    }
    // Face beads being pulled to their incenter skip collisions until they
    // arrive; crossing the polygon boundary is not enough (they still overlap
    // edge beads and would chatter).
    const float incenterArriveTol = 0.2f * diameter_;
    for (int fi = 0; fi < faces; ++fi) {
        const FaceBead& fb = neck_->faceBeads[static_cast<size_t>(fi)];
        std::vector<Vec2d> poly;
        const bool havePoly = FaceBoundaryPolygon(fb.face, &poly);
        const bool inside =
            !havePoly || PointInPolygon(Vec2d{static_cast<double>(fb.pos.x),
                                              static_cast<double>(fb.pos.y)},
                                        poly);
        bool pullToIncenter = !inside;
        if (r3Active_ && r3Faces.count(fb.face) != 0) {
            pullToIncenter = true;
        }
        if (R3HoneycombSpokesActive() && fb.face == r3Info_.face) {
            pullToIncenter = false;
        }
        if (!pullToIncenter) {
            continue;
        }
        if (havePoly) {
            const double precision = std::max(1e-4, 0.05 * static_cast<double>(radius_));
            const InscribedCircle circle = LargestInscribedCircle(poly, precision);
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
            const int body = beads + fi;
            const Vector2 p = PosOf(body);
            const Vector2 d = {static_cast<float>(target.x) - p.x,
                               static_cast<float>(target.y) - p.y};
            const float len = std::sqrt(d.x * d.x + d.y * d.y);
            if (len > 1e-6f) {
                // Rest length 0: force = stiffness * displacement toward the incenter.
                const float mag = params_->springStiffness * len;
                const Vector2 dir = {d.x / len, d.y / len};
                force[body].x += mag * dir.x;
                force[body].y += mag * dir.y;
            }
            if (len > incenterArriveTol) {
                skipFaceCollision[static_cast<size_t>(fi)] = 1;
            }
        }
        else if (!inside) {
            skipFaceCollision[static_cast<size_t>(fi)] = 1;
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
        if (r3Active_ && body >= beads) {
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

bool Physics::FaceBoundaryPolygon(int face, std::vector<Vec2d>* out) const
{
    if (out == nullptr || face < 0 || face >= static_cast<int>(loop_->faces.size()) ||
        face == exteriorFace_) {
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

void Physics::CollectR3InvolvedFaces(std::unordered_set<int>* out) const
{
    if (out == nullptr) {
        return;
    }
    out->clear();
    const auto addHalf = [&](int h) {
        const auto it = faceOfHalfEdge_.find(h);
        if (it != faceOfHalfEdge_.end() && it->second != exteriorFace_) {
            out->insert(it->second);
        }
    };
    if (r3Info_.face >= 0 && r3Info_.face != exteriorFace_) {
        out->insert(r3Info_.face);
    }
    for (int i = 0; i < 3; ++i) {
        addHalf(r3Info_.h[i]);
        addHalf(-r3Info_.h[i]);
        addHalf(r3Info_.a[i]);
        addHalf(-r3Info_.a[i]);
        addHalf(r3Info_.b[i]);
        addHalf(-r3Info_.b[i]);
    }
}

bool Physics::FaceBeadInsideFace(const FaceBead& fb) const
{
    std::vector<Vec2d> poly;
    if (!FaceBoundaryPolygon(fb.face, &poly)) {
        return true;
    }
    const Vec2d center{static_cast<double>(fb.pos.x), static_cast<double>(fb.pos.y)};
    return PointInPolygon(center, poly);
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
    // Loop edges (both ends the same crossing) must keep at least two edge
    // beads: reducing to one merges the +k/-k neighbor slots at that vertex
    // and subsequent reinserts corrupt the chain into a cycle.
    const bool loopEdge = EdgeIsLoop(edge);
    const int minBeads = loopEdge ? 3 : 2;
    if (edgeBeads < minBeads) {
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
    r3Meeting_[0] = r3Meeting_[1] = r3Meeting_[2] = -1;
    r3StallFrames_ = 0;
    r3SettleEdge_ = -1;
    r3SettleLength_ = 0.0f;
    r3SettleFrames_ = 0;
    r3Phase_ = R3Phase::PreHoneycomb;
    r3PreSub_ = R3PreSub::NeedRemove;
    r3HoneycombFaceFi_ = -1;
    for (int i = 0; i < 6; ++i) {
        r3HoneycombSpokes_[i] = -1;
    }
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
    const Vector2 fc = FaceBeadPosition(face);
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
        if (f == r3Info_.face || f == exteriorFace_) {
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
    std::unordered_set<int> faces;
    AppendIncidentFaces(edge, faces);
    for (const int f : faces) {
        if (f == exteriorFace_) {
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

bool Physics::InsertLargestGapIfNeeded()
{
    if (structureDirty_) {
        RebuildStructures();
    }
    const float gapThreshold = params_->gapSize * diameter_;
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

void Physics::StartR3PostRepair()
{
    if (autoActive_ || r3Active_) {
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
    if (++r3StallFrames_ > 800) {
        r3Active_ = false;
        r3Complete_ = false;
        ClearR3State();
        structureDirty_ = true;
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

    int n = 0;
    bool any = false;
    while (n++ < 8 && RemoveSharpBendIfAllowed()) {
        any = true;
    }
    Step();
    if (!any) {
        r3Active_ = false;
        r3Complete_ = false;
        ClearR3State();
        structureDirty_ = true;
    }
}

bool Physics::StartR3Move(int face)
{
    if (r3Active_ || autoActive_) {
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
        }
        return false;
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

        const int pick =
            candidates[static_cast<size_t>(std::rand() % static_cast<int>(candidates.size()))];
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
        int bestEdge = r3ExteriorEdges_[0];
        double bestTension = -1.0;
        for (const int k : r3ExteriorEdges_) {
            const double t = EdgeTension(k);
            if (t > bestTension) {
                bestTension = t;
                bestEdge = k;
            }
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
    if (r3HoneycombFaceFi_ < 0 ||
        r3HoneycombFaceFi_ >= static_cast<int>(neck_->faceBeads.size())) {
        return false;
    }
    for (int i = 0; i < 6; ++i) {
        if (r3HoneycombSpokes_[i] < 0) {
            return false;
        }
    }
    return true;
}

bool Physics::BindR3HoneycombSpokes()
{
    r3HoneycombFaceFi_ = -1;
    for (int i = 0; i < 6; ++i) {
        r3HoneycombSpokes_[i] = -1;
    }
    for (size_t fi = 0; fi < neck_->faceBeads.size(); ++fi) {
        if (neck_->faceBeads[fi].face == r3Info_.face) {
            r3HoneycombFaceFi_ = static_cast<int>(fi);
            break;
        }
    }
    if (r3HoneycombFaceFi_ < 0) {
        return false;
    }
    for (int i = 0; i < 3; ++i) {
        const int j = (i + 1) % 3;
        const std::vector<int> chain =
            WalkSideChain(r3Info_.vertexBead[i], r3Info_.h[i]);
        if (static_cast<int>(chain.size()) != 3 ||
            chain[0] != r3Info_.vertexBead[i] ||
            chain.back() != r3Info_.vertexBead[j] ||
            neck_->beads[static_cast<size_t>(chain[1])].isVertex ||
            neck_->beads[static_cast<size_t>(chain[1])].removed) {
            r3HoneycombFaceFi_ = -1;
            return false;
        }
        r3HoneycombSpokes_[2 * i] = r3Info_.vertexBead[i];
        r3HoneycombSpokes_[2 * i + 1] = chain[1];
    }
    return true;
}

bool Physics::R3HoneycombSpokesSettled() const
{
    if (!R3HoneycombSpokesActive()) {
        return false;
    }
    const int faceBody = BeadCount() + r3HoneycombFaceFi_;
    const float band = 0.35f * diameter_;
    const Vector2 pa = PosOf(faceBody);
    for (int i = 0; i < 6; ++i) {
        const int spoke = r3HoneycombSpokes_[i];
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

    // Unique exterior labels from the pre-move analysis (half-edge labels unchanged).
    std::vector<int> exteriors;
    for (int i = 0; i < 3; ++i) {
        for (const int h : {r3Info_.a[i], r3Info_.b[i]}) {
            const int k = std::abs(h);
            bool have = false;
            for (const int e : exteriors) {
                if (e == k) {
                    have = true;
                    break;
                }
            }
            if (!have) {
                exteriors.push_back(k);
            }
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
