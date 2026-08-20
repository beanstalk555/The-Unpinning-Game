#pragma once

#include <raylib.h>

#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "beads.h"
#include "multiloop.h"
#include "physics_params.h"
#include "r3_move.h"
#include "spatial_disk_grid.h"
#include "geometry.h"

// Millisecond accumulators for the debug timing overlay (reset each displayed frame).
struct PhysicsDebugTimings
{
    double frameTotalMs = 0.0;
    double stepMs = 0.0;
    int stepCount = 0;
    double integrateMs = 0.0;
    double collisionMs = 0.0;
    double rebuildStructuresMs = 0.0;
    double rebuildChainsMs = 0.0;
    double rebuildFaceGeometryMs = 0.0;
    double rebuildFaceBodiesMs = 0.0;
    double strikeHoneycombMs = 0.0;
    double beginRemovalRoundMs = 0.0;
    double autoStepOverheadMs = 0.0;
    int substepsTotal = 0;
};

// Spring-and-collision relaxation for a bead necklace.
//
// The simulation works directly on the necklace's bead positions, in grid
// coordinates (the same units as GridLayout, with the bead radius measured in
// those units). Springs of rest length equal to the bead diameter connect
// consecutive beads, so the rest configuration is exactly the tangent grid
// drawing (tension 0). Face beads float in their regions via hard-body
// collisions with the surrounding beads.
//
// Parameters are read from PhysicsParams each step (adjustable at runtime).
class Physics
{
public:
    Physics(const Multiloop& loop, Necklace& necklace, PhysicsParams& params,
            int exteriorFace);

    // Advance the simulation by one displayed frame.
    void Step();

    // Diagram tension, cached from the most recent Step()/AutoStep().
    // Stretch and compression are averaged over non-removed necklace beads.
    double StretchTension() const { return stretch_; }
    double CompressionTension() const { return compression_; }
    double TotalTension() const { return stretch_ + compression_; }

    // Cached tangle flag: unintended non-crossing at a vertex, or two
    // necklace edges whose center segments properly intersect. Refresh after
    // drag / auto-tighten / R3; false on a freshly built planar layout.
    bool Tangled() const { return tangled_; }
    void RefreshTangle();

    // Bend at an edge bead B sitting between neighbors A and C: the dot
    // product of the unit vectors BA and BC. +1 is folded back on itself,
    // 0 is a right angle, and -1 is perfectly straight.
    float Dot(int edgeBead) const;

    // ---- click and drag (grid coordinates) ----
    // Body ids: 0..beadCount-1 are necklace beads, then face beads.
    int PickBody(Vector2 gridPos, float maxGridDist) const;
    void BeginDrag(int body, Vector2 gridPos);
    void Drag(Vector2 gridPos);
    void EndDrag();
    bool Dragging() const { return dragBody_ >= 0; }

    // ---- auto-tighten ----
    void StartAutoTighten();
    void StopAutoTighten();
    bool AutoTightening() const { return autoActive_; }
    void AutoStep();

    // ---- Reidemeister 3 ----
    // Starts the animated move on a triangular face. Returns false when the
    // face does not qualify or another mode is running.
    bool StartR3Move(int face);
    // Reverts the cuts (if any) and abandons the move.
    void CancelR3Move();
    bool R3Active() const { return r3Active_; }
    bool R3Complete() const { return r3Complete_; }
    // Advance the pre-honeycomb / sliding animation by one displayed frame.
    void R3Step();
    // Valid once R3Complete(): the analyzed face and the meeting bead on each
    // side (the new crossing between X_i and X_{i+1} sits at meeting[i]).
    void GetR3Result(R3FaceInfo* info, int meeting[3]) const;
    // After FinalizeR3: top up exterior arcs sa–sf so monogons have ≥5 edge
    // beads and bigons have ≥2. Uses newLoop for topology; edits neck_.
    void EnsureR3ExteriorHoneycombMins(const Multiloop& newLoop);
    // After FinalizeR3 + a new Physics on the rewritten loop: fill oversized
    // gaps, then strip sharp bends without breaking honeycomb faces.
    void StartR3PostRepair();

    // Number of 1-/2-/3-gons currently held by a 6-spoke honeycomb truss.
    int HoneycombFortifyCount() const { return static_cast<int>(honeycombFaceFi_.size()); }

    // Excess area (polyArea minus boundary pie slices) for a face bead.
    double FaceExcessArea(const FaceBead& fb) const;

    // True when the face bead's center is inside its face polygon (even-odd raycast).
    bool FaceBeadInsideFace(const FaceBead& fb) const;

    // Same, looked up by combinatorial face index. Exterior face is infinite.
    double ExcessArea(int face) const;

    // PolyArea / (boundary bead count)^2. Exterior face is infinite.
    double IsoperimetricRatio(int face) const;

    // Debug: print shoelace polygon vertices (bead index, type, position) to stdout.
    void PrintShoelaceDebug(int face) const;

    void ResetDebugTimings();
    void SetDebugFrameTotalMs(double ms);
    const PhysicsDebugTimings& DebugTimings() const { return debugTimings_; }

    // Debug: collide every necklace bead with every other necklace bead
    // (all-pairs), ignoring face groups. Face beads still use per-face groups.
    void SetGlobalCollisions(bool enabled) { globalCollisions_ = enabled; }
    bool GlobalCollisions() const { return globalCollisions_; }

private:
    int BeadCount() const { return static_cast<int>(neck_->beads.size()); }
    Vector2 PosOf(int body) const;
    void SetPosOf(int body, Vector2 p);
    Vector2& VelOf(int body);

    void RebuildStructures();
    void RebuildSpringsFromChains();
    void RefreshEdgeSprings(int edge);
    void FlattenSpringsFromEdges();
    void RebuildFaceGeometry();
    void RebuildFaceBodiesFromCache();
    void RefreshFaceBody(int face);
    void RefreshFaceBodiesForEdge(int edge);
    void PruneRemovedFromFaceBodies();
    void RefreshEdgeChain(int edge);
    void WalkFaceBoundaryNoDots(int face);
    void AppendIncidentFaces(int edge, std::unordered_set<int>& out) const;
    void Integrate(float dt);
    void ComputeTension();
    int SubstepsForStep() const;
    void ApplyFaceCollisions(std::vector<Vector2>& force, int total,
                             const std::vector<char>& skipFaceBody);
    // Average spring tension along one edge's bead chain (stretch + compression).
    double EdgeTension(int edge) const;

    // Ordered beads of an edge: [tail vertex, edge beads..., head vertex].
    std::vector<int> WalkEdge(int edge) const;
    void RemoveEdgeBead(int bead);
    void ReinsertEdgeBead(int bead, int neighborA, int neighborC);
    // Returns the (possibly recycled) bead index of the inserted bead.
    int InsertEdgeBeadBetween(int neighborA, int neighborC, int edge);
    void ReplaceNeighbor(int node, int oldNeighbor, int newNeighbor);
    void ReplaceNeighborOnce(int node, int oldNeighbor, int newNeighbor);
    bool NodeReferences(int node, int target) const;
    int CountReferences(int node, int target) const;
    // True when +edge and -edge meet at the same crossing (1-gon kink).
    bool EdgeIsLoop(int edge) const;

    // Split consecutive pairs when center-to-center distance exceeds gapSize bead diameters.
    void MaintainBeadSpacingDuringDrag();

    bool BeginRemovalRound(); // remove a bead on each qualifying edge; false if none

    // Remove one bead from edge k (bend-first, else excess-area). Returns the
    // removed bead index, or -1. On success, neighborA/C are its former neighbors.
    int RemoveBeadFromEdge(int edge, bool bendPassOnly, int* neighborA = nullptr,
                           int* neighborC = nullptr, bool skipExcessCheck = false);
    double EdgeExcessScore(int edge) const;

    FaceBead* FaceBeadFor(int face);
    const FaceBead* FaceBeadFor(int face) const;
    Vector2 FaceBeadPosition(int face) const;

    // Mark edges bounding faces with at most six boundary beads as struck.
    // When boundariesCurrent is true, boundaryBeads are already fresh (full geometry rebuild).
    // Otherwise only refreshFaces get a dot-free boundary walk before the count check.
    void StrikeHoneycombEdges(bool boundariesCurrent,
                              const std::unordered_set<int>& refreshFaces);

    const Multiloop* loop_;
    Necklace* neck_;
    PhysicsParams* params_;
    int exteriorFace_ = -1;
    std::unordered_map<int, int> faceOfHalfEdge_;
    float radius_;
    float diameter_;

    std::vector<Vector2> vel_;                  // per necklace bead
    std::vector<Vector2> fvel_;                 // per face bead
    std::unordered_map<int, int> crossingOf_;   // half-edge -> crossing index

    // Rebuilt when topology changes (bead removal/reinsert).
    std::vector<std::vector<int>> chains_;      // chains_[k], k = 1..numEdges
    std::vector<std::pair<int, int>> springs_;  // body pairs, rest = diameter
    std::vector<std::vector<std::pair<int, int>>> edgeSprings_; // edgeSprings_[k]
    std::vector<std::vector<int>> faceBodies_;  // per face: [faceBody?, boundary...]
    std::vector<int> exteriorBoundary_;
    SpatialDiskGrid collisionGrid_;
    bool structureDirty_ = true;
    float sharpDotThreshold_ = -1.0f;

    double stretch_ = 0.0;
    double compression_ = 0.0;
    bool tangled_ = false;

    int dragBody_ = -1;
    Vector2 dragTarget_ = {0.0f, 0.0f};

    // Auto-tighten state machine.
    bool autoActive_ = false;
    struct PendingEdge
    {
        int edge = -1;
        float length = 0.0f;
        bool settled = false;
        int removedBead = -1;
        int neighborA = -1;
        int neighborC = -1;
    };
    std::vector<PendingEdge> pendingEdges_;
    std::vector<bool> struckEdges_;
    int framesWaited_ = 0;
    // True at auto-tighten start and after a failed reinsert; triggers full face geometry.
    bool needsGeometryEval_ = false;
    // Edges that lost or regained a bead in the previous removal round (for incremental updates).
    std::unordered_set<int> removedEdgesLastRound_;

    // ---- Reidemeister 3 state machine ----
    // Phase A: tighten the triangle to a honeycomb (one edge bead per side).
    // Phase A2: hex-spoke springs from the face bead to the 6-cycle.
    // Phase B: cut-and-slide terminals toward those three beads (meetings).
    // Phase C: after FinalizeR3, fill gaps then strip sharp bends.
    enum class R3Phase
    {
        PreHoneycomb,
        Fortify,
        Slide,
        PostRepair
    };
    enum class R3PreSub
    {
        NeedRemove,
        SettlingRemove,
        NeedInsert,
        SettlingInsert
    };

    // One sliding terminal: an edge bead cut from its old triangle vertex,
    // spring-tied to successive boundary beads along its slide path.
    struct R3Slider
    {
        int terminal = -1;
        int cutVertex = -1;       // old vertex bead it was cut from
        int cutVertexSlot = -1;   // slot in the vertex's neighborOfHalfEdge
        int cutNeighborSlot = -1; // slot in the terminal's neighbors
        std::vector<int> path;    // [old vertex, side beads...] in slide order
        int targetIdx = 0;        // current spring target = path[targetIdx]
    };

    void RebuildStructuresR3();
    std::vector<int> WalkFromTerminal(int terminal) const;
    std::vector<int> WalkSideChain(int vertexBead, int halfEdge) const;
    void AddBeadToFrozenGroups(int bead, int adjacent);
    void MaintainR3BranchSpacing();
    bool R3SliderSettled(const R3Slider& s) const;
    void AdvanceR3Slider(R3Slider& s);
    void R3PreHoneycombStep();
    bool BeginR3Slide();
    bool BindR3HoneycombSpokes();
    bool R3HoneycombSpokesActive() const;
    bool R3HoneycombSpokesSettled() const;
    void R3FortifyStep();
    void BeginFortifyOrSlide();
    int EdgeInteriorCount(int edge) const;
    float EdgeChainLength(int edge) const;
    int InsertBeadAtFarthestGap(int edge);
    int RemoveFarthestFromFaceCenter(int edge, int face);
    bool R3SettleTick(int edge, float lengthAtStart);
    // True when removing from edge would strip a neighbor face that is already
    // at the honeycomb floor (≤6 boundary beads). The clicked R3 face is ignored.
    bool NeighborHoneycombBlocksRemoval(int edge);
    bool HoneycombFloorBlocksRemoval(int edge);
    bool InsertLargestGapIfNeeded();
    bool RemoveSharpBendIfAllowed();
    bool RemoveMonogonExcessIfAllowed();
    void RebuildHoneycombSpokes();
    bool FaceHasHoneycombSpokes(int fi) const;
    int HoneycombSlotForFace(int face) const;
    void R3PostRepairStep();
    void ClearR3State();
    bool FaceBoundaryPolygon(int face, std::vector<Vec2d>* out) const;
    void CollectR3InvolvedFaces(std::unordered_set<int>* out) const;
    bool DetectTangle() const;

    bool r3Active_ = false;
    bool r3Complete_ = false;
    R3Phase r3Phase_ = R3Phase::PreHoneycomb;
    R3PreSub r3PreSub_ = R3PreSub::NeedRemove;
    R3FaceInfo r3Info_;
    // Sliders 2i / 2i+1 slide along side i: bTerminal[i] forward from X_i and
    // aTerminal[i+1] backward from X_{i+1}; they meet at the new crossing.
    R3Slider r3Sliders_[6];
    int r3Meeting_[3] = {-1, -1, -1};
    int r3StallFrames_ = 0;
    std::vector<std::pair<int, int>> r3CutEdges_; // (edge label, terminal bead)
    // Unique exterior edge labels |a[i]|, |b[i]| (sa–sf, possibly fewer than 6).
    std::vector<int> r3ExteriorEdges_;
    int r3SettleEdge_ = -1;
    float r3SettleLength_ = 0.0f;
    int r3SettleFrames_ = 0;
    // Packed 6-spoke truss for every 1-/2-/3-gon with 6 boundary beads.
    std::vector<int> honeycombFaceFi_;
    std::vector<int> honeycombSpokes_;

    bool debugTiming_ = false;
    PhysicsDebugTimings debugTimings_;
    bool globalCollisions_ = false;
};
