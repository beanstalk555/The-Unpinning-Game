#pragma once

#include <unordered_map>
#include <vector>

#include "beads.h"
#include "multiloop.h"
#include "permutation.h"

// Shared machinery for shortening Reidemeister moves (R1, R2): delete
// crossings and edges, splice the surviving arcs through the deleted
// crossings, and compact half-edge labels so the result is again a valid
// multiloop sigma (labels exactly +-1..+-m).

// Two half-edges (each the end of a surviving edge at a deleted crossing)
// whose arcs join once the crossing is gone.
struct SplicePair
{
    int u = 0;
    int v = 0;
};

struct RewriteResult
{
    Permutation sigma; // compact labels 1..numNewEdges
    int numNewEdges = 0;
    // Signed old half-edge (both signs of every surviving edge) -> signed new
    // label. The sign encodes whether the old edge runs with (+) or against
    // (-) the orientation of the merged edge it belongs to.
    std::unordered_map<int, int> oldToNew;
    std::vector<int> deletedEdges;
    // mergedArcs[m] (m = 1..numNewEdges; entry 0 unused): the old signed
    // half-edges traversed in +m order, one per constituent old edge.
    std::vector<std::vector<int>> mergedArcs;
    // Arcs that closed up without meeting a surviving crossing.
    std::vector<std::vector<int>> circleArcs;
    std::vector<int> circleStrand; // old strand id per circle
    // Old signed half-edge (both signs) on a circle -> circle index.
    std::unordered_map<int, int> circleOfHalfEdge;

    // New signed label for an old signed half-edge, or 0 if deleted / circle.
    int MapHalfEdge(int h) const;
    // Old edge label -> new edge label (0 if gone or on a circle).
    int MapEdge(int k) const;
};

bool DeleteCrossingsAndSplice(const Multiloop& old, const std::vector<int>& deletedCrossings,
                              const std::vector<int>& deletedEdges,
                              const std::vector<SplicePair>& splices, RewriteResult* out);

// Copies strand ids from the old loop onto the new one (max-overlap through
// the merged arcs), keeps numStrands at the old count so strands that became
// circles keep their color slot free.
void RemapStrandsThroughRewrite(const Multiloop& oldLoop, const RewriteResult& rw,
                                Multiloop* newLoop);

// Half-edges of old face `oldFace` that survive (mapped through oldToNew),
// excluding those on circles.
std::vector<int> MappedFaceHalfEdges(const Multiloop& oldLoop, const RewriteResult& rw,
                                     int oldFace);

// Greedy max-overlap matching of old faces onto new faces. `forbidden` new
// faces are never used. Returns -1 for old faces without any overlap.
std::vector<int> MatchFacesByOverlap(const Multiloop& oldLoop, const Multiloop& newLoop,
                                     const RewriteResult& rw, const std::vector<int>& oldFaces,
                                     const std::vector<int>& forbidden);

// Rewires the necklace once the move-specific surgery has relinked the beads
// at the deleted crossings: remaps vertex sigma cycles, rebuilds crossing
// indices, walks every new edge from its +m crossing to normalize edge labels
// and neighbor order, and turns closed arcs into circle rings (edge 0, strand
// from circleStrand). Face beads are untouched. Two-phase: returns false
// (leaving the necklace unchanged) when any walk is inconsistent.
bool ApplyRewriteToNecklace(Necklace& neck, const RewriteResult& rw, const Multiloop& newLoop);

// Refreshes phiCycle / vertexBeads / cached boundary of a face bead that was
// moved onto new face `face`.
void RebindFaceBead(FaceBead& fb, int face, const Multiloop& newLoop, const Necklace& neck);

// Ring of live circle beads through `bead` (edge 0), or empty if it is not a
// closed ring of edge beads.
std::vector<int> CircleRingThrough(const Necklace& neck, int bead);

// Resulting component list after a shortening move on `oldComponents`.
struct MoveComponentResult
{
    std::vector<MultiloopComponent> components;
    // Face-bead labels of a pending merge (two beads on the same new face);
    // -1 when no merge is needed.
    int mergeKeepLabel = -1;
    int mergeDropLabel = -1;
    // Old face index -> new face index (-1 if gone or split), for every old face.
    std::vector<int> faceMap;
    // New edge labels of the exterior arcs (for honeycomb minimum top-up).
    std::vector<int> exteriorArcs;
    // True when the moved component split into two pieces.
    bool disconnected = false;
};

struct MoveBookkeepingInput
{
    int movedComponent = -1;
    // Disconnecting R2: the old face whose half-edges end up in both pieces.
    int splitFace = -1;
    // Connected R2: the two old faces that merge (with the bigon) into one.
    int mergeFaceA = -1;
    int mergeFaceB = -1;
    // Old half-edges (surviving) whose images identify the two pieces of a
    // split face, or lie on the merged face. For R1 leave 0.
    int pieceHalfA = 0;
    int pieceHalfB = 0;
    // Old edge labels of the exterior arcs adjacent to the deleted crossings.
    std::vector<int> exteriorArcEdges;
};

// After ApplyRewriteToNecklace: moves face beads onto the new faces (or turns
// them into circle pads), erases beads whose region vanished or became
// exterior, and rebuilds the component list (pieces replace the moved
// component). Nesting uses bead positions: the piece whose ring or split-face
// polygon contains the other is the hole; holes of the split face that sit
// inside the new inner piece nest there (concentric circles). Combinatorial
// exterior-ownership is the fallback when geometry is ambiguous.
bool FinishMoveBookkeeping(const Multiloop& oldLoop, const Multiloop& newLoop,
                           const RewriteResult& rw, Necklace& neck,
                           const std::vector<MultiloopComponent>& oldComponents,
                           const MoveBookkeepingInput& in, MoveComponentResult* out);
