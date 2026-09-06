#pragma once

#include <memory>
#include <vector>

#include "beads.h"
#include "move_rewrite.h"
#include "multiloop.h"

// Combinatorial data for a shortening Reidemeister 1 move on a monogon.
//
// The face is the single half-edge e with sigma(e) = -e at crossing X, whose
// cycle is therefore (e, -e, p, q). Straight-through pairs are {e, p} and
// {-e, q}: the strand arrives on p, loops around e, and leaves on q, so
// removing the kink splices p to q.
struct R1FaceInfo
{
    int face = -1;
    int e = 0;
    int p = 0;
    int q = 0;
    int crossing = -1;
    int vertexBead = -1;
    int pTerminal = -1; // edge bead adjacent to X toward p
    int qTerminal = -1; // edge bead adjacent to X toward q
};

// Returns false unless the face is a monogon with the expected cycle shape.
// With a necklace, also fills vertexBead / terminals (false if a terminal is
// missing, shared, or not an edge bead).
bool AnalyzeR1Face(const Multiloop& loop, const Necklace* neck, int face, R1FaceInfo* out);

// Removes the kink: deletes X and edge |e|, splices p to q, and turns the
// vertex bead into an edge bead of the merged edge (the loop's remaining
// beads are removed). Face beads and the component list are updated through
// FinishMoveBookkeeping. Returns nullptr and leaves the necklace untouched
// on any inconsistency.
std::unique_ptr<Multiloop> FinalizeR1(const Multiloop& oldLoop, Necklace& neck,
                                      const R1FaceInfo& info,
                                      const std::vector<MultiloopComponent>& oldComponents,
                                      MoveComponentResult* out);
