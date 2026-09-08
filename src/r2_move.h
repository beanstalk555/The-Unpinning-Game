#pragma once

#include <memory>
#include <vector>

#include "beads.h"
#include "move_rewrite.h"
#include "multiloop.h"

// Combinatorial data for a shortening Reidemeister 2 move on a bigon.
//
// The face's phi cycle is (h0, h1): h0 runs from crossing X to crossing Y,
// h1 back from Y to X. With counterclockwise sigma cycles this forces
//
//     X = (h0, -h1, a, b)      Y = (h1, -h0, a', b')
//
// Straight-through pairs at X are {h0, a}, {-h1, b}; at Y {h1, a'},
// {-h0, b'}. Strand S1 enters on a, crosses along |h0| and leaves on b';
// strand S2 enters on b, crosses along |h1| and leaves on a'. The move
// splices a to b' and b to a'.
//
// Faces: L = face(a) and R = face(a') are the faces opposite the bigon at X
// and Y. The diagram disconnects exactly when L == R (that face is split
// between the two pieces); otherwise L, R and the bigon merge into one face.
struct R2FaceInfo
{
    int face = -1;
    int h0 = 0;
    int h1 = 0;
    int a = 0;
    int b = 0;
    int ap = 0; // a'
    int bp = 0; // b'
    int crossing[2] = {-1, -1};   // X, Y
    int vertexBead[2] = {-1, -1}; // X, Y
    int aTerminal = -1;
    int bTerminal = -1;
    int apTerminal = -1;
    int bpTerminal = -1;
    int faceL = -1;
    int faceR = -1;
    bool disconnects = false; // faceL == faceR
};

// Returns false unless the face is a bigon between two distinct crossings on
// two distinct edges. With a necklace, also fills vertex beads / terminals
// (false if the four terminals are not distinct edge beads).
bool AnalyzeR2Face(const Multiloop& loop, const Necklace* neck, int face, R2FaceInfo* out);

// When a disconnecting R2 nests an unloop in finite face L, the opposing
// host-face half-edge (the outer strand's XY-arc on L). 0 if the move does
// not nest through a finite face, or if `exteriorFace` is L itself.
int R2OpposingHostHalfEdge(const Multiloop& loop, const R2FaceInfo& info, int exteriorFace);

// Performs the move once each bigon side holds exactly one interior bead
// (m0 on |h0|, m1 on |h1|) and X, Y have been drawn tangent. New chains:
// aT - X - m1 - b'T and bT - m0 - Y - a'T. Face beads, nesting and the
// component list are updated through FinishMoveBookkeeping. Returns nullptr
// and leaves the necklace untouched on any inconsistency.
std::unique_ptr<Multiloop> FinalizeR2(const Multiloop& oldLoop, Necklace& neck,
                                      const R2FaceInfo& info,
                                      const std::vector<MultiloopComponent>& oldComponents,
                                      MoveComponentResult* out);
