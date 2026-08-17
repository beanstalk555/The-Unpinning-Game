#pragma once

// How to combine excess areas of the two faces bounding an edge.
enum class RatioCombineMode
{
    Min = 0,
    Max = 1,
    Avg = 2,
};

// Runtime-adjustable physics parameters (defaults match physics.cpp).
struct PhysicsParams
{
    float springStiffness = 500.0f;
    float collisionStiffness = 250.0f;
    float dragStiffness = 250.0f;
    float damping = 85.0f;
    float dt = 1.0f / 60.0f;
    int substeps = 20;
    float maxSpeed = 20.0f;

    double tensionThreshold = 3.323;
    int relaxStepsMax = 80;
    int autoStepsPerFrame = 30;
    // Beads with interior angle <= pi/3 + bendAngleEpsilon are removed first.
    float bendAngleEpsilon = 0.3f;

    // Edges qualify for the excess-area pass when the combined excess area of
    // their two bounding faces exceeds this (exterior face excluded).
    double minExcessArea = 0.538;
    RatioCombineMode ratioCombineMode = RatioCombineMode::Max;

    // Broad-phase bounding disk diameter as a multiple of the bead diameter.
    float collisionBoundsScale = 3.0f;

    // Center-to-center distance (in bead diameters) before inserting a bead while dragging.
    float gapSize = 1.2f;

    // Face-bead mass during Reidemeister 3 (edge beads stay mass 1).
    float r3FaceBeadMass = 10.0f;

    void ResetDefaults();
};
