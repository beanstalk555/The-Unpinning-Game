#include "physics_params.h"

void PhysicsParams::ResetDefaults()
{
    springStiffness = 500.0f;
    collisionStiffness = 250.0f;
    dragStiffness = 250.0f;
    damping = 80.0f;
    dt = 1.0f / 60.0f;
    substeps = 20;
    maxSpeed = 20.0f;
    tensionThreshold = 4.0;
    relaxStepsMax = 80;
    autoStepsPerFrame = 30;
    bendAngleEpsilon = 0.3f;
    minExcessArea = 0.538;
    ratioCombineMode = RatioCombineMode::Max;
    collisionBoundsScale = 3.0f;
    gapSize = 1.2f;
    r3FortifyStiffness = 500.0f;
    r3FaceBeadMass = 10.0f;
}
