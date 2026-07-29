#pragma once

#include <raylib.h>

#include "physics_params.h"

// Floating "Physics parameters" panel (raylib supports one OS window; this
// behaves as a separate draggable tool window drawn on top of the main view).
struct PhysicsTuningUI
{
    bool open = false;
    Vector2 pos = {0.0f, 0.0f};
    bool posInitialized = false;
    int activeSlider = -1;
    bool draggingWindow = false;
    Vector2 dragOffset = {0.0f, 0.0f};
    bool ratioDropdownOpen = false;
};

// True when the panel is open and the mouse is over it (blocks canvas input).
bool PhysicsTuningBlocksInput(const PhysicsTuningUI& ui, Vector2 mouse);

// Draw and handle the tuning panel. Call near the end of each frame.
void FramePhysicsTuning(PhysicsTuningUI& ui, PhysicsParams& params, Font font,
                        Vector2 mouse, bool click, int screenWidth, int screenHeight);
