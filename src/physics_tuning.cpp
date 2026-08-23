#include "physics_tuning.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

namespace
{
    const Color kPanelFill = {252, 252, 255, 250};
    const Color kTitleBar = {70, 110, 230, 255};
    const Color kBorder = {175, 177, 186, 255};
    const Color kTextDark = {40, 42, 48, 255};
    const Color kTextDim = {130, 132, 140, 255};
    const Color kAccent = {70, 110, 230, 255};

    constexpr float kPanelW = 400.0f;
    constexpr float kTitleH = 32.0f;
    constexpr float kRowH = 46.0f;
    constexpr float kPad = 14.0f;
    constexpr float kResetH = 30.0f;
    constexpr float kExtraRowH = 34.0f;
    constexpr int kSliderCount = 18;

    float Clampf(float v, float lo, float hi)
    {
        return v < lo ? lo : (v > hi ? hi : v);
    }

    float TextWidth(Font font, const std::string& text, float size)
    {
        return MeasureTextEx(font, text.c_str(), size, 0.0f).x;
    }

    void DrawTextS(Font font, const std::string& text, float x, float y, float size, Color color)
    {
        DrawTextEx(font, text.c_str(), Vector2{x, y}, size, 0.0f, color);
    }

    struct SliderSpec
    {
        const char* label;
        float minVal;
        float maxVal;
        bool isInt;
        bool isDouble;
    };

    const SliderSpec kSpecs[kSliderCount] = {
        {"Spring stiffness", 1.0f, 2000.0f, false, false},
        {"Collision stiffness", 1.0f, 500.0f, false, false},
        {"Drag stiffness", 1.0f, 500.0f, false, false},
        {"Damping", 0.0f, 100.0f, false, false},
        {"Timestep (1/dt)", 10.0f, 480.0f, false, false},
        {"Substeps", 1.0f, 100.0f, true, false},
        {"Max speed", 0.5f, 50.0f, false, false},
        {"Tension threshold", 0.1f, 100.0f, false, true},
        {"Relax steps max", 10.0f, 500.0f, true, false},
        {"Auto steps / frame", 1.0f, 120.0f, true, false},
        {"Bend angle epsilon", 0.0f, 0.6f, false, false},
        {"Min excess area", 0.0f, 5.0f, false, true},
        {"Collision bounds (x bead)", 1.0f, 8.0f, false, false},
        {"Gap size (x bead)", 1.0f, 1.5f, false, false},
        {"Honeycomb fortify stiffness", 1.0f, 2000.0f, false, false},
        {"R3 face bead mass", 1.0f, 1000.0f, false, false},
        {"R3 cycle repel stiffness", 0.0f, 20.0f, false, false},
        {"R3 post-repair gap size", 1.0f, 1.2f, false, false},
    };

    const char* kRatioModeLabels[] = {"minimum", "maximum", "average"};

    float ReadParam(const PhysicsParams& params, int index)
    {
        switch (index) {
        case 0: return params.springStiffness;
        case 1: return params.collisionStiffness;
        case 2: return params.dragStiffness;
        case 3: return params.damping;
        case 4: return 1.0f / params.dt;
        case 5: return static_cast<float>(params.substeps);
        case 6: return params.maxSpeed;
        case 7: return static_cast<float>(params.tensionThreshold);
        case 8: return static_cast<float>(params.relaxStepsMax);
        case 9: return static_cast<float>(params.autoStepsPerFrame);
        case 10: return params.bendAngleEpsilon;
        case 11: return static_cast<float>(params.minExcessArea);
        case 12: return params.collisionBoundsScale;
        case 13: return params.gapSize;
        case 14: return params.honeycombFortifyStiffness;
        case 15: return params.r3FaceBeadMass;
        case 16: return params.r3CycleRepelStiffness;
        case 17: return params.r3PostRepairGapSize;
        default: return 0.0f;
        }
    }

    void WriteParam(PhysicsParams& params, int index, float value)
    {
        switch (index) {
        case 0: params.springStiffness = value; break;
        case 1: params.collisionStiffness = value; break;
        case 2: params.dragStiffness = value; break;
        case 3: params.damping = value; break;
        case 4: params.dt = 1.0f / value; break;
        case 5: params.substeps = static_cast<int>(std::lround(value)); break;
        case 6: params.maxSpeed = value; break;
        case 7: params.tensionThreshold = value; break;
        case 8: params.relaxStepsMax = static_cast<int>(std::lround(value)); break;
        case 9: params.autoStepsPerFrame = static_cast<int>(std::lround(value)); break;
        case 10: params.bendAngleEpsilon = value; break;
        case 11: params.minExcessArea = value; break;
        case 12: params.collisionBoundsScale = value; break;
        case 13: params.gapSize = value; break;
        case 14: params.honeycombFortifyStiffness = value; break;
        case 15: params.r3FaceBeadMass = value; break;
        case 16: params.r3CycleRepelStiffness = value; break;
        case 17: params.r3PostRepairGapSize = value; break;
        default: break;
        }
    }

    std::string FormatValue(int index, float value)
    {
        char buf[32];
        if (kSpecs[index].isInt) {
            std::snprintf(buf, sizeof(buf), "%d", static_cast<int>(std::lround(value)));
        }
        else if (index == 4) {
            std::snprintf(buf, sizeof(buf), "%.0f Hz", value);
        }
        else if (kSpecs[index].isDouble || index == 3) {
            std::snprintf(buf, sizeof(buf), "%.3f", value);
        }
        else if (index == 17) {
            std::snprintf(buf, sizeof(buf), "%.2f", value);
        }
        else {
            std::snprintf(buf, sizeof(buf), "%.1f", value);
        }
        return buf;
    }

    constexpr float kScrollW = 12.0f;
    constexpr float kMinPanelH = kTitleH + 160.0f;

    float ContentHeight(const PhysicsTuningUI& ui)
    {
        float h = kPad + kSliderCount * kRowH + 2.0f * (20.0f + kExtraRowH) + kResetH +
                  kPad;
        if (ui.ratioDropdownOpen) {
            h += 3.0f * 26.0f;
        }
        return h;
    }

    Rectangle PanelRect(const PhysicsTuningUI& ui)
    {
        return Rectangle{ui.pos.x, ui.pos.y, kPanelW, ui.viewH};
    }
}

bool PhysicsTuningBlocksInput(const PhysicsTuningUI& ui, Vector2 mouse)
{
    return ui.open && CheckCollisionPointRec(mouse, PanelRect(ui));
}

void FramePhysicsTuning(PhysicsTuningUI& ui, PhysicsParams& params, Font font,
                        Vector2 mouse, bool click, int screenWidth, int screenHeight)
{
    if (!ui.open) {
        return;
    }

    if (!ui.posInitialized) {
        ui.pos = Vector2{static_cast<float>(screenWidth) - kPanelW - 16.0f, 120.0f};
        ui.posInitialized = true;
    }

    const float contentH = ContentHeight(ui);
    const float maxH = std::max(kMinPanelH, static_cast<float>(screenHeight) - 16.0f);
    ui.viewH = Clampf(kTitleH + contentH, kMinPanelH, maxH);

    ui.pos.x = Clampf(ui.pos.x, 0.0f, static_cast<float>(screenWidth) - kPanelW);
    ui.pos.y = Clampf(ui.pos.y, 0.0f, std::max(0.0f, static_cast<float>(screenHeight) - ui.viewH));

    const Rectangle panel = PanelRect(ui);
    const Rectangle titleBar = {panel.x, panel.y, panel.width, kTitleH};
    const Rectangle closeBtn = {panel.x + panel.width - 28.0f, panel.y + 6.0f, 22.0f, 20.0f};
    const Rectangle body = {panel.x, panel.y + kTitleH, panel.width, panel.height - kTitleH};
    const bool needsScroll = contentH > body.height + 0.5f;
    const float scrollMax = needsScroll ? contentH - body.height : 0.0f;
    const Rectangle scrollTrack = {panel.x + panel.width - kPad - kScrollW, body.y + 6.0f,
                                   kScrollW, std::max(12.0f, body.height - 12.0f)};

    if (click) {
        if (CheckCollisionPointRec(mouse, closeBtn)) {
            ui.open = false;
            ui.draggingWindow = false;
            ui.draggingScroll = false;
            ui.activeSlider = -1;
            ui.ratioDropdownOpen = false;
            return;
        }
        if (CheckCollisionPointRec(mouse, titleBar) && !CheckCollisionPointRec(mouse, closeBtn)) {
            ui.draggingWindow = true;
            ui.dragOffset = {mouse.x - ui.pos.x, mouse.y - ui.pos.y};
        }
        else if (needsScroll && CheckCollisionPointRec(mouse, scrollTrack)) {
            ui.draggingScroll = true;
            ui.activeSlider = -1;
        }
    }
    if (!IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
        ui.draggingWindow = false;
        ui.draggingScroll = false;
        ui.activeSlider = -1;
    }
    if (ui.draggingWindow) {
        ui.pos.x = mouse.x - ui.dragOffset.x;
        ui.pos.y = mouse.y - ui.dragOffset.y;
    }

    if (needsScroll && CheckCollisionPointRec(mouse, panel) && !ui.draggingWindow &&
        ui.activeSlider < 0) {
        ui.scrollY -= GetMouseWheelMove() * 40.0f;
    }
    if (ui.draggingScroll && scrollMax > 0.0f) {
        const float thumbH = std::max(24.0f, scrollTrack.height * (body.height / contentH));
        const float travel = std::max(1.0f, scrollTrack.height - thumbH);
        const float t = Clampf((mouse.y - scrollTrack.y - thumbH * 0.5f) / travel, 0.0f, 1.0f);
        ui.scrollY = t * scrollMax;
    }
    ui.scrollY = Clampf(ui.scrollY, 0.0f, scrollMax);

    DrawRectangle(static_cast<int>(panel.x) + 4, static_cast<int>(panel.y) + 4,
                  static_cast<int>(panel.width), static_cast<int>(panel.height),
                  Color{0, 0, 0, 35});
    DrawRectangleRec(panel, kPanelFill);
    DrawRectangleLinesEx(panel, 2.0f, kBorder);

    DrawRectangleRec(titleBar, kTitleBar);
    DrawTextS(font, "Physics parameters", titleBar.x + 12.0f, titleBar.y + 7.0f, 18.0f, WHITE);
    DrawRectangleRec(closeBtn, Color{255, 255, 255, 40});
    DrawTextS(font, "X", closeBtn.x + 6.0f, closeBtn.y + 1.0f, 16.0f, WHITE);

    const float trackX = panel.x + kPad;
    const float trackW = panel.width - 2.0f * kPad - (needsScroll ? (kScrollW + 6.0f) : 0.0f);
    const float labelSize = 15.0f;
    const float valueSize = 14.0f;
    const bool inBody = CheckCollisionPointRec(mouse, body) &&
                        !(needsScroll && CheckCollisionPointRec(mouse, scrollTrack));

    BeginScissorMode(static_cast<int>(body.x), static_cast<int>(body.y),
                     static_cast<int>(body.width), static_cast<int>(body.height));

    float y = body.y + kPad - ui.scrollY;

    for (int i = 0; i < kSliderCount; ++i) {
        const SliderSpec& spec = kSpecs[i];
        const float value = ReadParam(params, i);
        const std::string valueText = FormatValue(i, value);

        DrawTextS(font, spec.label, trackX, y, labelSize, kTextDark);
        const float vw = TextWidth(font, valueText, valueSize);
        DrawTextS(font, valueText, trackX + trackW - vw, y, valueSize, kTextDim);

        const float trackY = y + 22.0f;
        const Rectangle trackHit = {trackX - 4.0f, trackY - 10.0f, trackW + 8.0f, 24.0f};
        if (click && inBody && CheckCollisionPointRec(mouse, trackHit)) {
            ui.activeSlider = i;
            ui.ratioDropdownOpen = false;
        }
        if (ui.activeSlider == i && IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            const float t = Clampf((mouse.x - trackX) / trackW, 0.0f, 1.0f);
            const float newVal = spec.minVal + t * (spec.maxVal - spec.minVal);
            WriteParam(params, i, newVal);
        }

        const float t = (ReadParam(params, i) - spec.minVal) / (spec.maxVal - spec.minVal);
        DrawRectangle(static_cast<int>(trackX), static_cast<int>(trackY) - 2,
                      static_cast<int>(trackW), 4, kBorder);
        DrawRectangle(static_cast<int>(trackX), static_cast<int>(trackY) - 2,
                      static_cast<int>(trackW * t), 4, kAccent);
        DrawCircleV(Vector2{trackX + trackW * t, trackY}, 7.0f, kAccent);
        DrawCircleV(Vector2{trackX + trackW * t, trackY}, 3.5f, WHITE);

        y += kRowH;
    }

    {
        DrawTextS(font, "Edge excess combine:", trackX, y + 4.0f, labelSize, kTextDark);
        y += 20.0f;
        const int mode = static_cast<int>(params.ratioCombineMode);
        const std::string shown = kRatioModeLabels[mode];
        const Rectangle box = {trackX, y, trackW, 26.0f};
        if (click && inBody && CheckCollisionPointRec(mouse, box)) {
            ui.ratioDropdownOpen = !ui.ratioDropdownOpen;
        }
        DrawRectangleRec(box, WHITE);
        DrawRectangleLinesEx(box, 2.0f, ui.ratioDropdownOpen ? kAccent : kBorder);
        DrawTextS(font, shown, box.x + 8.0f, box.y + 4.0f, 15.0f, kTextDark);
        const float ax = box.x + box.width - 18.0f;
        const float ay = box.y + (ui.ratioDropdownOpen ? 16.0f : 10.0f);
        DrawTriangle(Vector2{ax - 5.0f, ay}, Vector2{ax + 5.0f, ay},
                     Vector2{ax, ay + (ui.ratioDropdownOpen ? -5.0f : 5.0f)}, kTextDim);
        y += kExtraRowH;

        if (ui.ratioDropdownOpen) {
            const float itemH = 26.0f;
            const Rectangle list = {box.x, box.y + box.height, box.width, 3.0f * itemH};
            DrawRectangleRec(list, WHITE);
            DrawRectangleLinesEx(list, 2.0f, kAccent);
            for (int m = 0; m < 3; ++m) {
                const Rectangle item = {list.x, list.y + m * itemH, list.width, itemH};
                const bool hover = inBody && CheckCollisionPointRec(mouse, item);
                if (m == mode) {
                    DrawRectangleRec(item, kAccent);
                }
                else if (hover) {
                    DrawRectangleRec(item, Color{222, 228, 246, 255});
                }
                DrawTextS(font, kRatioModeLabels[m], item.x + 8.0f, item.y + 4.0f, 15.0f,
                          m == mode ? WHITE : kTextDark);
                if (hover && click) {
                    params.ratioCombineMode = static_cast<RatioCombineMode>(m);
                    ui.ratioDropdownOpen = false;
                }
            }
        }
    }

    const Rectangle resetBtn = {trackX, y, trackW, kResetH};
    const bool resetHover = inBody && CheckCollisionPointRec(mouse, resetBtn);
    DrawRectangleRounded(resetBtn, 0.3f, 6,
                         resetHover ? ColorBrightness(kAccent, -0.12f) : kAccent);
    const std::string resetLabel = "Reset defaults";
    const float rw = TextWidth(font, resetLabel, 16.0f);
    DrawTextS(font, resetLabel, resetBtn.x + (resetBtn.width - rw) / 2.0f,
              resetBtn.y + 6.0f, 16.0f, WHITE);
    if (click && resetHover) {
        params.ResetDefaults();
    }

    EndScissorMode();

    if (needsScroll) {
        DrawRectangleRec(scrollTrack, Color{220, 222, 228, 255});
        const float thumbH = std::max(24.0f, scrollTrack.height * (body.height / contentH));
        const float travel = std::max(1.0f, scrollTrack.height - thumbH);
        const float thumbY = scrollTrack.y + travel * (scrollMax > 0.0f ? ui.scrollY / scrollMax : 0.0f);
        const Rectangle thumb = {scrollTrack.x + 1.0f, thumbY, scrollTrack.width - 2.0f, thumbH};
        DrawRectangleRounded(thumb, 0.4f, 4, ui.draggingScroll ? kAccent : Color{140, 144, 156, 255});
    }

    if (click && !CheckCollisionPointRec(mouse, panel)) {
        ui.ratioDropdownOpen = false;
    }
}
