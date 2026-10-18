#include "app.h"

#include <raylib.h>
#include <rlgl.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <deque>
#include <exception>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "beads.h"
#include "drawing.h"
#include "geometry.h"
#include "multiloop.h"
#include "orthogonal.h"
#include "permutation.h"
#include "physics.h"
#include "physics_params.h"
#include "physics_tuning.h"
#include "move_rewrite.h"
#include "r0_move.h"
#include "r1_move.h"
#include "r2_move.h"
#include "r3_move.h"
#include "random_multiloop.h"
#include "unpinning_game.h"

namespace
{
    // Raylib's high-DPI centering can park the window off-screen. Give the
    // restored window a size that fits, with the title bar on-screen, then
    // open maximized into the work area.
    void OpenWindowMaximized()
    {
        const int mon = GetCurrentMonitor();
        const Vector2 origin = GetMonitorPosition(mon);
        const int monX = static_cast<int>(origin.x);
        const int monY = static_cast<int>(origin.y);
        const int monW = GetMonitorWidth(mon);
        const int monH = GetMonitorHeight(mon);
        if (monW > 0 && monH > 0) {
            const int topMargin = 48;
            const int bottomMargin = 72;
            const int maxH = monH > topMargin + bottomMargin ? monH - topMargin - bottomMargin : monH;
            int winW = 1460;
            int winH = 850;
            if (winW > monW) {
                winW = monW;
            }
            if (winH > maxH) {
                winH = maxH;
            }
            SetWindowSize(winW, winH);
            const int x = monX + (monW - winW) / 2;
            const int y = monY + topMargin + (maxH > winH ? (maxH - winH) / 2 : 0);
            SetWindowPosition(x, y);
        }
        MaximizeWindow();
    }

    const Color kBackground = {245, 245, 248, 255};
    const Color kHighlightFill = {255, 220, 40, 90};
    const Color kTextDark = {40, 42, 48, 255};
    const Color kTextDim = {130, 132, 140, 255};
    const Color kAccent = {70, 110, 230, 255};
    const Color kAccentDark = {55, 90, 200, 255};
    const Color kSelection = {70, 110, 230, 70};
    const Color kBorder = {175, 177, 186, 255};
    const Color kErrorRed = {200, 40, 50, 255};
    const Color kHover = {222, 228, 246, 255};
    const Color kBeadOutline = {45, 45, 50, 255};
    const Color kFaceBeadFill = {88, 162, 70, 255};
    const Color kFaceBeadEdge = {50, 105, 40, 255};
    const Color kGoldInner = {255, 228, 120, 165};
    const Color kGoldOuter = {176, 120, 18, 130};
    const Color kSilverInner = {228, 232, 238, 165};
    const Color kSilverOuter = {118, 124, 136, 125};
    const Color kBronzeInner = {214, 154, 88, 165};
    const Color kBronzeOuter = {138, 72, 28, 130};

    // Distinct strand colors: 8 hand-picked ones, then a golden-angle walk
    // around the hue circle with varying saturation/value (distinct for 50+).
    Color StrandColor(int strand)
    {
        static const Color base[] = {
            {204, 36, 50, 255},   // red
            {38, 70, 200, 255},   // blue
            {22, 140, 90, 255},   // green
            {235, 140, 25, 255},  // orange
            {140, 48, 190, 255},  // purple
            {110, 75, 38, 255},   // brown
            {0, 150, 160, 255},   // teal
            {225, 70, 165, 255},  // pink
        };
        constexpr int kNumBase = static_cast<int>(sizeof(base) / sizeof(base[0]));
        if (strand < kNumBase) {
            return base[strand];
        }
        const int i = strand - kNumBase;
        const float hue = std::fmod(i * 137.508f + 24.0f, 360.0f);
        const float saturation = ((i / 3) % 2 == 0) ? 0.85f : 0.55f;
        const float value = 0.50f + 0.13f * (i % 3);
        return ColorFromHSV(hue, saturation, value);
    }

    // ----------------------------------------------------------------
    // Font helpers.
    // ----------------------------------------------------------------

    Font gFont;
    Font gFontBold;
    bool gFontLoaded = false;
    bool gFontBoldLoaded = false;

    // ASCII plus the infinity sign used by the move-allowance dropdown.
    void FillUiCodepoints(int* codepoints, int* count)
    {
        int n = 0;
        for (int cp = 32; cp <= 126; ++cp) {
            codepoints[n++] = cp;
        }
        codepoints[n++] = 0x221E;
        *count = n;
    }

    void LoadUiFont()
    {
        const char* candidates[] = {
            "C:/Windows/Fonts/segoeui.ttf", // Segoe UI: friendly, modern, sans-serif
            "C:/Windows/Fonts/calibri.ttf",
            "C:/Windows/Fonts/arial.ttf",
            "/System/Library/Fonts/Supplemental/Arial.ttf", // macOS
            "/Library/Fonts/Arial.ttf",
            "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf", // Linux, if installed
            "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf", // Debian/Ubuntu
            "/usr/share/fonts/liberation-sans/LiberationSans-Regular.ttf", // Fedora
            "/usr/share/fonts/TTF/LiberationSans-Regular.ttf", // Arch
            "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
            "/usr/share/fonts/dejavu-sans-fonts/DejaVuSans.ttf",
            "/usr/share/fonts/dejavu/DejaVuSans.ttf",
            "/usr/share/fonts/TTF/DejaVuSans.ttf",
            "/usr/share/fonts/truetype/noto/NotoSans-Regular.ttf",
            "/usr/share/fonts/noto/NotoSans-Regular.ttf",
            "/usr/share/fonts/google-noto/NotoSans-Regular.ttf",
            "/usr/share/fonts/truetype/ubuntu/Ubuntu-R.ttf",
            "/usr/share/fonts/truetype/freefont/FreeSans.ttf",
        };
        const char* boldCandidates[] = {
            "C:/Windows/Fonts/segoeuib.ttf",
            "C:/Windows/Fonts/calibrib.ttf",
            "C:/Windows/Fonts/arialbd.ttf",
            "/System/Library/Fonts/Supplemental/Arial Bold.ttf", // macOS
            "/Library/Fonts/Arial Bold.ttf",
            "/usr/share/fonts/truetype/msttcorefonts/Arial_Bold.ttf",
            "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
            "/usr/share/fonts/liberation-sans/LiberationSans-Bold.ttf",
            "/usr/share/fonts/TTF/LiberationSans-Bold.ttf",
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
            "/usr/share/fonts/dejavu-sans-fonts/DejaVuSans-Bold.ttf",
            "/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf",
            "/usr/share/fonts/TTF/DejaVuSans-Bold.ttf",
            "/usr/share/fonts/truetype/noto/NotoSans-Bold.ttf",
            "/usr/share/fonts/noto/NotoSans-Bold.ttf",
            "/usr/share/fonts/google-noto/NotoSans-Bold.ttf",
            "/usr/share/fonts/truetype/ubuntu/Ubuntu-B.ttf",
            "/usr/share/fonts/truetype/freefont/FreeSansBold.ttf",
        };
        int codepoints[128];
        int codepointCount = 0;
        FillUiCodepoints(codepoints, &codepointCount);
        for (const char* path : candidates) {
            if (FileExists(path)) {
                gFont = LoadFontEx(path, 64, codepoints, codepointCount);
                if (gFont.texture.id != 0) {
                    SetTextureFilter(gFont.texture, TEXTURE_FILTER_BILINEAR);
                    gFontLoaded = true;
                    break;
                }
            }
        }
        if (!gFontLoaded) {
            gFont = GetFontDefault();
        }
        for (const char* path : boldCandidates) {
            if (FileExists(path)) {
                gFontBold = LoadFontEx(path, 64, codepoints, codepointCount);
                if (gFontBold.texture.id != 0) {
                    SetTextureFilter(gFontBold.texture, TEXTURE_FILTER_BILINEAR);
                    gFontBoldLoaded = true;
                    break;
                }
            }
        }
        if (!gFontBoldLoaded) {
            gFontBold = gFont;
        }
    }

    void UnloadUiFont()
    {
        if (gFontBoldLoaded) {
            UnloadFont(gFontBold);
            gFontBoldLoaded = false;
        }
        if (gFontLoaded) {
            UnloadFont(gFont);
            gFontLoaded = false;
        }
    }

    float TextWidthFont(const Font& font, const std::string& text, float size)
    {
        return MeasureTextEx(font, text.c_str(), size, 0.0f).x;
    }

    float TextWidth(const std::string& text, float size)
    {
        return TextWidthFont(gFont, text, size);
    }

    void DrawTextFont(const Font& font, const std::string& text, float x, float y, float size,
                      Color color)
    {
        DrawTextEx(font, text.c_str(), Vector2{x, y}, size, 0.0f, color);
    }

    std::vector<std::string> BuildPhysicsTimingLines(const PhysicsDebugTimings& t, bool autoTightening)
    {
        std::vector<std::string> lines;
        char buf[96];
        std::snprintf(buf, sizeof(buf), "Physics timing (%.1f ms frame):", t.frameTotalMs);
        lines.push_back(buf);
        if (autoTightening) {
            std::snprintf(buf, sizeof(buf), "  auto-step: %.2f ms", t.autoStepOverheadMs);
        }
        else {
            std::snprintf(buf, sizeof(buf), "  step: %.2f ms", t.stepMs);
        }
        lines.push_back(buf);
        std::snprintf(buf, sizeof(buf), "  steps: %d  substeps: %d", t.stepCount, t.substepsTotal);
        lines.push_back(buf);
        std::snprintf(buf, sizeof(buf), "  integrate: %.2f ms", t.integrateMs);
        lines.push_back(buf);
        std::snprintf(buf, sizeof(buf), "  collisions: %.2f ms", t.collisionMs);
        lines.push_back(buf);
        if (t.beginRemovalRoundMs > 0.001) {
            std::snprintf(buf, sizeof(buf), "  removal round: %.2f ms", t.beginRemovalRoundMs);
            lines.push_back(buf);
        }
        if (t.strikeHoneycombMs > 0.001) {
            std::snprintf(buf, sizeof(buf), "    honeycomb walk: %.2f ms", t.strikeHoneycombMs);
            lines.push_back(buf);
        }
        if (t.rebuildStructuresMs > 0.001) {
            std::snprintf(buf, sizeof(buf), "  rebuild structures: %.2f ms", t.rebuildStructuresMs);
            lines.push_back(buf);
        }
        if (t.rebuildFaceGeometryMs > 0.001) {
            std::snprintf(buf, sizeof(buf), "    face geometry: %.2f ms", t.rebuildFaceGeometryMs);
            lines.push_back(buf);
        }
        if (t.rebuildChainsMs > 0.001) {
            std::snprintf(buf, sizeof(buf), "  rebuild chains: %.2f ms", t.rebuildChainsMs);
            lines.push_back(buf);
        }
        if (t.rebuildFaceBodiesMs > 0.001) {
            std::snprintf(buf, sizeof(buf), "  rebuild face bodies: %.2f ms", t.rebuildFaceBodiesMs);
            lines.push_back(buf);
        }
        return lines;
    }

    void DrawTextS(const std::string& text, float x, float y, float size, Color color)
    {
        DrawTextFont(gFont, text, x, y, size, color);
    }

    void DrawDebugTextPanel(float x, float y, const std::vector<std::string>& lines, float fontSize,
                            float lineH, float pad, Color fill, float maxHeight = 0.0f)
    {
        float boxW = 0.0f;
        for (const std::string& line : lines) {
            boxW = std::max(boxW, TextWidth(line, fontSize));
        }
        boxW += 2.0f * pad;
        int visibleLines = static_cast<int>(lines.size());
        float boxH = static_cast<float>(visibleLines) * lineH + 2.0f * pad;
        if (maxHeight > 0.0f && boxH > maxHeight) {
            const float inner = maxHeight - 2.0f * pad;
            visibleLines = std::max(1, static_cast<int>(inner / lineH));
            boxH = static_cast<float>(visibleLines) * lineH + 2.0f * pad;
        }
        DrawRectangleRec({x, y, boxW, boxH}, fill);
        DrawRectangleLinesEx({x, y, boxW, boxH}, 2.0f, kAccent);
        BeginScissorMode(static_cast<int>(x), static_cast<int>(y), static_cast<int>(boxW),
                         static_cast<int>(boxH));
        float ty = y + pad;
        for (int i = 0; i < visibleLines && i < static_cast<int>(lines.size()); ++i) {
            DrawTextS(lines[static_cast<size_t>(i)], x + pad, ty, fontSize, kTextDark);
            ty += lineH;
        }
        EndScissorMode();
    }

    float DebugPanelHeight(int lineCount, float lineH, float pad, float maxHeight = 0.0f)
    {
        float h = static_cast<float>(lineCount) * lineH + 2.0f * pad;
        if (maxHeight > 0.0f) {
            h = std::min(h, maxHeight);
        }
        return h;
    }

    // ----------------------------------------------------------------
    // Small geometry helpers (screen space).
    // ----------------------------------------------------------------

    float Clampf(float v, float lo, float hi)
    {
        return v < lo ? lo : (v > hi ? hi : v);
    }

    float Dist(Vector2 a, Vector2 b)
    {
        const float dx = a.x - b.x;
        const float dy = a.y - b.y;
        return std::sqrt(dx * dx + dy * dy);
    }

    // The point on the boundary of the bead at `center` toward `other`.
    Vector2 TangentPoint(Vector2 center, float radius, Vector2 other)
    {
        const float d = Dist(center, other);
        if (d < 1e-6f) {
            return center;
        }
        return Vector2{center.x + (other.x - center.x) * radius / d,
                       center.y + (other.y - center.y) * radius / d};
    }

    // Draws the arc inside the bead (center, R) joining boundary points a and
    // b, meeting the boundary at right angles: the arc of the unique circle
    // through a and b orthogonal to the bead, or the chord when a and b are
    // (nearly) antipodal.
    void DrawOrthogonalArc(Vector2 center, float radius, Vector2 a, Vector2 b,
                           float thickness, Color color)
    {
        const Vector2 mid = {(a.x + b.x) * 0.5f, (a.y + b.y) * 0.5f};
        const float midDist = Dist(mid, center);
        if (midDist < radius * 0.04f) {
            DrawLineEx(a, b, thickness, color);
            DrawCircleV(a, thickness / 2.0f, color);
            DrawCircleV(b, thickness / 2.0f, color);
            return;
        }

        // The orthogonal circle's center m lies on the ray from the bead
        // center through the chord midpoint, at distance R^2 / |mid - c|.
        const Vector2 dir = {(mid.x - center.x) / midDist, (mid.y - center.y) / midDist};
        const float t = radius * radius / midDist;
        const Vector2 m = {center.x + dir.x * t, center.y + dir.y * t};
        const float rho = Dist(m, a);

        const float ang0 = std::atan2(a.y - m.y, a.x - m.x);
        const float ang1 = std::atan2(b.y - m.y, b.x - m.x);
        const float twoPi = 6.2831853f;
        float sweep = std::fmod(ang1 - ang0 + 2.0f * twoPi, twoPi);
        // Pick the arc that stays inside the bead.
        const float midAngle = ang0 + sweep * 0.5f;
        const Vector2 arcMid = {m.x + rho * std::cos(midAngle), m.y + rho * std::sin(midAngle)};
        if (Dist(arcMid, center) > radius) {
            sweep -= twoPi;
        }

        const int steps = std::max(8, static_cast<int>(std::fabs(sweep) * 12.0f));
        Vector2 prev = a;
        DrawCircleV(a, thickness / 2.0f, color);
        for (int i = 1; i <= steps; ++i) {
            const float ang = ang0 + sweep * static_cast<float>(i) / static_cast<float>(steps);
            const Vector2 p = {m.x + rho * std::cos(ang), m.y + rho * std::sin(ang)};
            DrawLineEx(prev, p, thickness, color);
            DrawCircleV(p, thickness / 2.0f, color);
            prev = p;
        }
    }

    std::string CycleToString(const std::vector<int>& cycle)
    {
        std::string out = "(";
        for (size_t i = 0; i < cycle.size(); ++i) {
            if (i > 0) {
                out += ',';
            }
            out += std::to_string(cycle[i]);
        }
        out += ')';
        return out;
    }

    // ----------------------------------------------------------------
    // Text box widget.
    // ----------------------------------------------------------------

    struct TextBoxState
    {
        std::string text;
        int caret = 0;
        int anchor = 0;
        bool focused = false;
        bool selecting = false;
        float scrollX = 0.0f;
        double lastClickTime = -1.0;
    };

    // Handles input and drawing for one text box. Returns true when Enter is
    // pressed while the box is focused. When disabled, the box is drawn grayed
    // and ignores all input.
    bool TextBoxFrame(TextBoxState& tb, Rectangle box, float fontSize, Vector2 mouse,
                      bool click, bool digitsOnly, bool enabled = true)
    {
        if (!enabled) {
            tb.focused = false;
            tb.selecting = false;
            DrawRectangleRec(box, Color{225, 226, 230, 255});
            DrawRectangleLinesEx(box, 2.0f, kBorder);
            BeginScissorMode(static_cast<int>(box.x) + 4, static_cast<int>(box.y) + 2,
                             static_cast<int>(box.width) - 8, static_cast<int>(box.height) - 4);
            DrawTextS(tb.text, box.x + 8.0f, box.y + (box.height - fontSize) / 2.0f - 1.0f,
                      fontSize, kTextDim);
            EndScissorMode();
            return false;
        }
        const bool ctrl = IsKeyDown(KEY_LEFT_CONTROL) || IsKeyDown(KEY_RIGHT_CONTROL);
        const bool shift = IsKeyDown(KEY_LEFT_SHIFT) || IsKeyDown(KEY_RIGHT_SHIFT);
        const auto pressed = [](int key) { return IsKeyPressed(key) || IsKeyPressedRepeat(key); };
        const auto accepts = [&](char c) {
            return digitsOnly ? (c >= '0' && c <= '9') : (c >= 32 && c < 127);
        };
        const auto computePrefix = [&](std::vector<float>& prefix) {
            prefix.assign(tb.text.size() + 1, 0.0f);
            for (size_t i = 0; i < tb.text.size(); ++i) {
                prefix[i + 1] = prefix[i] + TextWidth(std::string(1, tb.text[i]), fontSize);
            }
        };
        std::vector<float> prefix;
        computePrefix(prefix);

        const auto caretFromMouse = [&]() {
            const float relX = mouse.x - (box.x + 8.0f - tb.scrollX);
            int best = 0;
            float bestDist = 1e9f;
            for (size_t i = 0; i < prefix.size(); ++i) {
                const float d = std::fabs(prefix[i] - relX);
                if (d < bestDist) {
                    bestDist = d;
                    best = static_cast<int>(i);
                }
            }
            return best;
        };

        // ---- input ----
        if (click) {
            if (CheckCollisionPointRec(mouse, box)) {
                const double now = GetTime();
                tb.focused = true;
                if (now - tb.lastClickTime < 0.35) {
                    // Double click: select all.
                    tb.anchor = 0;
                    tb.caret = static_cast<int>(tb.text.size());
                    tb.selecting = false;
                }
                else {
                    tb.selecting = true;
                    tb.caret = caretFromMouse();
                    tb.anchor = tb.caret;
                }
                tb.lastClickTime = now;
            }
            else {
                tb.focused = false;
            }
        }
        if (tb.selecting) {
            if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
                tb.caret = caretFromMouse();
            }
            else {
                tb.selecting = false;
            }
        }

        bool enterPressed = false;
        if (tb.focused) {
            const auto eraseSelection = [&]() {
                const int lo = std::min(tb.caret, tb.anchor);
                const int hi = std::max(tb.caret, tb.anchor);
                if (lo != hi) {
                    tb.text.erase(lo, hi - lo);
                    tb.caret = lo;
                    tb.anchor = lo;
                }
            };
            const bool hasSelection = tb.caret != tb.anchor;

            int c;
            while ((c = GetCharPressed()) != 0) {
                if (accepts(static_cast<char>(c)) && !ctrl) {
                    eraseSelection();
                    tb.text.insert(tb.text.begin() + tb.caret, static_cast<char>(c));
                    ++tb.caret;
                    tb.anchor = tb.caret;
                }
            }
            if (ctrl && IsKeyPressed(KEY_A)) {
                tb.anchor = 0;
                tb.caret = static_cast<int>(tb.text.size());
            }
            if (ctrl && IsKeyPressed(KEY_C) && hasSelection) {
                const int lo = std::min(tb.caret, tb.anchor);
                const int hi = std::max(tb.caret, tb.anchor);
                SetClipboardText(tb.text.substr(lo, hi - lo).c_str());
            }
            if (ctrl && IsKeyPressed(KEY_V)) {
                const char* clip = GetClipboardText();
                if (clip != nullptr) {
                    eraseSelection();
                    std::string cleaned;
                    for (const char* p = clip; *p != '\0'; ++p) {
                        if (accepts(*p)) {
                            cleaned.push_back(*p);
                        }
                    }
                    tb.text.insert(tb.caret, cleaned);
                    tb.caret += static_cast<int>(cleaned.size());
                    tb.anchor = tb.caret;
                }
            }
            if (pressed(KEY_BACKSPACE)) {
                if (tb.caret != tb.anchor) {
                    eraseSelection();
                }
                else if (tb.caret > 0) {
                    tb.text.erase(tb.caret - 1, 1);
                    --tb.caret;
                    tb.anchor = tb.caret;
                }
            }
            if (pressed(KEY_DELETE)) {
                if (tb.caret != tb.anchor) {
                    eraseSelection();
                }
                else if (tb.caret < static_cast<int>(tb.text.size())) {
                    tb.text.erase(tb.caret, 1);
                }
            }
            if (pressed(KEY_LEFT)) {
                tb.caret = std::max(0, tb.caret - 1);
                if (!shift) {
                    tb.anchor = tb.caret;
                }
            }
            if (pressed(KEY_RIGHT)) {
                tb.caret = std::min(static_cast<int>(tb.text.size()), tb.caret + 1);
                if (!shift) {
                    tb.anchor = tb.caret;
                }
            }
            if (IsKeyPressed(KEY_HOME)) {
                tb.caret = 0;
                if (!shift) {
                    tb.anchor = tb.caret;
                }
            }
            if (IsKeyPressed(KEY_END)) {
                tb.caret = static_cast<int>(tb.text.size());
                if (!shift) {
                    tb.anchor = tb.caret;
                }
            }
            if (IsKeyPressed(KEY_ENTER) || IsKeyPressed(KEY_KP_ENTER)) {
                enterPressed = true;
            }
            if (IsKeyPressed(KEY_ESCAPE)) {
                tb.focused = false;
            }

            // Text may have changed: recompute widths, keep the caret in view.
            computePrefix(prefix);
            tb.caret = std::min(tb.caret, static_cast<int>(tb.text.size()));
            tb.anchor = std::min(tb.anchor, static_cast<int>(tb.text.size()));
            const float innerWidth = box.width - 16.0f;
            const float caretX = prefix[tb.caret];
            if (caretX - tb.scrollX < 0.0f) {
                tb.scrollX = caretX;
            }
            if (caretX - tb.scrollX > innerWidth) {
                tb.scrollX = caretX - innerWidth;
            }
            tb.scrollX = Clampf(tb.scrollX, 0.0f, std::max(0.0f, prefix.back() - innerWidth));
        }

        // ---- draw ----
        DrawRectangleRec(box, WHITE);
        DrawRectangleLinesEx(box, 2.0f, tb.focused ? kAccent : kBorder);
        BeginScissorMode(static_cast<int>(box.x) + 4, static_cast<int>(box.y) + 2,
                         static_cast<int>(box.width) - 8, static_cast<int>(box.height) - 4);
        const float textX = box.x + 8.0f - tb.scrollX;
        const float textY = box.y + (box.height - fontSize) / 2.0f - 1.0f;
        if (tb.caret != tb.anchor) {
            const int lo = std::min(tb.caret, tb.anchor);
            const int hi = std::max(tb.caret, tb.anchor);
            DrawRectangle(static_cast<int>(textX + prefix[lo]), static_cast<int>(textY) - 2,
                          static_cast<int>(prefix[hi] - prefix[lo]),
                          static_cast<int>(fontSize) + 5, kSelection);
        }
        DrawTextS(tb.text, textX, textY, fontSize, kTextDark);
        if (tb.focused && (tb.selecting || static_cast<int>(GetTime() * 2.0) % 2 == 0)) {
            DrawRectangle(static_cast<int>(textX + prefix[tb.caret]), static_cast<int>(textY) - 2,
                          2, static_cast<int>(fontSize) + 5, kTextDark);
        }
        EndScissorMode();
        return enterPressed;
    }

    // ----------------------------------------------------------------
    // Application state.
    // ----------------------------------------------------------------

    struct App
    {
        TextBoxState sigmaBox;
        TextBoxState countBox;
        float density = 0.8f;
        bool sliderDragging = false;
        bool forceOneStrand = false;
        bool forceSimple = false;
        bool dropdownOpen = false;
        float dropdownScroll = 0.0f;
        std::mt19937 rng{std::random_device{}()};

        std::string error;
        std::unique_ptr<Multiloop> loop;
        std::vector<MultiloopComponent> components;
        std::vector<std::string> faceLabels;
        int selectedFace = -1;
        bool exteriorHighlighted = false;
        std::unique_ptr<GridLayout> layout;
        std::unique_ptr<Necklace> necklace;
        std::unique_ptr<Physics> physics;
        PhysicsParams physicsParams;
        PhysicsTuningUI tuningUI;
        bool debug = false;
        bool showBeads = false;
        bool highlightReidemeister = false;
        bool showDebugTiming = true;
        bool showDebugExcessAreas = false;
        bool showGlobalCollisions = false;
        int debugFaceBeadIndex = -1;
        int pickClickBody = -1;
        Vector2 pickClickPos = {0.0f, 0.0f};
        bool pickPressedOnCanvas = false;
        // The 1-face case (a simple closed curve) has no crossings, so it
        // bypasses the permutation pipeline and is drawn directly.
        bool unknotMode = false;
        bool wasAutoTightening = false;
        bool wasR3Active = false;

        // Camera in grid coordinates: `viewCell` pixels per grid unit, look-at
        // is the grid point mapped to the canvas center, `viewAngle` is
        // counterclockwise about that point (y up). Animated on fit.
        const void* viewLayoutKey = nullptr;
        bool viewReady = false;
        float viewCell = 1.0f;
        float viewLookX = 0.0f;
        float viewLookY = 0.0f;
        float viewAngle = 0.0f;
        int viewGridHeight = 0;
        bool viewAnimating = false;
        float viewAnimT = 0.0f;
        float viewCell0 = 1.0f;
        float viewCell1 = 1.0f;
        float viewLookX0 = 0.0f;
        float viewLookX1 = 0.0f;
        float viewLookY0 = 0.0f;
        float viewLookY1 = 0.0f;
        float viewAngle0 = 0.0f;
        float viewAngle1 = 0.0f;
        bool viewFitPending = false;

        void RequestViewFit() { viewFitPending = true; }

        static constexpr int kR3UndoLimit = 10;
        struct R3UndoSnapshot
        {
            std::unique_ptr<Multiloop> loop;
            std::vector<MultiloopComponent> components;
            GridLayout layout;
            Necklace necklace;
            bool unknotMode = false;
            int selectedFace = -1;
            bool exteriorHighlighted = false;
            float viewCell = 1.0f;
            float viewLookX = 0.0f;
            float viewLookY = 0.0f;
            float viewAngle = 0.0f;
            int viewGridHeight = 0;
        };
        std::deque<R3UndoSnapshot> r3Undo;

        void ClearR3Undo() { r3Undo.clear(); }

        R3UndoSnapshot CaptureSnapshot() const
        {
            R3UndoSnapshot snap;
            if (!layout || !necklace) {
                return snap;
            }
            if (loop) {
                snap.loop.reset(new Multiloop(*loop));
            }
            snap.components = components;
            snap.layout = *layout;
            snap.necklace = *necklace;
            snap.unknotMode = unknotMode;
            snap.selectedFace = selectedFace;
            snap.exteriorHighlighted = exteriorHighlighted;
            snap.viewCell = viewCell;
            snap.viewLookX = viewLookX;
            snap.viewLookY = viewLookY;
            snap.viewAngle = viewAngle;
            snap.viewGridHeight = viewGridHeight;
            return snap;
        }

        static R3UndoSnapshot CloneSnapshot(const R3UndoSnapshot& src)
        {
            R3UndoSnapshot dst;
            if (src.loop) {
                dst.loop.reset(new Multiloop(*src.loop));
            }
            dst.components = src.components;
            dst.layout = src.layout;
            dst.necklace = src.necklace;
            dst.unknotMode = src.unknotMode;
            dst.selectedFace = src.selectedFace;
            dst.exteriorHighlighted = src.exteriorHighlighted;
            dst.viewCell = src.viewCell;
            dst.viewLookX = src.viewLookX;
            dst.viewLookY = src.viewLookY;
            dst.viewAngle = src.viewAngle;
            dst.viewGridHeight = src.viewGridHeight;
            return dst;
        }

        void ApplySnapshot(R3UndoSnapshot snap)
        {
            if (physics && physics->AutoTightening()) {
                physics->StopAutoTighten();
            }
            if (!snap.unknotMode && !snap.loop) {
                return;
            }
            physics.reset();
            loop = std::move(snap.loop);
            components = std::move(snap.components);
            layout.reset(new GridLayout(std::move(snap.layout)));
            necklace.reset(new Necklace(std::move(snap.necklace)));
            unknotMode = snap.unknotMode;
            selectedFace = snap.selectedFace;
            exteriorHighlighted = snap.exteriorHighlighted;
            debugFaceBeadIndex = -1;
            dropdownOpen = false;
            error.clear();
            wasAutoTightening = false;
            wasR3Active = false;
            faceLabels.clear();
            if (unknotMode) {
                if (CircleCount() <= 1) {
                    faceLabels = {"()", "()'"};
                }
            }
            else if (loop) {
                for (const std::vector<int>& face : loop->faces) {
                    faceLabels.push_back(CycleToString(face));
                }
            }
            SyncSigmaText();
            viewLayoutKey = layout.get();
            viewReady = true;
            viewCell = snap.viewCell;
            viewLookX = snap.viewLookX;
            viewLookY = snap.viewLookY;
            viewAngle = snap.viewAngle;
            viewGridHeight = snap.viewGridHeight;
            viewAnimating = false;
            try {
                if (unknotMode) {
                    physics.reset(new Physics(*necklace, physicsParams));
                }
                else if (loop) {
                    physics.reset(new Physics(*loop, *necklace, physicsParams, PhysicsExteriors()));
                }
                if (physics) {
                    physics->SetNesting(components);
                }
            }
            catch (const std::exception& e) {
                physics.reset();
                error = e.what();
            }
        }

        void PushR3Undo()
        {
            if (!layout || !necklace) {
                return;
            }
            r3Undo.push_back(CaptureSnapshot());
            while (static_cast<int>(r3Undo.size()) > kR3UndoLimit) {
                r3Undo.pop_front();
            }
        }

        void UndoR3()
        {
            if (r3Undo.empty()) {
                return;
            }
            R3UndoSnapshot snap = std::move(r3Undo.back());
            r3Undo.pop_back();
            ApplySnapshot(std::move(snap));
        }

        R3UndoSnapshot gameStart_;
        bool gameStartValid_ = false;
        bool pendingGameFit = false;
        UnpinningGame game;

        void RememberGameStart()
        {
            gameStart_ = CaptureSnapshot();
            gameStartValid_ = layout && necklace;
        }

        void RestoreGameStart()
        {
            if (!gameStartValid_) {
                return;
            }
            ApplySnapshot(CloneSnapshot(gameStart_));
        }

        // Free-hand drawing mode.
        bool drawing = false;
        std::vector<std::vector<Vector2>> drawCurves; // completed closed curves
        std::vector<Vector2> drawCurrent;             // in-progress open polyline
        std::vector<Vector2> drawCrossings;           // crossings of committed geometry
        std::string drawWarning;
        Vector2 drawWarningPos = {0.0f, 0.0f};
        double drawWarningUntil = 0.0;
        double lastCanvasClick = -1.0;
        // Snapshot for Cancel.
        std::string savedSigmaText;
        int savedFace = -1;
        bool savedUnknot = false;

        App()
        {
            countBox.text = "10";
            StartRandomOddMatch();
        }

        void StartRandomOddMatch()
        {
            RandomMultiloop random;
            try {
                random = MakeOddOneStrandMultiloop(rng());
            }
            catch (const std::exception& e) {
                error = e.what();
                return;
            }
            sigmaBox.text = random.sigma.ToString();
            SubmitSigma();
            if (!loop) {
                return;
            }
            SelectFaceContaining(random.exteriorHalfEdge);
            if (!layout || !necklace || !physics) {
                return;
            }
            PinEveryFaceBead();
            BeginUnpinningMatch();
            if (!physics->MoveActive() && !physics->AutoTightening()) {
                physics->StartAutoTighten();
            }
        }

        std::unordered_set<int> PhysicsExteriors() const
        {
            std::unordered_set<int> faces;
            for (const MultiloopComponent& component : components) {
                if (component.kind == MultiloopPart::Crossings && component.exteriorFace >= 0) {
                    faces.insert(component.exteriorFace);
                }
            }
            if (faces.empty() && selectedFace >= 0) {
                faces.insert(selectedFace);
            }
            return faces;
        }

        struct HighlightSnapshot
        {
            std::unordered_map<int, bool> byFace;
            std::unordered_map<int, bool> byCircle;
            std::unordered_map<int, bool> pinnedByFace;
            std::unordered_map<int, bool> pinnedByCircle;
            bool exterior = false;
        };

        static int CirclePadStrand(const Necklace& neck, const FaceBead& fb)
        {
            for (const int bead : fb.boundaryBeads) {
                if (bead >= 0 && bead < static_cast<int>(neck.beads.size()) &&
                    !neck.beads[static_cast<size_t>(bead)].removed) {
                    return neck.beads[static_cast<size_t>(bead)].strand;
                }
            }
            return -1;
        }

        HighlightSnapshot CaptureHighlights() const
        {
            HighlightSnapshot snap;
            snap.exterior = exteriorHighlighted;
            if (!necklace) {
                return snap;
            }
            for (const FaceBead& fb : necklace->faceBeads) {
                if (fb.face >= 0) {
                    snap.byFace[fb.face] = snap.byFace[fb.face] || fb.highlighted;
                    snap.pinnedByFace[fb.face] = snap.pinnedByFace[fb.face] || fb.pinned;
                }
                else {
                    const int strand = CirclePadStrand(*necklace, fb);
                    if (strand >= 0) {
                        snap.byCircle[strand] = snap.byCircle[strand] || fb.highlighted;
                        snap.pinnedByCircle[strand] =
                            snap.pinnedByCircle[strand] || fb.pinned;
                    }
                }
            }
            return snap;
        }

        void RestoreHighlights(const HighlightSnapshot& snap, const std::unordered_set<int>& oldExt,
                               const std::unordered_set<int>& newExt)
        {
            if (necklace) {
                for (FaceBead& fb : necklace->faceBeads) {
                    fb.highlighted = false;
                    fb.pinned = false;
                    if (fb.face >= 0) {
                        const auto it = snap.byFace.find(fb.face);
                        if (it != snap.byFace.end() && it->second) {
                            fb.highlighted = true;
                        }
                        if (oldExt.count(fb.face) != 0 && snap.exterior) {
                            fb.highlighted = true;
                        }
                        const auto pinIt = snap.pinnedByFace.find(fb.face);
                        if (pinIt != snap.pinnedByFace.end() && pinIt->second) {
                            fb.pinned = true;
                        }
                    }
                    else {
                        const int strand = CirclePadStrand(*necklace, fb);
                        const auto it = snap.byCircle.find(strand);
                        if (strand >= 0 && it != snap.byCircle.end() && it->second) {
                            fb.highlighted = true;
                        }
                        const auto pinIt = snap.pinnedByCircle.find(strand);
                        if (strand >= 0 && pinIt != snap.pinnedByCircle.end() &&
                            pinIt->second) {
                            fb.pinned = true;
                        }
                    }
                }
            }
            bool newExtWasFiniteHighlight = false;
            bool newExtWasOldExt = false;
            for (const int f : newExt) {
                const auto it = snap.byFace.find(f);
                if (it != snap.byFace.end() && it->second) {
                    newExtWasFiniteHighlight = true;
                }
                if (oldExt.count(f) != 0) {
                    newExtWasOldExt = true;
                }
            }
            exteriorHighlighted = newExtWasFiniteHighlight || (newExtWasOldExt && snap.exterior);
        }

        int CircleCount() const
        {
            int n = 0;
            for (const MultiloopComponent& component : components) {
                if (component.kind == MultiloopPart::Circle) {
                    ++n;
                }
            }
            return n;
        }

        int CrossingComponentCount() const
        {
            int n = 0;
            for (const MultiloopComponent& component : components) {
                if (component.kind == MultiloopPart::Crossings) {
                    ++n;
                }
            }
            return n;
        }

        int DrawnStrandCount() const
        {
            if (!loop) {
                return CircleCount();
            }
            int extra = 0;
            for (const MultiloopComponent& component : components) {
                if (component.kind == MultiloopPart::Circle && component.strand < 0) {
                    ++extra;
                }
            }
            return loop->numStrands + extra;
        }

        void SyncSigmaText()
        {
            sigmaBox.text = FormatSigmaNotation(loop ? loop->sigma : Permutation(), components);
            sigmaBox.caret = 0;
            sigmaBox.anchor = 0;
            sigmaBox.scrollX = 0.0f;
        }

        void RebuildPacked()
        {
            layout.reset();
            necklace.reset();
            physics.reset();
            error.clear();
            exteriorHighlighted = false;
            // Packing side by side un-nests every component.
            ClearNesting(components);
            try {
                const Permutation empty;
                const Multiloop dummy(empty);
                const Multiloop& src = loop ? *loop : dummy;
                PackedGridLayout packed = ComputePackedGridLayout(src, components);
                necklace.reset(new Necklace(BuildPackedNecklace(src, packed, components)));
                layout.reset(new GridLayout(std::move(packed.layout)));
                if (!loop) {
                    physics.reset(new Physics(*necklace, physicsParams));
                }
                else {
                    physics.reset(new Physics(*loop, *necklace, physicsParams, PhysicsExteriors()));
                }
            }
            catch (const std::exception& e) {
                layout.reset();
                necklace.reset();
                physics.reset();
                error = e.what();
            }
            if (layout && necklace) {
                RequestViewFit();
            }
        }

        void EnterCircles(std::vector<MultiloopComponent> comps)
        {
            loop.reset();
            components = std::move(comps);
            unknotMode = true;
            dropdownOpen = false;
            dropdownScroll = 0.0f;
            debugFaceBeadIndex = -1;
            ClearR3Undo();
            error.clear();
            const int n = static_cast<int>(components.size());
            if (n <= 1) {
                layout.reset(new GridLayout(UnknotGridLayout()));
                necklace.reset(new Necklace(BuildUnknotNecklace()));
                faceLabels = {"()", "()'"};
                selectedFace = 0;
                exteriorHighlighted = false;
                sigmaBox.text = "()";
                sigmaBox.caret = 0;
                sigmaBox.anchor = 0;
                sigmaBox.scrollX = 0.0f;
                try {
                    physics.reset(new Physics(*necklace, physicsParams));
                }
                catch (const std::exception& e) {
                    physics.reset();
                    error = e.what();
                }
                RequestViewFit();
                return;
            }
            faceLabels.clear();
            selectedFace = -1;
            SyncSigmaText();
            RebuildPacked();
        }

        void EnterCirclesMode(int n)
        {
            std::vector<MultiloopComponent> comps;
            comps.reserve(static_cast<size_t>(std::max(0, n)));
            for (int i = 0; i < n; ++i) {
                MultiloopComponent circle;
                circle.kind = MultiloopPart::Circle;
                comps.push_back(std::move(circle));
            }
            EnterCircles(std::move(comps));
        }

        void SubmitSigma()
        {
            error.clear();
            loop.reset();
            layout.reset();
            necklace.reset();
            physics.reset();
            faceLabels.clear();
            selectedFace = -1;
            dropdownOpen = false;
            dropdownScroll = 0.0f;
            unknotMode = false;
            components.clear();
            ClearR3Undo();
            try {
                const SigmaNotation parsed = ParseSigmaNotation(sigmaBox.text);
                if (parsed.components.empty()) {
                    error = "sigma is empty";
                    return;
                }
                int circles = 0;
                int crossings = 0;
                for (const MultiloopComponent& component : parsed.components) {
                    if (component.kind == MultiloopPart::Circle) {
                        ++circles;
                    }
                    else {
                        ++crossings;
                    }
                }
                if (crossings == 0) {
                    EnterCirclesMode(circles);
                    return;
                }
                const std::string message = ValidateMultiloopSigma(parsed.sigma);
                if (!message.empty()) {
                    error = message;
                    return;
                }
                loop.reset(new Multiloop(parsed.sigma));
                components = parsed.components;
                ApplyDefaultExteriors(*loop, components);
                for (const std::vector<int>& face : loop->faces) {
                    faceLabels.push_back(CycleToString(face));
                }
                SyncSigmaText();
                selectedFace = -1;
                for (const MultiloopComponent& component : components) {
                    if (component.kind == MultiloopPart::Crossings && component.exteriorFace >= 0) {
                        selectedFace = component.exteriorFace;
                        break;
                    }
                }
                RebuildPacked();
            }
            catch (const std::exception& e) {
                error = e.what();
            }
        }

        void EnterUnknotMode()
        {
            EnterCirclesMode(1);
        }

        void SelectFace(int index)
        {
            if (unknotMode) {
                selectedFace = index;
                return; // both faces give the same picture
            }
            if (!loop || index < 0 || index >= static_cast<int>(loop->faces.size())) {
                return;
            }
            const HighlightSnapshot snap = CaptureHighlights();
            const std::unordered_set<int> oldExt = PhysicsExteriors();
            selectedFace = index;
            for (MultiloopComponent& component : components) {
                if (ComponentOwnsFace(component, *loop, index)) {
                    component.exteriorFace = index;
                    break;
                }
            }
            ClearR3Undo();
            RebuildPacked();
            RestoreHighlights(snap, oldExt, PhysicsExteriors());
        }

        void ResetToGrid()
        {
            if (unknotMode) {
                if (physics && physics->AutoTightening()) {
                    physics->StopAutoTighten();
                }
                debugFaceBeadIndex = -1;
                const HighlightSnapshot snap = CaptureHighlights();
                const std::unordered_set<int> oldExt = PhysicsExteriors();
                if (CircleCount() <= 1) {
                    necklace.reset(new Necklace(BuildUnknotNecklace()));
                    physics.reset();
                    try {
                        physics.reset(new Physics(*necklace, physicsParams));
                    }
                    catch (const std::exception& e) {
                        physics.reset();
                        error = e.what();
                    }
                    RestoreHighlights(snap, oldExt, PhysicsExteriors());
                    RequestViewFit();
                    return;
                }
                RebuildPacked();
                RestoreHighlights(snap, oldExt, PhysicsExteriors());
                return;
            }
            if (!loop) {
                return;
            }
            if (physics && physics->AutoTightening()) {
                physics->StopAutoTighten();
            }
            debugFaceBeadIndex = -1;

            std::unordered_map<int, int> labelByFace;
            const HighlightSnapshot snap = CaptureHighlights();
            const std::unordered_set<int> oldExt = PhysicsExteriors();
            if (necklace) {
                for (const FaceBead& fb : necklace->faceBeads) {
                    if (fb.face >= 0 && fb.label > 0) {
                        labelByFace[fb.face] = fb.label;
                    }
                }
            }

            RebuildPacked();

            if (necklace && !labelByFace.empty()) {
                for (FaceBead& fb : necklace->faceBeads) {
                    const auto it = labelByFace.find(fb.face);
                    if (it != labelByFace.end()) {
                        fb.label = it->second;
                    }
                }
            }
            RestoreHighlights(snap, oldExt, PhysicsExteriors());
        }

        // Swap in the rewritten multiloop once the R3 slide has converged.
        // The necklace keeps its bead positions; only topology and labels
        // change, so the picture is undisturbed.
        bool CompleteR3Move()
        {
            if (!physics || !loop || !necklace) {
                return false;
            }
            R3FaceInfo info;
            int meeting[3] = {-1, -1, -1};
            physics->GetR3Result(&info, meeting);
            const std::vector<int> oldExteriors = ExteriorFaceIndices(components);
            std::vector<int> newExteriors;
            std::unique_ptr<Multiloop> newLoop =
                FinalizeR3(*loop, *necklace, info, meeting, oldExteriors, &newExteriors);
            if (!newLoop) {
                physics->CancelR3Move();
                return false;
            }
            physics->EnsureR3ExteriorHoneycombMins(*newLoop);
            physics.reset();
            // Nesting parents follow finite remnants of their faces.
            {
                const std::unordered_set<int> skipExt(newExteriors.begin(), newExteriors.end());
                for (MultiloopComponent& component : components) {
                    if (component.parentFace >= 0) {
                        component.parentFace =
                            RemapFaceByOverlap(*loop, *newLoop, component.parentFace, &skipExt);
                    }
                }
            }
            loop = std::move(newLoop);
            debugFaceBeadIndex = -1;
            {
                size_t k = 0;
                for (MultiloopComponent& component : components) {
                    if (component.kind != MultiloopPart::Crossings) {
                        continue;
                    }
                    if (k < newExteriors.size()) {
                        if (component.exteriorFace == selectedFace) {
                            selectedFace = newExteriors[k];
                        }
                        component.exteriorFace = newExteriors[k];
                    }
                    ++k;
                }
            }
            if (!newExteriors.empty() && selectedFace < 0) {
                selectedFace = newExteriors.front();
            }
            faceLabels.clear();
            for (const std::vector<int>& face : loop->faces) {
                faceLabels.push_back(CycleToString(face));
            }
            SyncSigmaText();
            try {
                physics.reset(new Physics(*loop, *necklace, physicsParams, PhysicsExteriors()));
                physics->SetNesting(components);
                physics->StartPostRepair();
            }
            catch (const std::exception& e) {
                physics.reset();
                error = e.what();
                return false;
            }
            return true;
        }

        // True when some component sits inside `face` (a region with holes).
        bool FaceHasNested(int face) const
        {
            if (face < 0) {
                return false;
            }
            for (const MultiloopComponent& component : components) {
                if (component.parentFace == face) {
                    return true;
                }
            }
            return false;
        }

        // True when some component sits inside the circle of `strand`.
        bool CircleHasNested(int strand) const
        {
            if (strand < 0) {
                return false;
            }
            for (const MultiloopComponent& component : components) {
                if (component.parentCircle == strand) {
                    return true;
                }
            }
            return false;
        }

        enum class MoveHint
        {
            None,
            R0,
            R1,
            R2,
            R3
        };

        int PadStrandOf(const FaceBead& fb) const
        {
            if (!necklace) {
                return -1;
            }
            for (const int bead : fb.boundaryBeads) {
                if (bead < 0 || bead >= static_cast<int>(necklace->beads.size())) {
                    continue;
                }
                const Bead& b = necklace->beads[static_cast<size_t>(bead)];
                if (!b.removed && !b.isVertex) {
                    return b.strand;
                }
            }
            return -1;
        }

        // Legal R0/R1/R2/R3 faces only (Analyze*, not mere degree).
        MoveHint FaceBeadMoveHint(const FaceBead& fb) const
        {
            if (fb.pinned) {
                return MoveHint::None;
            }
            const bool circlePad = unknotMode || fb.face < 0;
            if (circlePad && fb.holeBeads.empty() && !CircleHasNested(PadStrandOf(fb))) {
                return MoveHint::R0;
            }
            if (!loop || unknotMode || fb.face < 0 || FaceHasNested(fb.face)) {
                return MoveHint::None;
            }
            {
                R1FaceInfo info;
                if (AnalyzeR1Face(*loop, nullptr, fb.face, &info)) {
                    return MoveHint::R1;
                }
            }
            {
                R2FaceInfo info;
                if (AnalyzeR2Face(*loop, nullptr, fb.face, &info)) {
                    return MoveHint::R2;
                }
            }
            {
                R3FaceInfo info;
                if (AnalyzeR3Face(*loop, nullptr, fb.face, &info)) {
                    return MoveHint::R3;
                }
            }
            return MoveHint::None;
        }

        // Start R0/R1/R2/R3 from an unpinned region. No-op if the move is
        // invalid (busy physics, nested components, wrong degree, etc.).
        bool TryStartReidemeisterFromFaceBead(int fi)
        {
            if (!physics || !necklace || fi < 0 ||
                fi >= static_cast<int>(necklace->faceBeads.size())) {
                return false;
            }
            const FaceBead& fb = necklace->faceBeads[static_cast<size_t>(fi)];
            const bool idle = !physics->AutoTightening() && !physics->MoveActive();
            if (!idle) {
                return false;
            }
            const MoveHint hint = FaceBeadMoveHint(fb);
            if (hint == MoveHint::None) {
                return false;
            }
            PushR3Undo();
            bool started = false;
            if (hint == MoveHint::R0) {
                started = physics->StartR0Move(fi);
            }
            else if (hint == MoveHint::R1) {
                started = physics->StartR1Move(fb.face);
            }
            else if (hint == MoveHint::R2) {
                started = physics->StartR2Move(fb.face);
            }
            else {
                started = physics->StartR3Move(fb.face);
            }
            if (!started && !r3Undo.empty()) {
                r3Undo.pop_back();
            }
            return started;
        }

        void TryStartReidemeisterAt(Vector2 gridPos)
        {
            if (!physics || !necklace) {
                return;
            }
            const int region = physics->PickRegionAt(gridPos);
            if (region < 0 || region >= static_cast<int>(necklace->faceBeads.size())) {
                return;
            }
            TryStartReidemeisterFromFaceBead(region);
        }

        // Pins every face bead, including circle pads, without an undo snapshot.
        void PinEveryFaceBead()
        {
            if (!necklace) {
                return;
            }
            for (FaceBead& fb : necklace->faceBeads) {
                fb.pinned = true;
            }
        }

        void SetAllPinned(bool value)
        {
            if (!necklace) {
                return;
            }
            bool any = false;
            for (const FaceBead& fb : necklace->faceBeads) {
                if (fb.pinned != value) {
                    any = true;
                    break;
                }
            }
            if (!any) {
                return;
            }
            PushR3Undo();
            for (FaceBead& fb : necklace->faceBeads) {
                fb.pinned = value;
            }
        }

        bool AnyPinned() const
        {
            if (!necklace) {
                return false;
            }
            for (const FaceBead& fb : necklace->faceBeads) {
                if (fb.pinned) {
                    return true;
                }
            }
            return false;
        }

        // An unpinned region whose click would start R0, R1, R2, or R3.
        bool HasPlayableMove() const
        {
            if (!necklace) {
                return false;
            }
            for (const FaceBead& fb : necklace->faceBeads) {
                if (FaceBeadMoveHint(fb) != MoveHint::None) {
                    return true;
                }
            }
            return false;
        }

        std::string RegionLogName(const FaceBead& fb) const
        {
            if (fb.face >= 0) {
                return "region " + std::to_string(fb.label);
            }
            const int strand = necklace ? CirclePadStrand(*necklace, fb) : -1;
            if (strand >= 0) {
                return "circle strand " + std::to_string(strand);
            }
            return "circle " + std::to_string(fb.label);
        }

        std::string DescribeStartingPins() const
        {
            if (!necklace) {
                return "none";
            }
            std::vector<std::pair<int, std::string>> items;
            for (const FaceBead& fb : necklace->faceBeads) {
                if (!fb.pinned) {
                    continue;
                }
                if (fb.face >= 0) {
                    items.push_back({fb.label, "region " + std::to_string(fb.label)});
                }
                else {
                    const int strand = CirclePadStrand(*necklace, fb);
                    const int key = 100000 + (strand >= 0 ? strand : fb.label);
                    items.push_back({key, RegionLogName(fb)});
                }
            }
            if (items.empty()) {
                return "none";
            }
            std::sort(items.begin(), items.end());
            std::string out;
            for (size_t i = 0; i < items.size(); ++i) {
                if (i > 0) {
                    out += ", ";
                }
                out += items[i].second;
            }
            return out;
        }

        std::string InfiniteRegionText() const
        {
            if (selectedFace >= 0 && selectedFace < static_cast<int>(faceLabels.size())) {
                return faceLabels[static_cast<size_t>(selectedFace)];
            }
            if (unknotMode) {
                return "()";
            }
            return "(none)";
        }

        static UnpinningGame::MoveKind HintToKind(MoveHint hint)
        {
            if (hint == MoveHint::R1) {
                return UnpinningGame::MoveKind::R1;
            }
            if (hint == MoveHint::R2) {
                return UnpinningGame::MoveKind::R2;
            }
            if (hint == MoveHint::R3) {
                return UnpinningGame::MoveKind::R3;
            }
            return UnpinningGame::MoveKind::R0;
        }

        void BeginUnpinningMatch()
        {
            if (!layout || !necklace) {
                return;
            }
            RememberGameStart();
            ClearR3Undo();
            dropdownOpen = false;
            tuningUI.open = false;
            game.Begin(sigmaBox.text, InfiniteRegionText(), DescribeStartingPins());
            pendingGameFit = true;
        }

        void ReplayUnpinningMatch()
        {
            RestoreGameStart();
            PinEveryFaceBead();
            ClearR3Undo();
            dropdownOpen = false;
            tuningUI.open = false;
            game.Begin(sigmaBox.text, InfiniteRegionText(), DescribeStartingPins());
            RequestViewFit();
        }

        void ConfirmExitToSandbox()
        {
            game.MarkEarlyExit();
            ClearR3Undo();
            dropdownOpen = false;
            RequestViewFit();
        }

        void BackToSandboxFromGame()
        {
            RestoreGameStart();
            ClearR3Undo();
            game.Deactivate();
            dropdownOpen = false;
            RequestViewFit();
        }

        void HandleGameClick(Vector2 gridMouse, float pickDist)
        {
            if (!game.Active() || game.ExitPrompt() || game.Over() || !physics || !necklace) {
                return;
            }
            if (physics->MoveActive() || physics->AutoTightening() || physics->Dragging()) {
                return;
            }
            const int body = physics->PickBody(gridMouse, pickDist);
            const int beadCount = static_cast<int>(necklace->beads.size());
            if (body >= beadCount) {
                const int fi = body - beadCount;
                if (fi >= 0 && fi < static_cast<int>(necklace->faceBeads.size())) {
                    FaceBead& fb = necklace->faceBeads[static_cast<size_t>(fi)];
                    if (fb.pinned && game.AcceptsPinRemoval()) {
                        const std::string name = RegionLogName(fb);
                        fb.pinned = false;
                        game.OnPinRemoved(name);
                        return;
                    }
                }
            }
            if (!game.AcceptsMove()) {
                return;
            }
            const int region = physics->PickRegionAt(gridMouse);
            if (region < 0 || region >= static_cast<int>(necklace->faceBeads.size())) {
                return;
            }
            const FaceBead& fb = necklace->faceBeads[static_cast<size_t>(region)];
            const MoveHint hint = FaceBeadMoveHint(fb);
            if (hint == MoveHint::None) {
                return;
            }
            const std::string name = RegionLogName(fb);
            if (!TryStartReidemeisterFromFaceBead(region)) {
                return;
            }
            game.OnMoveStarted(HintToKind(hint), name);
        }

        void ResolveTangleReset()
        {
            game.OnResetPressed();
            ResetToGrid();
            bool started = false;
            if (physics && !physics->MoveActive() && !physics->AutoTightening()) {
                physics->StartAutoTighten();
                started = physics->AutoTightening();
            }
            if (!started && physics) {
                physics->RefreshTangle();
                game.OnSettled(physics->Tangled(), AnyPinned(), GetTime());
            }
        }

        // Face of `newLoop` sharing the most half-edges with `oldFace` of
        // `oldLoop`. Faces in `skip` (new exteriors) are ignored so a parent
        // that split onto the unbounded region keeps a finite remnant.
        static int RemapFaceByOverlap(const Multiloop& oldLoop, const Multiloop& newLoop,
                                      int oldFace, const std::unordered_set<int>* skip = nullptr)
        {
            if (oldFace < 0 || oldFace >= static_cast<int>(oldLoop.faces.size())) {
                return -1;
            }
            const std::vector<int>& cycle = oldLoop.faces[static_cast<size_t>(oldFace)];
            int best = -1;
            size_t bestOverlap = 0;
            for (size_t f = 0; f < newLoop.faces.size(); ++f) {
                if (skip != nullptr && skip->count(static_cast<int>(f)) != 0) {
                    continue;
                }
                size_t overlap = 0;
                for (const int h : newLoop.faces[f]) {
                    if (std::find(cycle.begin(), cycle.end(), h) != cycle.end()) {
                        ++overlap;
                    }
                }
                if (overlap > bestOverlap) {
                    bestOverlap = overlap;
                    best = static_cast<int>(f);
                }
            }
            return best;
        }

        // Shared tail of CompleteR1Move / CompleteR2Move: the necklace has been
        // rewired and `result` carries the new component list.
        void CompleteShorteningMove(std::unique_ptr<Multiloop> newLoop,
                                    const MoveComponentResult& result)
        {
            physics->EnsureExteriorHoneycombMins(*newLoop, result.exteriorArcs);
            physics.reset();
            int newSelected = -1;
            if (selectedFace >= 0 && selectedFace < static_cast<int>(result.faceMap.size())) {
                newSelected = result.faceMap[static_cast<size_t>(selectedFace)];
            }
            loop = std::move(newLoop);
            components = result.components;
            debugFaceBeadIndex = -1;
            if (newSelected < 0) {
                for (const MultiloopComponent& component : components) {
                    if (component.kind == MultiloopPart::Crossings && component.exteriorFace >= 0) {
                        newSelected = component.exteriorFace;
                        break;
                    }
                }
            }
            if (newSelected >= static_cast<int>(loop->faces.size())) {
                newSelected = -1;
            }
            selectedFace = newSelected;
            faceLabels.clear();
            for (const std::vector<int>& face : loop->faces) {
                faceLabels.push_back(CycleToString(face));
            }
            SyncSigmaText();
            try {
                // With no crossings left this is the all-circles multiloop; the
                // necklace stays where it is (Reset packs side by side).
                physics.reset(new Physics(*loop, *necklace, physicsParams, PhysicsExteriors()));
                physics->SetNesting(components);
                if (result.mergeKeepLabel >= 0 && result.mergeDropLabel >= 0) {
                    physics->SetFaceMergePair(result.mergeKeepLabel, result.mergeDropLabel);
                }
                physics->StartPostRepair();
            }
            catch (const std::exception& e) {
                physics.reset();
                error = e.what();
            }
        }

        bool CompleteR1Move()
        {
            if (!physics || !loop || !necklace) {
                return false;
            }
            R1FaceInfo info;
            physics->GetR1Result(&info);
            MoveComponentResult result;
            std::unique_ptr<Multiloop> newLoop =
                FinalizeR1(*loop, *necklace, info, components, &result);
            if (!newLoop) {
                physics->CancelMove();
                return false;
            }
            CompleteShorteningMove(std::move(newLoop), result);
            return static_cast<bool>(physics);
        }

        bool CompleteR2Move()
        {
            if (!physics || !loop || !necklace) {
                return false;
            }
            R2FaceInfo info;
            physics->GetR2Result(&info);
            MoveComponentResult result;
            std::unique_ptr<Multiloop> newLoop =
                FinalizeR2(*loop, *necklace, info, components, &result);
            if (!newLoop) {
                physics->CancelMove();
                return false;
            }
            CompleteShorteningMove(std::move(newLoop), result);
            return static_cast<bool>(physics);
        }

        // Empty canvas after the last circle pops. Leaves the undo stack intact.
        void ClearDiagramKeepUndo()
        {
            loop.reset();
            components.clear();
            unknotMode = false;
            necklace.reset();
            physics.reset();
            layout.reset();
            faceLabels.clear();
            selectedFace = -1;
            sigmaBox.text.clear();
            sigmaBox.caret = 0;
            sigmaBox.anchor = 0;
            sigmaBox.scrollX = 0.0f;
            debugFaceBeadIndex = -1;
            error.clear();
            dropdownOpen = false;
            wasAutoTightening = false;
            wasR3Active = false;
        }

        void RefreshUnknotLabels()
        {
            if (CircleCount() <= 1) {
                faceLabels = {"()", "()'"};
                selectedFace = 0;
            }
            else {
                faceLabels.clear();
                selectedFace = -1;
            }
        }

        enum class MoveFinish
        {
            Failed,
            Done,
            Cleared
        };

        MoveFinish CompleteR0Move()
        {
            if (!physics || !necklace) {
                return MoveFinish::Failed;
            }
            R0MoveInfo info;
            physics->GetR0Result(&info);
            if (!FinalizeR0(*necklace, info, &components)) {
                physics->CancelMove();
                return MoveFinish::Failed;
            }
            physics.reset();
            debugFaceBeadIndex = -1;
            if (components.empty()) {
                ClearDiagramKeepUndo();
                return MoveFinish::Cleared;
            }
            SyncSigmaText();
            bool allCircles = true;
            for (const MultiloopComponent& component : components) {
                if (component.kind != MultiloopPart::Circle) {
                    allCircles = false;
                    break;
                }
            }
            try {
                if (unknotMode || !loop || allCircles || loop->numEdges == 0) {
                    physics.reset(new Physics(*necklace, physicsParams));
                    physics->SetNesting(components);
                    unknotMode = true;
                    loop.reset();
                    RefreshUnknotLabels();
                }
                else {
                    physics.reset(new Physics(*loop, *necklace, physicsParams, PhysicsExteriors()));
                    physics->SetNesting(components);
                    physics->StartPostRepair();
                }
            }
            catch (const std::exception& e) {
                physics.reset();
                error = e.what();
                return MoveFinish::Failed;
            }
            return MoveFinish::Done;
        }

        void SelectFaceContaining(int halfEdge)
        {
            for (size_t f = 0; f < loop->faces.size(); ++f) {
                for (const int h : loop->faces[f]) {
                    if (h == halfEdge) {
                        SelectFace(static_cast<int>(f));
                        return;
                    }
                }
            }
        }

        void GenerateRandom()
        {
            error.clear();
            ClearR3Undo();
            int count = 0;
            try {
                count = std::stoi(countBox.text);
            }
            catch (const std::exception&) {
                error = "Enter a valid number of faces";
                return;
            }
            count = std::max(1, count);
            countBox.text = std::to_string(count);

            if (count == 1) {
                // A simple closed curve: no crossings.
                EnterUnknotMode();
                return;
            }
            if (count == 2) {
                if (forceSimple) {
                    // Two overlapping simple circles (the vesica): two
                    // crossings, no self-intersections, outer face on edge 1.
                    sigmaBox.text = "(1,2,3,4)(-4,-3,-2,-1)";
                    SubmitSigma();
                    if (loop) {
                        SelectFaceContaining(1);
                    }
                    return;
                }
                // A figure-eight curl: one crossing, two monogon faces, with
                // the bigon (-1,-2) as the infinite region.
                sigmaBox.text = "(1,-1,2,-2)";
                SubmitSigma();
                if (loop) {
                    SelectFaceContaining(-1);
                }
                return;
            }

            try {
                const RandomMultiloop random =
                    GenerateRandomMultiloop(count, density, forceOneStrand, forceSimple, rng());
                sigmaBox.text = random.sigma.ToString();
                SubmitSigma();
                if (loop) {
                    // Default infinite region: the unbounded region of the
                    // medial graph.
                    SelectFaceContaining(random.exteriorHalfEdge);
                }
            }
            catch (const std::exception& e) {
                error = e.what();
            }
        }

        // ---- free-hand drawing ----

        // A single tolerance governs both placement rules.
        static constexpr float kTol = 5.0f;
        static constexpr float kCloseTol = 18.0f;        // snap radius for the start vertex
        static constexpr double kDoubleClickTime = 0.1; // max seconds between double-click presses

        void StartDrawing()
        {
            savedSigmaText = sigmaBox.text;
            savedFace = selectedFace;
            savedUnknot = unknotMode;
            drawing = true;
            drawCurves.clear();
            drawCurrent.clear();
            drawCrossings.clear();
            drawWarning.clear();
            loop.reset();
            layout.reset();
            necklace.reset();
            physics.reset();
            faceLabels.clear();
            selectedFace = -1;
            unknotMode = false;
            components.clear();
            dropdownOpen = false;
            error.clear();
        }

        void CancelDrawing()
        {
            drawing = false;
            drawCurves.clear();
            drawCurrent.clear();
            drawCrossings.clear();
            drawWarning.clear();
            if (!savedSigmaText.empty()) {
                sigmaBox.text = savedSigmaText;
                SubmitSigma();
                if (loop && savedFace >= 0 && savedFace < static_cast<int>(faceLabels.size())) {
                    SelectFace(savedFace);
                }
            }
            else if (savedUnknot) {
                EnterUnknotMode();
            }
        }

        void FinishDrawing()
        {
            // Only closed curves contribute; an unfinished open polyline is
            // ignored. Failures keep drawing mode active so the user can fix
            // the diagram and try again.
            if (drawCurves.empty()) {
                Warn("Draw and close at least one loop first.", GetMousePosition());
                return;
            }
            const DrawnMultiloop drawn = BuildMultiloopFromCurves(drawCurves);
            if (!drawn.error.empty()) {
                Warn(drawn.error, GetMousePosition());
                return;
            }

            drawing = false;
            drawWarning.clear();
            drawCurrent.clear();
            drawCrossings.clear();
            try {
                if (drawn.isUnknot) {
                    EnterUnknotMode();
                    return;
                }
                int crossings = 0;
                for (const MultiloopComponent& component : drawn.components) {
                    if (component.kind == MultiloopPart::Crossings) {
                        ++crossings;
                    }
                }
                if (crossings == 0) {
                    EnterCircles(drawn.components);
                    return;
                }
                loop.reset(new Multiloop(drawn.sigma));
                ApplyDrawnStrandColors(*loop, drawn);
                components = drawn.components;
                ApplyDefaultExteriors(*loop, components);
                unknotMode = false;
                ClearR3Undo();
                faceLabels.clear();
                for (const std::vector<int>& face : loop->faces) {
                    faceLabels.push_back(CycleToString(face));
                }
                SyncSigmaText();
                if (drawn.exteriorHalfEdge != 0) {
                    selectedFace = FaceIndexContainingHalfEdge(*loop, drawn.exteriorHalfEdge);
                }
                if (selectedFace < 0) {
                    for (const MultiloopComponent& component : components) {
                        if (component.kind == MultiloopPart::Crossings &&
                            component.exteriorFace >= 0) {
                            selectedFace = component.exteriorFace;
                            break;
                        }
                    }
                }
                RebuildPacked();
                if (!error.empty()) {
                    drawing = true;
                    Warn(error, GetMousePosition());
                }
            }
            catch (const std::exception& e) {
                drawing = true;
                error = e.what();
                Warn(error, GetMousePosition());
            }
        }

        // Removes the most recently drawn point. When the current polyline is
        // empty, the last completed loop is reopened so its points can be
        // edited too.
        void DeleteLastPoint()
        {
            if (drawCurrent.empty() && !drawCurves.empty()) {
                drawCurrent = drawCurves.back();
                drawCurves.pop_back();
            }
            if (!drawCurrent.empty()) {
                drawCurrent.pop_back();
            }
            drawWarning.clear();
            RecomputeDrawCrossings();
        }

        // All committed segments (completed curves plus the in-progress
        // polyline) as endpoint pairs.
        std::vector<std::pair<Vector2, Vector2>> CommittedSegments() const
        {
            std::vector<std::pair<Vector2, Vector2>> segs;
            for (const std::vector<Vector2>& curve : drawCurves) {
                for (size_t i = 0; i < curve.size(); ++i) {
                    segs.push_back({curve[i], curve[(i + 1) % curve.size()]});
                }
            }
            for (size_t i = 0; i + 1 < drawCurrent.size(); ++i) {
                segs.push_back({drawCurrent[i], drawCurrent[i + 1]});
            }
            return segs;
        }

        void RecomputeDrawCrossings()
        {
            drawCrossings.clear();
            std::vector<std::pair<Vector2, Vector2>> segs;
            std::vector<int> owner; // curve index, -1 for the open polyline
            for (size_t c = 0; c < drawCurves.size(); ++c) {
                const std::vector<Vector2>& curve = drawCurves[c];
                for (size_t i = 0; i < curve.size(); ++i) {
                    segs.push_back({curve[i], curve[(i + 1) % curve.size()]});
                    owner.push_back(static_cast<int>(c));
                }
            }
            for (size_t i = 0; i + 1 < drawCurrent.size(); ++i) {
                segs.push_back({drawCurrent[i], drawCurrent[i + 1]});
                owner.push_back(-1);
            }
            for (size_t i = 0; i < segs.size(); ++i) {
                for (size_t j = i + 1; j < segs.size(); ++j) {
                    // Skip segments that meet at a shared polyline vertex.
                    if (owner[i] == -1 && owner[j] == -1 &&
                        static_cast<int>(j) - static_cast<int>(i) == 1) {
                        continue;
                    }
                    if (owner[i] >= 0 && owner[i] == owner[j]) {
                        const int m = static_cast<int>(drawCurves[owner[i]].size());
                        const int diff = std::abs(static_cast<int>(i) - static_cast<int>(j));
                        if (diff == 1 || diff == m - 1) {
                            continue;
                        }
                    }
                    Vector2 pt;
                    if (SegmentsProperlyCross(segs[i].first, segs[i].second,
                                              segs[j].first, segs[j].second, &pt)) {
                        drawCrossings.push_back(pt);
                    }
                }
            }
        }

        // Every corner of the committed geometry (polyline vertices). A new
        // segment must keep clear of these as well as of crossings.
        std::vector<Vector2> CommittedCorners() const
        {
            std::vector<Vector2> corners;
            for (const std::vector<Vector2>& curve : drawCurves) {
                for (const Vector2& p : curve) {
                    corners.push_back(p);
                }
            }
            for (const Vector2& p : drawCurrent) {
                corners.push_back(p);
            }
            return corners;
        }

        // Rule (b): a new segment a->b may not pass within tolerance of any
        // existing corner or crossing (its own endpoints excepted).
        bool SegmentPlaceable(Vector2 a, Vector2 b, std::string* reason) const
        {
            const auto grazes = [&](Vector2 v) {
                if (std::hypot(v.x - a.x, v.y - a.y) < 1e-3f ||
                    std::hypot(v.x - b.x, v.y - b.y) < 1e-3f) {
                    return false; // shared endpoint, not a violation
                }
                return DistPointToSegment(v, a, b) < kTol;
            };
            for (const Vector2& corner : CommittedCorners()) {
                if (grazes(corner)) {
                    *reason = "Segment passes too close to an existing corner.";
                    return false;
                }
            }
            for (const Vector2& c : drawCrossings) {
                if (grazes(c)) {
                    *reason = "Segment passes too close to an existing crossing.";
                    return false;
                }
            }
            return true;
        }

        // Rule (a): a new point may not be placed within tolerance of an
        // existing segment.
        bool PointPlaceable(Vector2 p, std::string* reason) const
        {
            for (const std::pair<Vector2, Vector2>& s : CommittedSegments()) {
                if (DistPointToSegment(p, s.first, s.second) < kTol) {
                    *reason = "Too close to an existing strand.";
                    return false;
                }
            }
            return true;
        }

        void Warn(const std::string& message, Vector2 at)
        {
            drawWarning = message;
            drawWarningPos = at;
            drawWarningUntil = GetTime() + 2.5;
        }

        // Closes the in-progress polyline back to its start vertex. Returns
        // false (with a warning) when the closing segment is illegal.
        bool CloseCurrentLoop(Vector2 warnAt)
        {
            if (drawCurrent.size() < 3) {
                Warn("Need at least 3 points to close a loop.", warnAt);
                return false;
            }
            std::string reason;
            if (!SegmentPlaceable(drawCurrent.back(), drawCurrent.front(), &reason)) {
                Warn(reason, warnAt);
                return false;
            }
            drawCurves.push_back(drawCurrent);
            drawCurrent.clear();
            RecomputeDrawCrossings();
            return true;
        }

        void HandleDrawClick(Vector2 p)
        {
            if (drawCurrent.empty()) {
                std::string reason;
                if (!PointPlaceable(p, &reason)) {
                    Warn(reason, p);
                    return;
                }
                drawCurrent.push_back(p);
                return;
            }
            const Vector2 last = drawCurrent.back();
            // Close by clicking near the start vertex.
            if (drawCurrent.size() >= 3 &&
                std::hypot(p.x - drawCurrent.front().x, p.y - drawCurrent.front().y) <= kCloseTol) {
                CloseCurrentLoop(p);
                return;
            }
            std::string reason;
            if (!PointPlaceable(p, &reason)) {
                Warn(reason, p);
                return;
            }
            if (!SegmentPlaceable(last, p, &reason)) {
                Warn(reason, p);
                return;
            }
            drawCurrent.push_back(p);
            RecomputeDrawCrossings();
        }

        // A double-click closes the current loop. The first click of the pair
        // already added a point at (about) the cursor, so drop it first.
        void HandleDoubleClick(Vector2 p)
        {
            if (!drawCurrent.empty() &&
                std::hypot(drawCurrent.back().x - p.x, drawCurrent.back().y - p.y) < kCloseTol * 1.5f) {
                drawCurrent.pop_back();
                RecomputeDrawCrossings();
            }
            CloseCurrentLoop(p);
        }
    };

    // ----------------------------------------------------------------
    // Drawing the multiloop.
    // ----------------------------------------------------------------

    // Half-edge label near a crossing, kept within a box of side cell/3
    // centered on the crossing.
    void DrawHalfEdgeLabel(Vector2 vertex, Vector2 towards, int halfEdge, float cell)
    {
        Vector2 dir = {towards.x - vertex.x, towards.y - vertex.y};
        const float len = std::sqrt(dir.x * dir.x + dir.y * dir.y);
        if (len > 0.0f) {
            dir.x /= len;
            dir.y /= len;
        }
        const Vector2 perp = {dir.y, -dir.x};

        const std::string text = std::to_string(halfEdge);
        const float boxHalf = cell / 4.5f;
        float fontSize = Clampf(cell * 0.15f, 7.0f, 16.0f);
        float textWidth = TextWidth(text, fontSize);
        if (textWidth + 4.0f > 2.0f * boxHalf) {
            fontSize = std::max(6.0f, fontSize * (2.0f * boxHalf - 4.0f) / (textWidth + 4.0f));
            textWidth = TextWidth(text, fontSize);
        }

        Vector2 pos = {vertex.x + dir.x * cell * 0.17f + perp.x * cell * 0.11f,
                       vertex.y + dir.y * cell * 0.17f + perp.y * cell * 0.11f};
        const float halfW = textWidth / 2.0f + 2.0f;
        const float halfH = fontSize / 2.0f + 1.0f;
        if (boxHalf > halfW) {
            pos.x = Clampf(pos.x, vertex.x - boxHalf + halfW, vertex.x + boxHalf - halfW);
        }
        else {
            pos.x = vertex.x + dir.x * (boxHalf - halfW);
        }
        if (boxHalf > halfH) {
            pos.y = Clampf(pos.y, vertex.y - boxHalf + halfH, vertex.y + boxHalf - halfH);
        }

        DrawRectangle(static_cast<int>(pos.x - halfW), static_cast<int>(pos.y - halfH),
                      static_cast<int>(2.0f * halfW), static_cast<int>(2.0f * halfH),
                      Color{255, 255, 255, 235});
        DrawTextS(text, pos.x - textWidth / 2.0f, pos.y - fontSize / 2.0f, fontSize, BLACK);
    }

    // Mapping between grid coordinates (y up) and screen pixels, shared by
    // rendering and bead picking so that they stay in sync. Rotation is
    // counterclockwise about the look-at; screen y is down.
    struct View
    {
        float cell = 1.0f;
        float lookX = 0.0f;
        float lookY = 0.0f;
        float angle = 0.0f;
        float cx = 0.0f;
        float cy = 0.0f;
        float cosA = 1.0f;
        float sinA = 0.0f;

        Vector2 ToScreen(const BeadPos& g) const
        {
            const float dx = g.x - lookX;
            const float dy = g.y - lookY;
            const float rx = dx * cosA - dy * sinA;
            const float ry = dx * sinA + dy * cosA;
            return Vector2{cx + rx * cell, cy - ry * cell};
        }
        Vector2 ToGrid(Vector2 s) const
        {
            const float rx = (s.x - cx) / cell;
            const float ry = -(s.y - cy) / cell;
            const float dx = rx * cosA + ry * sinA;
            const float dy = -rx * sinA + ry * cosA;
            return Vector2{lookX + dx, lookY + dy};
        }
    };

    constexpr float kPiF = 3.14159265358979323846f;

    float EaseInOutCubic(float t)
    {
        t = Clampf(t, 0.0f, 1.0f);
        if (t < 0.5f) {
            return 4.0f * t * t * t;
        }
        const float u = -2.0f * t + 2.0f;
        return 1.0f - (u * u * u) * 0.5f;
    }

    // AABB size is π-periodic; pick the lift of `candidate` nearest `current`.
    float NearestPiPeriodic(float current, float candidate)
    {
        float d = candidate - current;
        d -= kPiF * std::round(d / kPiF);
        return current + d;
    }

    void LayoutCamera(const GridLayout& layout, Rectangle canvas, float* cell, float* lookX,
                      float* lookY)
    {
        const float margin = 45.0f;
        const float cellW = (canvas.width - 2 * margin) / std::max(1, layout.width);
        const float cellH = (canvas.height - 2 * margin) / std::max(1, layout.height);
        *cell = Clampf(std::min(cellW, cellH), 5.0f, 110.0f);
        *lookX = 0.5f * static_cast<float>(layout.width);
        *lookY = 0.5f * static_cast<float>(layout.height);
    }

    bool BeadCamera(const GridLayout& layout, const Necklace& necklace, Rectangle canvas,
                    float currentAngle, bool searchRotation, float* cell, float* lookX,
                    float* lookY, float* angle)
    {
        std::vector<Vector2> pts;
        pts.reserve(necklace.beads.size() + necklace.faceBeads.size());
        for (const Bead& bead : necklace.beads) {
            if (bead.removed) {
                continue;
            }
            pts.push_back(Vector2{bead.pos.x, bead.pos.y});
        }
        for (const FaceBead& fb : necklace.faceBeads) {
            pts.push_back(Vector2{fb.pos.x, fb.pos.y});
        }
        if (pts.empty()) {
            LayoutCamera(layout, canvas, cell, lookX, lookY);
            *angle = searchRotation ? 0.0f : currentAngle;
            return false;
        }

        const float pad = std::max(4.0f * necklace.radius, 0.75f);
        const float margin = 45.0f;
        const float availW = std::max(canvas.width - 2.0f * margin, 1.0f);
        const float availH = std::max(canvas.height - 2.0f * margin, 1.0f);

        const auto evaluate = [&](float theta, float* outLookX, float* outLookY) {
            const float c = std::cos(theta);
            const float s = std::sin(theta);
            float minX = pts[0].x * c - pts[0].y * s;
            float maxX = minX;
            float minY = pts[0].x * s + pts[0].y * c;
            float maxY = minY;
            for (size_t i = 1; i < pts.size(); ++i) {
                const float rx = pts[i].x * c - pts[i].y * s;
                const float ry = pts[i].x * s + pts[i].y * c;
                minX = std::min(minX, rx);
                maxX = std::max(maxX, rx);
                minY = std::min(minY, ry);
                maxY = std::max(maxY, ry);
            }
            minX -= pad;
            maxX += pad;
            minY -= pad;
            maxY += pad;
            const float w = std::max(maxX - minX, 0.25f);
            const float h = std::max(maxY - minY, 0.25f);
            if (outLookX != nullptr && outLookY != nullptr) {
                const float midRx = 0.5f * (minX + maxX);
                const float midRy = 0.5f * (minY + maxY);
                *outLookX = midRx * c + midRy * s;
                *outLookY = -midRx * s + midRy * c;
            }
            return std::min(availW / w, availH / h);
        };

        if (!searchRotation) {
            const float scale = evaluate(currentAngle, lookX, lookY);
            *cell = Clampf(scale, 3.0f, 250.0f);
            *angle = currentAngle;
            return true;
        }

        constexpr int kCoarse = 180;
        float bestTheta = 0.0f;
        float bestScale = -1.0f;
        for (int i = 0; i < kCoarse; ++i) {
            const float theta = kPiF * static_cast<float>(i) / static_cast<float>(kCoarse);
            const float scale = evaluate(theta, nullptr, nullptr);
            if (scale > bestScale) {
                bestScale = scale;
                bestTheta = theta;
            }
        }
        const float window = kPiF / static_cast<float>(kCoarse);
        constexpr int kFine = 20;
        for (int i = 0; i <= kFine; ++i) {
            const float t = -1.0f + 2.0f * static_cast<float>(i) / static_cast<float>(kFine);
            const float theta = bestTheta + t * window;
            const float scale = evaluate(theta, nullptr, nullptr);
            if (scale > bestScale) {
                bestScale = scale;
                bestTheta = theta;
            }
        }

        evaluate(bestTheta, lookX, lookY);
        *cell = Clampf(bestScale, 3.0f, 250.0f);
        *angle = NearestPiPeriodic(currentAngle, bestTheta);
        return true;
    }

    void BeginViewFit(App& app, Rectangle canvas, bool searchRotation = true)
    {
        if (!app.layout) {
            return;
        }
        float cell = app.viewCell;
        float lookX = app.viewLookX;
        float lookY = app.viewLookY;
        float angle = searchRotation ? 0.0f : app.viewAngle;
        if (app.necklace) {
            BeadCamera(*app.layout, *app.necklace, canvas, app.viewAngle, searchRotation, &cell,
                       &lookX, &lookY, &angle);
        }
        else {
            LayoutCamera(*app.layout, canvas, &cell, &lookX, &lookY);
        }
        const float dCell = std::fabs(cell - app.viewCell);
        const float dLook = std::hypot(lookX - app.viewLookX, lookY - app.viewLookY);
        const float dAngle = std::fabs(angle - app.viewAngle);
        if (dCell < 0.15f && dLook < 0.02f && dAngle < 0.002f) {
            app.viewCell = cell;
            app.viewLookX = lookX;
            app.viewLookY = lookY;
            app.viewAngle = angle;
            app.viewAnimating = false;
            return;
        }
        app.viewCell0 = app.viewCell;
        app.viewLookX0 = app.viewLookX;
        app.viewLookY0 = app.viewLookY;
        app.viewAngle0 = app.viewAngle;
        app.viewCell1 = cell;
        app.viewLookX1 = lookX;
        app.viewLookY1 = lookY;
        app.viewAngle1 = angle;
        app.viewAnimT = 0.0f;
        app.viewAnimating = true;
    }

    View MakeView(App& app, const GridLayout& layout, Rectangle canvas)
    {
        const void* key = &layout;
        if (app.viewFitPending) {
            BeginViewFit(app, canvas);
            app.viewFitPending = false;
            app.viewLayoutKey = key;
            app.viewReady = true;
            app.viewGridHeight = layout.height;
        }
        else if (!app.viewReady || app.viewLayoutKey != key) {
            LayoutCamera(layout, canvas, &app.viewCell, &app.viewLookX, &app.viewLookY);
            app.viewAngle = 0.0f;
            app.viewGridHeight = layout.height;
            app.viewLayoutKey = key;
            app.viewReady = true;
            app.viewAnimating = false;
        }
        app.viewGridHeight = layout.height;

        if (app.viewAnimating) {
            const float dt = std::min(GetFrameTime(), 0.05f);
            app.viewAnimT += dt / 0.50f;
            const float u = EaseInOutCubic(app.viewAnimT);
            app.viewCell = app.viewCell0 + (app.viewCell1 - app.viewCell0) * u;
            app.viewLookX = app.viewLookX0 + (app.viewLookX1 - app.viewLookX0) * u;
            app.viewLookY = app.viewLookY0 + (app.viewLookY1 - app.viewLookY0) * u;
            app.viewAngle = app.viewAngle0 + (app.viewAngle1 - app.viewAngle0) * u;
            if (app.viewAnimT >= 1.0f) {
                app.viewCell = app.viewCell1;
                app.viewLookX = app.viewLookX1;
                app.viewLookY = app.viewLookY1;
                app.viewAngle = app.viewAngle1;
                app.viewAnimating = false;
            }
        }

        View view;
        view.cell = app.viewCell;
        view.lookX = app.viewLookX;
        view.lookY = app.viewLookY;
        view.angle = app.viewAngle;
        view.cx = canvas.x + canvas.width * 0.5f;
        view.cy = canvas.y + canvas.height * 0.5f;
        view.cosA = std::cos(app.viewAngle);
        view.sinA = std::sin(app.viewAngle);
        return view;
    }

    Color LerpColor(Color a, Color b, float t)
    {
        t = Clampf(t, 0.0f, 1.0f);
        Color out;
        out.r = static_cast<unsigned char>(a.r + (static_cast<int>(b.r) - a.r) * t);
        out.g = static_cast<unsigned char>(a.g + (static_cast<int>(b.g) - a.g) * t);
        out.b = static_cast<unsigned char>(a.b + (static_cast<int>(b.b) - a.b) * t);
        out.a = static_cast<unsigned char>(a.a + (static_cast<int>(b.a) - a.a) * t);
        return out;
    }

    void FillPolygonScreen(const std::vector<Vector2>& pts, Color inner, Color outer,
                           Vector2 center)
    {
        if (pts.size() < 3) {
            return;
        }
        float minY = pts[0].y;
        float maxY = pts[0].y;
        float rMaxSq = 0.0f;
        const bool gradient = inner.r != outer.r || inner.g != outer.g || inner.b != outer.b ||
                              inner.a != outer.a;
        for (const Vector2& p : pts) {
            minY = std::min(minY, p.y);
            maxY = std::max(maxY, p.y);
            if (gradient) {
                const float dx = p.x - center.x;
                const float dy = p.y - center.y;
                rMaxSq = std::max(rMaxSq, dx * dx + dy * dy);
            }
        }
        const int y0 = static_cast<int>(std::floor(minY));
        const int y1 = static_cast<int>(std::ceil(maxY));
        if (y1 <= y0) {
            return;
        }
        const float rMax = std::sqrt(rMaxSq);
        const bool useGrad = gradient && rMax >= 1.0f;
        const auto colorAt = [&](float x, float y) {
            if (!useGrad) {
                return inner;
            }
            const float dx = x - center.x;
            const float dy = y - center.y;
            const float t = 1.0f - std::sqrt(dx * dx + dy * dy) / rMax;
            return LerpColor(outer, inner, t);
        };
        const auto emitVert = [&](float x, float y, Color c) {
            rlColor4ub(c.r, c.g, c.b, c.a);
            rlVertex2f(x, y);
        };
        const auto emitQuad = [&](float left, float right, float top, float bot) {
            const Color cLT = colorAt(left, top);
            const Color cRT = colorAt(right, top);
            const Color cRB = colorAt(right, bot);
            const Color cLB = colorAt(left, bot);
            emitVert(left, top, cLT);
            emitVert(right, top, cRT);
            emitVert(right, bot, cRB);
            emitVert(left, top, cLT);
            emitVert(right, bot, cRB);
            emitVert(left, bot, cLB);
        };
        // Even-odd scanline. Optional radial tint: split spans at center.x so
        // the middle of the polygon interpolates to the inner (lighter) color.
        rlDisableBackfaceCulling();
        rlBegin(RL_TRIANGLES);
        if (!useGrad) {
            rlColor4ub(inner.r, inner.g, inner.b, inner.a);
        }
        std::vector<float> xs;
        xs.reserve(pts.size());
        const int n = static_cast<int>(pts.size());
        for (int y = y0; y < y1; ++y) {
            const float scan = static_cast<float>(y) + 0.5f;
            xs.clear();
            for (int i = 0, j = n - 1; i < n; j = i++) {
                const float yA = pts[static_cast<size_t>(j)].y;
                const float yB = pts[static_cast<size_t>(i)].y;
                if ((yA > scan) == (yB > scan)) {
                    continue;
                }
                const float dy = yB - yA;
                if (std::fabs(dy) < 1e-8f) {
                    continue;
                }
                const float xA = pts[static_cast<size_t>(j)].x;
                const float xB = pts[static_cast<size_t>(i)].x;
                xs.push_back((xB - xA) * (scan - yA) / dy + xA);
            }
            std::sort(xs.begin(), xs.end());
            const float top = static_cast<float>(y);
            const float bot = static_cast<float>(y + 1);
            for (size_t k = 0; k + 1 < xs.size(); k += 2) {
                const float left = xs[k];
                const float right = xs[k + 1];
                if (right - left < 0.25f) {
                    continue;
                }
                if (useGrad && left < center.x && center.x < right) {
                    emitQuad(left, center.x, top, bot);
                    emitQuad(center.x, right, top, bot);
                }
                else if (useGrad) {
                    emitQuad(left, right, top, bot);
                }
                else {
                    rlVertex2f(left, top);
                    rlVertex2f(right, top);
                    rlVertex2f(right, bot);
                    rlVertex2f(left, top);
                    rlVertex2f(right, bot);
                    rlVertex2f(left, bot);
                }
            }
        }
        rlEnd();
        rlDrawRenderBatchActive();
        rlEnableBackfaceCulling();
    }

    void FillPolygonScreen(const std::vector<Vector2>& pts, Color color)
    {
        FillPolygonScreen(pts, color, color, Vector2{0.0f, 0.0f});
    }

    void DrawLayout(const App& app, const View& view, Rectangle canvas)
    {
        const Necklace& necklace = *app.necklace;

        const float cell = view.cell;
        const auto toScreen = [&](const BeadPos& g) { return view.ToScreen(g); };
        const float beadRadius = necklace.radius * cell;
        const float thickness = Clampf(cell * 0.08f, 1.25f, 4.5f);

        if (app.physics) {
            struct FillRegion
            {
                std::vector<Vector2> outer;
                std::vector<std::vector<Vector2>> holes;
                std::vector<Vec2d> outerGrid;
                Vec2d centroid{0.0, 0.0};
                Vector2 screenCenter{0.0f, 0.0f};
                double area = 0.0;
                bool highlighted = false;
                App::MoveHint moveHint = App::MoveHint::None;
            };
            std::vector<FillRegion> regions;
            regions.reserve(necklace.faceBeads.size() + 1);
            const auto appendRegion = [&](const FaceBead& fb, bool highlighted,
                                          App::MoveHint moveHint) {
                std::vector<Vec2d> poly;
                std::vector<std::vector<Vec2d>> holes;
                if (!app.physics->FaceBeadBoundaryPolygon(fb, &poly, &holes)) {
                    return;
                }
                FillRegion region;
                region.highlighted = highlighted;
                region.moveHint = moveHint;
                region.area = ShoelaceArea(poly);
                region.outerGrid = poly;
                region.centroid = Vec2d{static_cast<double>(fb.pos.x),
                                        static_cast<double>(fb.pos.y)};
                region.screenCenter =
                    view.ToScreen(BeadPos{fb.pos.x, fb.pos.y});
                region.outer.reserve(poly.size());
                for (const Vec2d& p : poly) {
                    region.outer.push_back(
                        view.ToScreen(BeadPos{static_cast<float>(p.x), static_cast<float>(p.y)}));
                }
                for (const std::vector<Vec2d>& hole : holes) {
                    std::vector<Vector2> hs;
                    hs.reserve(hole.size());
                    for (const Vec2d& p : hole) {
                        hs.push_back(
                            view.ToScreen(BeadPos{static_cast<float>(p.x), static_cast<float>(p.y)}));
                    }
                    if (hs.size() >= 3) {
                        region.holes.push_back(std::move(hs));
                    }
                }
                regions.push_back(std::move(region));
            };
            for (const FaceBead& fb : necklace.faceBeads) {
                appendRegion(fb, fb.highlighted,
                             app.highlightReidemeister ? app.FaceBeadMoveHint(fb)
                                                       : App::MoveHint::None);
            }
            if (const FaceBead* ghost = app.physics->ShorteningHighlightBead()) {
                App::MoveHint hint = App::MoveHint::None;
                if (app.highlightReidemeister) {
                    if (app.physics->R0Active()) {
                        hint = App::MoveHint::R0;
                    }
                    else if (app.physics->R1Active()) {
                        hint = App::MoveHint::R1;
                    }
                    else if (app.physics->R2Active()) {
                        hint = App::MoveHint::R2;
                    }
                }
                appendRegion(*ghost, ghost->highlighted, hint);
            }
            const auto contains = [&](size_t outerI, size_t innerI) {
                const FillRegion& outer = regions[outerI];
                const FillRegion& inner = regions[innerI];
                return inner.area + 1e-12 < outer.area &&
                       PointInPolygon(inner.centroid, outer.outerGrid);
            };
            if (app.exteriorHighlighted) {
                DrawRectangleRec(canvas, kHighlightFill);
            }
            // Punch every finite/circle interior so exterior yellow is one
            // layer, then paint highlighted faces (punching nested/sibling
            // interiors that sit inside them), then punch unhighlighted
            // siblings that do not contain a highlighted region.
            for (const FillRegion& region : regions) {
                FillPolygonScreen(region.outer, kBackground);
            }
            std::vector<size_t> highlighted;
            highlighted.reserve(regions.size());
            for (size_t i = 0; i < regions.size(); ++i) {
                if (regions[i].highlighted) {
                    highlighted.push_back(i);
                }
            }
            std::sort(highlighted.begin(), highlighted.end(), [&](size_t a, size_t b) {
                return regions[a].area > regions[b].area;
            });
            for (const size_t i : highlighted) {
                FillPolygonScreen(regions[i].outer, kHighlightFill);
                for (size_t j = 0; j < regions.size(); ++j) {
                    if (j != i && contains(i, j)) {
                        FillPolygonScreen(regions[j].outer, kBackground);
                    }
                }
                for (const std::vector<Vector2>& hole : regions[i].holes) {
                    FillPolygonScreen(hole, kBackground);
                }
            }
            for (size_t u = 0; u < regions.size(); ++u) {
                if (regions[u].highlighted) {
                    continue;
                }
                bool coversHighlight = false;
                for (const size_t h : highlighted) {
                    if (contains(u, h)) {
                        coversHighlight = true;
                        break;
                    }
                }
                if (!coversHighlight) {
                    FillPolygonScreen(regions[u].outer, kBackground);
                }
            }
            if (app.highlightReidemeister) {
                std::vector<size_t> metal;
                metal.reserve(regions.size());
                for (size_t i = 0; i < regions.size(); ++i) {
                    if (!regions[i].highlighted && regions[i].moveHint != App::MoveHint::None) {
                        metal.push_back(i);
                    }
                }
                std::sort(metal.begin(), metal.end(), [&](size_t a, size_t b) {
                    return regions[a].area > regions[b].area;
                });
                const auto metalColors = [](App::MoveHint hint, Color* inner, Color* outer) {
                    if (hint == App::MoveHint::R2) {
                        *inner = kGoldInner;
                        *outer = kGoldOuter;
                    }
                    else if (hint == App::MoveHint::R3) {
                        *inner = kBronzeInner;
                        *outer = kBronzeOuter;
                    }
                    else {
                        *inner = kSilverInner;
                        *outer = kSilverOuter;
                    }
                };
                for (const size_t i : metal) {
                    Color inner;
                    Color outer;
                    metalColors(regions[i].moveHint, &inner, &outer);
                    FillPolygonScreen(regions[i].outer, inner, outer, regions[i].screenCenter);
                    for (size_t j = 0; j < regions.size(); ++j) {
                        if (j != i && contains(i, j)) {
                            FillPolygonScreen(regions[j].outer, kBackground);
                        }
                    }
                    for (const std::vector<Vector2>& hole : regions[i].holes) {
                        FillPolygonScreen(hole, kBackground);
                    }
                }
            }
        }

        // Vertex and edge beads underneath the multiloop (optional).
        const float outline = Clampf(beadRadius * 0.10f, 1.0f, 3.0f);
        if (app.showBeads) {
            for (const Bead& bead : necklace.beads) {
                if (bead.removed) {
                    continue;
                }
                const Vector2 c = toScreen(bead.pos);
                DrawCircleV(c, beadRadius, WHITE);
                DrawRing(c, beadRadius - outline, beadRadius, 0.0f, 360.0f, 48, kBeadOutline);
            }
        }
        // Face beads are always labeled; the disk is drawn only when pinned.
        const Color kOutsideFill = {210, 55, 50, 255};
        const Color kOutsideEdge = {140, 25, 22, 255};
        const Color kExcessHighlight = {255, 120, 0, 255};
        for (const FaceBead& faceBead : necklace.faceBeads) {
            const Vector2 c = toScreen(faceBead.pos);
            const bool outside = app.debug && app.physics &&
                                 !app.physics->FaceBeadInsideFace(faceBead);
            if (faceBead.pinned) {
                DrawCircleV(c, beadRadius, outside ? kOutsideFill : kFaceBeadFill);
                Color ringColor = outside ? kOutsideEdge : kFaceBeadEdge;
                float ringThickness = outline;
                if (!outside && app.debug && app.physics &&
                    app.physics->FaceExcessArea(faceBead) > app.physicsParams.minExcessArea) {
                    ringColor = kExcessHighlight;
                    ringThickness = Clampf(outline * 2.5f, 2.5f, 6.0f);
                }
                DrawRing(c, beadRadius - ringThickness, beadRadius, 0.0f, 360.0f, 48, ringColor);
            }
            const std::string text = std::to_string(faceBead.label);
            float fontSize = Clampf(beadRadius * 1.15f, 8.0f, 40.0f);
            float textWidth = TextWidthFont(gFontBold, text, fontSize);
            if (textWidth > beadRadius * 1.5f) {
                fontSize *= beadRadius * 1.5f / textWidth;
                textWidth = TextWidthFont(gFontBold, text, fontSize);
            }
            DrawTextFont(gFontBold, text, c.x - textWidth / 2.0f, c.y - fontSize / 2.0f,
                         fontSize, outside ? Color{80, 10, 10, 255} : Color{15, 35, 12, 255});
        }

        const auto drawEdgeBeadArc = [&](const Bead& bead, Color color) {
            if (bead.neighbors[0] < 0 || bead.neighbors[1] < 0) {
                return;
            }
            const Vector2 c = toScreen(bead.pos);
            const Vector2 n0 = toScreen(necklace.beads[bead.neighbors[0]].pos);
            const Vector2 n1 = toScreen(necklace.beads[bead.neighbors[1]].pos);
            DrawOrthogonalArc(c, beadRadius, TangentPoint(c, beadRadius, n0),
                              TangentPoint(c, beadRadius, n1), thickness, color);
        };

        if (app.unknotMode || !app.loop) {
            for (const Bead& bead : necklace.beads) {
                if (bead.removed || bead.isVertex) {
                    continue;
                }
                const int strand = (bead.strand >= 0) ? bead.strand : 0;
                drawEdgeBeadArc(bead, StrandColor(strand));
            }
            return;
        }

        const Multiloop& loop = *app.loop;

        std::vector<std::vector<Vector2>> homotopyPolys;
        std::vector<int> homotopyStrands;
        const bool homotopyActive = app.physics && app.physics->R3HomotopyActive();
        if (homotopyActive) {
            app.physics->R3HomotopyPolylines(&homotopyPolys, &homotopyStrands);
        }

        // The multiloop itself: orthogonal circle arcs within each bead,
        // drawn strand by strand so higher-index strands always lie on top.
        for (int strand = 0; strand < loop.numStrands; ++strand) {
            const Color color = StrandColor(strand);
            for (size_t bi = 0; bi < necklace.beads.size(); ++bi) {
                const Bead& bead = necklace.beads[bi];
                if (bead.removed) {
                    continue;
                }
                if (homotopyActive &&
                    app.physics->R3HomotopyHidesBead(static_cast<int>(bi))) {
                    continue;
                }
                const Vector2 c = toScreen(bead.pos);
                if (bead.isVertex) {
                    if (bead.sigmaCycle.size() < 4) {
                        continue;
                    }
                    // Slots cut by an R3 slide are -1; skip those arcs.
                    for (int j = 0; j < 2; ++j) {
                        const int edge = std::abs(bead.sigmaCycle[static_cast<size_t>(j)]);
                        if (edge < 1 || edge > loop.numEdges ||
                            loop.strandOfEdge[static_cast<size_t>(edge)] != strand) {
                            continue;
                        }
                        if (bead.neighborOfHalfEdge[j] < 0 || bead.neighborOfHalfEdge[j + 2] < 0) {
                            continue;
                        }
                        const Vector2 n0 = toScreen(necklace.beads[bead.neighborOfHalfEdge[j]].pos);
                        const Vector2 n1 = toScreen(necklace.beads[bead.neighborOfHalfEdge[j + 2]].pos);
                        DrawOrthogonalArc(c, beadRadius, TangentPoint(c, beadRadius, n0),
                                          TangentPoint(c, beadRadius, n1), thickness, color);
                    }
                    // During an R3 slide the old triangle vertices keep exactly
                    // two links; draw the through-arc so the boundary is unbroken.
                    int intact[4];
                    int intactCount = 0;
                    for (int s = 0; s < 4; ++s) {
                        if (bead.neighborOfHalfEdge[s] >= 0) {
                            intact[intactCount++] = s;
                        }
                    }
                    if (intactCount == 2) {
                        const int edge = std::abs(bead.sigmaCycle[static_cast<size_t>(intact[0])]);
                        if (edge >= 1 && edge <= loop.numEdges &&
                            loop.strandOfEdge[static_cast<size_t>(edge)] == strand) {
                            const Vector2 n0 =
                                toScreen(necklace.beads[bead.neighborOfHalfEdge[intact[0]]].pos);
                            const Vector2 n1 =
                                toScreen(necklace.beads[bead.neighborOfHalfEdge[intact[1]]].pos);
                            DrawOrthogonalArc(c, beadRadius, TangentPoint(c, beadRadius, n0),
                                              TangentPoint(c, beadRadius, n1), thickness, color);
                        }
                    }
                }
                else {
                    if (bead.edge < 1 || bead.edge > loop.numEdges ||
                        loop.strandOfEdge[bead.edge] != strand) {
                        continue;
                    }
                    if (bead.neighbors[0] < 0 || bead.neighbors[1] < 0) {
                        continue; // cut terminal during an R3 slide
                    }
                    const Vector2 n0 = toScreen(necklace.beads[bead.neighbors[0]].pos);
                    const Vector2 n1 = toScreen(necklace.beads[bead.neighbors[1]].pos);
                    DrawOrthogonalArc(c, beadRadius, TangentPoint(c, beadRadius, n0),
                                      TangentPoint(c, beadRadius, n1), thickness, color);
                }
            }
            if (homotopyActive) {
                for (size_t a = 0; a < homotopyPolys.size(); ++a) {
                    if (homotopyStrands[a] != strand) {
                        continue;
                    }
                    const std::vector<Vector2>& pts = homotopyPolys[a];
                    for (size_t k = 0; k + 1 < pts.size(); ++k) {
                        const Vector2 s0 = view.ToScreen(BeadPos{pts[k].x, pts[k].y});
                        const Vector2 s1 = view.ToScreen(BeadPos{pts[k + 1].x, pts[k + 1].y});
                        DrawLineEx(s0, s1, thickness, color);
                    }
                }
            }
        }

        for (const Bead& bead : necklace.beads) {
            if (bead.removed || bead.isVertex || bead.edge != 0) {
                continue;
            }
            const int strand = (bead.strand >= 0) ? bead.strand : 0;
            drawEdgeBeadArc(bead, StrandColor(strand));
        }

        if (app.debug) {
            // Vertex permutation: half-edge labels at each crossing bead.
            for (const Bead& bead : necklace.beads) {
                if (!bead.isVertex) {
                    continue;
                }
                const Vector2 vertex = toScreen(bead.pos);
                for (size_t i = 0; i < bead.sigmaCycle.size(); ++i) {
                    const int neighbor = bead.neighborOfHalfEdge[i];
                    if (neighbor < 0) {
                        continue;
                    }
                    const Vector2 towards = toScreen(necklace.beads[neighbor].pos);
                    DrawHalfEdgeLabel(vertex, towards, bead.sigmaCycle[i], cell);
                }
            }
            // Face permutation labels are shown in the lower-left when a face
            // bead is clicked (see Frame).
        }
    }

    // Renders the free-hand drawing canvas: committed curves, the in-progress
    // polyline, a live rubber-band segment to the cursor, and crossing dots.
    void DrawDrawingCanvas(const App& app, Rectangle canvas, Vector2 mouse)
    {
        const float thickness = 3.0f;
        const auto drawPolyline = [&](const std::vector<Vector2>& pts, bool closed, Color color) {
            for (size_t i = 0; i + 1 < pts.size(); ++i) {
                DrawLineEx(pts[i], pts[i + 1], thickness, color);
            }
            if (closed && pts.size() >= 2) {
                DrawLineEx(pts.back(), pts.front(), thickness, color);
            }
            for (const Vector2& p : pts) {
                DrawCircleV(p, thickness, color);
            }
        };

        for (size_t c = 0; c < app.drawCurves.size(); ++c) {
            drawPolyline(app.drawCurves[c], true, StrandColor(static_cast<int>(c)));
        }

        const bool inCanvas = CheckCollisionPointRec(mouse, canvas);
        if (!app.drawCurrent.empty()) {
            const Color color = StrandColor(static_cast<int>(app.drawCurves.size()));
            drawPolyline(app.drawCurrent, false, color);

            const Vector2 last = app.drawCurrent.back();
            // Closing indicator at the start vertex.
            const Vector2 start = app.drawCurrent.front();
            if (app.drawCurrent.size() >= 3) {
                DrawRing(start, App::kCloseTol - 2.0f, App::kCloseTol, 0.0f, 360.0f, 32,
                         Color{120, 180, 120, 200});
            }
            if (inCanvas) {
                const bool closing = app.drawCurrent.size() >= 3 &&
                    std::hypot(mouse.x - start.x, mouse.y - start.y) <= App::kCloseTol;
                const Vector2 target = closing ? start : mouse;
                std::string reason;
                const bool ok = closing
                    ? app.SegmentPlaceable(last, start, &reason)
                    : (app.PointPlaceable(mouse, &reason) && app.SegmentPlaceable(last, mouse, &reason));
                const Color band = ok ? Color{90, 170, 90, 220} : Color{200, 70, 70, 220};
                DrawLineEx(last, target, thickness, band);
            }
        }
        else if (inCanvas) {
            std::string reason;
            const bool ok = app.PointPlaceable(mouse, &reason);
            DrawCircleV(mouse, thickness + 1.0f,
                        ok ? Color{90, 170, 90, 220} : Color{200, 70, 70, 220});
        }

        for (const Vector2& c : app.drawCrossings) {
            DrawCircleV(c, 4.5f, Color{30, 30, 35, 255});
            DrawCircleV(c, 2.0f, WHITE);
        }
    }

    std::vector<std::string> WrapText(const Font& font, const std::string& text, float size,
                                      float maxWidth)
    {
        std::vector<std::string> lines;
        std::string line;
        std::string word;
        for (size_t i = 0; i <= text.size(); ++i) {
            const char c = (i < text.size()) ? text[i] : ' ';
            if (c != ' ') {
                word.push_back(c);
                continue;
            }
            if (word.empty()) {
                continue;
            }
            const std::string trial = line.empty() ? word : line + " " + word;
            if (TextWidthFont(font, trial, size) <= maxWidth) {
                line = trial;
            }
            else {
                if (!line.empty()) {
                    lines.push_back(line);
                }
                line = word;
            }
            word.clear();
        }
        if (!line.empty()) {
            lines.push_back(line);
        }
        return lines;
    }

    // ----------------------------------------------------------------
    // Per-frame UI.
    // ----------------------------------------------------------------

    void Frame(App& app)
    {
        const int screenWidth = GetScreenWidth();
        const int screenHeight = GetScreenHeight();
        const Vector2 mouse = GetMousePosition();
        const bool click = IsMouseButtonPressed(MOUSE_BUTTON_LEFT);
        const bool rightClick = IsMouseButtonPressed(MOUSE_BUTTON_RIGHT);

        app.game.Tick(GetTime());
        const bool inGame = app.game.Active();
        const bool gameOver = inGame && app.game.Over();
        if (inGame && IsKeyPressed(KEY_ESCAPE)) {
            if (app.game.RulesOpen()) {
                app.game.CloseRules();
            }
            else if (app.game.ExitPrompt()) {
                app.game.DismissExit();
            }
            else if (!gameOver) {
                app.game.AskExit();
            }
        }
        if (inGame && app.pendingGameFit) {
            app.pendingGameFit = false;
            app.RequestViewFit();
        }

        const float pad = 10.0f;
        const float row0 = pad;
        const float row1 = pad + 40.0f;
        const float row2 = row1 + 42.0f;
        const float row3 = row2 + 38.0f;
        const float row4 = row3 + 38.0f;
        const float gameBarH = 72.0f;
        const float canvasTop = gameOver ? 0.0f : (inGame ? gameBarH : (row4 + 36.0f));
        const Rectangle canvas = {0.0f, canvasTop, static_cast<float>(screenWidth),
                                  static_cast<float>(screenHeight) - canvasTop};
        if (IsWindowResized() && app.layout) {
            BeginViewFit(app, canvas);
        }

        BeginDrawing();
        ClearBackground(kBackground);

        const bool drawing = app.drawing;
        const Color kDisabledFill = {205, 207, 214, 255};

        const auto bigButton = [&](float bx, const std::string& text, bool enabled,
                                   Color baseColor) -> float {
            const float w = TextWidth(text, 18.0f) + 24.0f;
            const Rectangle r = {bx, row0, w, 30.0f};
            const bool hover = enabled && CheckCollisionPointRec(mouse, r);
            const Color fill = !enabled ? kDisabledFill
                                        : (hover ? ColorBrightness(baseColor, -0.15f) : baseColor);
            DrawRectangleRounded(r, 0.35f, 6, fill);
            DrawTextS(text, r.x + 12.0f, r.y + 5.0f, 18.0f, enabled ? WHITE : kBackground);
            return w;
        };
        const auto buttonClicked = [&](float bx, float w) {
            return click && CheckCollisionPointRec(mouse, Rectangle{bx, row0, w, 30.0f});
        };

        Rectangle dropdownBox = {0.0f, 0.0f, 0.0f, 0.0f};
        Rectangle allowanceBox = {0.0f, 0.0f, 0.0f, 0.0f};
        bool uiHoldsClick = false;

        if (!inGame) {
        // ---- Row 0: draw / generate ----
        float x = pad;
        if (!drawing) {
            const float w = bigButton(x, "Draw multiloop", true, kAccent);
            if (buttonClicked(x, w)) {
                app.StartDrawing();
            }
            x += w + 10.0f;
        }
        else {
            const float wd = bigButton(x, "Done", true, kAccent);
            if (buttonClicked(x, wd)) {
                app.FinishDrawing();
            }
            x += wd + 8.0f;
            const float wc = bigButton(x, "Cancel", true, kErrorRed);
            if (buttonClicked(x, wc)) {
                app.CancelDrawing();
            }
            x += wc + 10.0f;
        }
        DrawTextS("or", x, row0 + 6.0f, 18.0f, drawing ? kTextDim : kTextDark);
        x += TextWidth("or", 18.0f) + 10.0f;
        const float genBoxLeft = x;
        {
            const float w = bigButton(x, "Generate random multiloop", !drawing, kAccent);
            if (!drawing && buttonClicked(x, w)) {
                app.GenerateRandom();
            }
            x += w + 12.0f;

            const Color labelColor = drawing ? kTextDim : kTextDark;
            DrawTextS("with approximately", x, row0 + 6.0f, 18.0f, labelColor);
            x += TextWidth("with approximately", 18.0f) + 8.0f;
            const Rectangle countRect = {x, row0 + 1.0f, 56.0f, 28.0f};
            if (TextBoxFrame(app.countBox, countRect, 17.0f, mouse, click, true, !drawing) && !drawing) {
                app.GenerateRandom();
            }
            x += 56.0f + 8.0f;
            const std::string facesText = "regions (>50 may cause performance issues).";
            DrawTextS(facesText, x, row0 + 6.0f, 18.0f, labelColor);
            x += TextWidth(facesText, 18.0f) + 30.0f;

            // Density slider: treelike (0) .. meshlike (1).
            DrawTextS("treelike", x, row0 + 7.0f, 17.0f, kTextDim);
            x += TextWidth("treelike", 17.0f) + 10.0f;
            const float trackWidth = 170.0f;
            const Rectangle trackHit = {x - 8.0f, row0 + 2.0f, trackWidth + 16.0f, 30.0f};
            if (!drawing && click && CheckCollisionPointRec(mouse, trackHit)) {
                app.sliderDragging = true;
            }
            if (!IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
                app.sliderDragging = false;
            }
            if (app.sliderDragging && !drawing) {
                app.density = Clampf((mouse.x - x) / trackWidth, 0.0f, 1.0f);
            }
            const float trackY = row0 + 16.0f;
            const Color trackFill = drawing ? kDisabledFill : kAccent;
            DrawRectangle(static_cast<int>(x), static_cast<int>(trackY) - 2,
                          static_cast<int>(trackWidth), 4, kBorder);
            DrawRectangle(static_cast<int>(x), static_cast<int>(trackY) - 2,
                          static_cast<int>(trackWidth * app.density), 4, trackFill);
            DrawCircleV(Vector2{x + trackWidth * app.density, trackY}, 8.0f, trackFill);
            DrawCircleV(Vector2{x + trackWidth * app.density, trackY}, 4.0f, WHITE);
            x += trackWidth + 10.0f;
            DrawTextS("meshlike", x, row0 + 7.0f, 17.0f, kTextDim);
            x += TextWidth("meshlike", 17.0f) + 30.0f;
        }

        const auto checkbox = [&](float cx, float cy, const char* label, bool* value, bool enabled) {
            const Rectangle box = {cx, cy + 5.0f, 20.0f, 20.0f};
            const Rectangle hit = {cx, cy, 28.0f + TextWidth(label, 18.0f), 30.0f};
            if (enabled && click && CheckCollisionPointRec(mouse, hit)) {
                *value = !*value;
            }
            DrawRectangleRec(box, enabled ? WHITE : Color{225, 226, 230, 255});
            DrawRectangleLinesEx(box, 2.0f, kBorder);
            if (*value) {
                DrawRectangle(static_cast<int>(box.x) + 4, static_cast<int>(box.y) + 4, 12, 12,
                              enabled ? kAccent : Color{170, 172, 180, 255});
            }
            DrawTextS(label, box.x + 28.0f, cy + 6.0f, 18.0f, enabled ? kTextDark : kTextDim);
            return hit.width;
        };

        // Strand options, to the right of the slider (mutually exclusive).
        {
            const bool wasOneStrand = app.forceOneStrand;
            x += checkbox(x, row0, "force one strand", &app.forceOneStrand, !drawing) + 16.0f;
            if (app.forceOneStrand && !wasOneStrand) {
                app.forceSimple = false;
            }
            const bool wasSimple = app.forceSimple;
            x += checkbox(x, row0, "force simple", &app.forceSimple, !drawing) + 16.0f;
            if (app.forceSimple && !wasSimple) {
                app.forceOneStrand = false;
            }
        }

        // Thin outline grouping the random-generation controls, matching the
        // button height and clearing the last checkbox label.
        {
            const Rectangle group = {genBoxLeft, row0, x - 8.0f - genBoxLeft, 30.0f};
            DrawRectangleRoundedLinesEx(group, 0.5f, 8, 1.5f,
                                        drawing ? Color{198, 200, 208, 255} : kAccent);
        }

        // ---- Row 1: sigma entry ----
        const float sigmaLabelWidth = TextWidth("Vertex permutation in cycle notation:", 20.0f);
        DrawTextS("Vertex permutation in cycle notation:", pad, row1 + 7.0f, 20.0f,
                  drawing ? kTextDim : kTextDark);
        const Rectangle sigmaRect = {pad + sigmaLabelWidth + 12.0f, row1,
                                     screenWidth - (pad + sigmaLabelWidth + 12.0f) - pad, 34.0f};
        if (TextBoxFrame(app.sigmaBox, sigmaRect, 19.0f, mouse, click, false, !drawing)) {
            app.SubmitSigma();
        }

        // ---- Row 2: infinite face dropdown ----
        const std::string facePrompt = "Choose infinite face:";
        DrawTextS(facePrompt, pad, row2 + 4.0f, 20.0f, drawing ? kTextDim : kTextDark);
        dropdownBox = {pad + TextWidth(facePrompt, 20.0f) + 12.0f, row2, 360.0f, 30.0f};
        if (!drawing && click && CheckCollisionPointRec(mouse, dropdownBox) && !app.faceLabels.empty()) {
            app.dropdownOpen = !app.dropdownOpen;
            if (app.dropdownOpen && app.selectedFace >= 0) {
                app.dropdownScroll = static_cast<float>(std::max(0, app.selectedFace - 3));
            }
        }
        if (IsKeyPressed(KEY_ESCAPE)) {
            app.dropdownOpen = false;
            app.game.CloseAllowanceMenu();
        }
        {
            DrawRectangleRec(dropdownBox, app.faceLabels.empty() ? kBackground : WHITE);
            DrawRectangleLinesEx(dropdownBox, 2.0f, app.dropdownOpen ? kAccent : kBorder);
            std::string shown;
            Color shownColor = kTextDark;
            if (app.faceLabels.empty()) {
                shown = "(enter sigma first)";
                shownColor = kTextDim;
            }
            else if (app.selectedFace < 0) {
                shown = "(select a face)";
                shownColor = kTextDim;
            }
            else {
                shown = app.faceLabels[app.selectedFace];
            }
            BeginScissorMode(static_cast<int>(dropdownBox.x) + 4, static_cast<int>(dropdownBox.y),
                             static_cast<int>(dropdownBox.width) - 28, static_cast<int>(dropdownBox.height));
            DrawTextS(shown, dropdownBox.x + 8.0f, dropdownBox.y + 5.0f, 17.0f, shownColor);
            EndScissorMode();
            // Arrow.
            const float ax = dropdownBox.x + dropdownBox.width - 18.0f;
            const float ay = dropdownBox.y + (app.dropdownOpen ? 18.0f : 12.0f);
            const float dy = app.dropdownOpen ? -6.0f : 6.0f;
            DrawTriangle(Vector2{ax - 6.0f, ay}, Vector2{ax + 6.0f, ay}, Vector2{ax, ay + dy},
                         kTextDim);
        }
        // Strand count, display checkboxes, then the error text, to the
        // right of the dropdown.
        float rx = dropdownBox.x + dropdownBox.width + 16.0f;
        if (app.loop || app.unknotMode) {
            const int circleCount = app.CircleCount();
            const int crossingComps = app.CrossingComponentCount();
            const int numStrands = app.unknotMode ? std::max(1, circleCount)
                                                  : app.DrawnStrandCount();
            const int numRegions =
                app.unknotMode
                    ? std::max(1, circleCount)
                    : (static_cast<int>(app.loop->faces.size()) - crossingComps + circleCount);
            const std::string strandText = "(" + std::to_string(numStrands) +
                "-strand multiloop with " + std::to_string(numRegions) + " regions)";
            DrawTextS(strandText, rx, row2 + 6.0f, 18.0f, kTextDim);
            rx += TextWidth(strandText, 18.0f) + 16.0f;
        }
        rx += checkbox(rx, row2, "show beads", &app.showBeads, !drawing) + 16.0f;
#ifdef DEBUG_UI
        rx += checkbox(rx, row2, "debug", &app.debug, !drawing) + 16.0f;
        if (app.debug) {
            rx += checkbox(rx, row2, "timers", &app.showDebugTiming, !drawing) + 16.0f;
            rx += checkbox(rx, row2, "excess areas", &app.showDebugExcessAreas, !drawing) + 16.0f;
            rx += checkbox(rx, row2, "global collisions", &app.showGlobalCollisions, !drawing) +
                  16.0f;
        }
#endif
        if (!app.error.empty()) {
            DrawTextS(app.error, rx, row2 + 6.0f, 18.0f, kErrorRed);
        }

        // ---- Row 3: auto-tighten ----
        {
            const bool active = app.physics && app.physics->AutoTightening();
            const bool enabled = !drawing && app.physics != nullptr;
            const std::string label = active ? "Stop auto-tighten" : "Auto-tighten";
            const float w = TextWidth(label, 18.0f) + 24.0f;
            const Rectangle r = {pad, row3, w, 30.0f};
            const bool hover = enabled && CheckCollisionPointRec(mouse, r);
            const Color base = active ? kErrorRed : kAccent;
            const Color fill = !enabled ? kDisabledFill
                                        : (hover ? ColorBrightness(base, -0.15f) : base);
            DrawRectangleRounded(r, 0.35f, 6, fill);
            DrawTextS(label, r.x + 12.0f, r.y + 5.0f, 18.0f, enabled ? WHITE : kBackground);
            if (enabled && click && CheckCollisionPointRec(mouse, r)) {
                if (active) {
                    app.physics->StopAutoTighten();
                }
                else {
                    app.physics->StartAutoTighten();
                }
            }
            const std::string resetLabel = "Reset";
            const float resetW = TextWidth(resetLabel, 18.0f) + 24.0f;
            const float resetX = pad + w + 8.0f;
            const Rectangle resetR = {resetX, row3, resetW, 30.0f};
            const bool resetEnabled = enabled;
            const bool resetHover = resetEnabled && CheckCollisionPointRec(mouse, resetR);
            DrawRectangleRounded(resetR, 0.35f, 6,
                                 !resetEnabled ? kDisabledFill
                                               : (resetHover ? ColorBrightness(kAccent, -0.15f)
                                                             : kAccent));
            DrawTextS(resetLabel, resetR.x + 12.0f, resetR.y + 5.0f, 18.0f,
                      resetEnabled ? WHITE : kBackground);
            if (resetEnabled && click && CheckCollisionPointRec(mouse, resetR)) {
                app.ResetToGrid();
            }
            const std::string fitLabel = "Zoom to fit";
            const float fitW = TextWidth(fitLabel, 18.0f) + 24.0f;
            const float fitX = resetX + resetW + 8.0f;
            const Rectangle fitR = {fitX, row3, fitW, 30.0f};
            const bool fitEnabled = !drawing && app.layout && app.necklace;
            const bool fitHover = fitEnabled && CheckCollisionPointRec(mouse, fitR);
            DrawRectangleRounded(fitR, 0.35f, 6,
                                 !fitEnabled ? kDisabledFill
                                             : (fitHover ? ColorBrightness(kAccent, -0.15f)
                                                         : kAccent));
            DrawTextS(fitLabel, fitR.x + 12.0f, fitR.y + 5.0f, 18.0f,
                      fitEnabled ? WHITE : kBackground);
            if (fitEnabled && click && CheckCollisionPointRec(mouse, fitR)) {
                BeginViewFit(app, canvas);
            }
            const std::string undoLabel = "Undo";
            const float undoW = TextWidth(undoLabel, 18.0f) + 24.0f;
            const float undoX = fitX + fitW + 8.0f;
            const Rectangle undoR = {undoX, row3, undoW, 30.0f};
            const bool undoEnabled = !drawing && !app.r3Undo.empty();
            const bool undoHover = undoEnabled && CheckCollisionPointRec(mouse, undoR);
            DrawRectangleRounded(undoR, 0.35f, 6,
                                 !undoEnabled ? kDisabledFill
                                              : (undoHover ? ColorBrightness(kAccent, -0.15f)
                                                           : kAccent));
            DrawTextS(undoLabel, undoR.x + 12.0f, undoR.y + 5.0f, 18.0f,
                      undoEnabled ? WHITE : kBackground);
            if (undoEnabled && click && CheckCollisionPointRec(mouse, undoR)) {
                app.UndoR3();
            }
            const std::string pinAllLabel = "Pin all";
            const float pinAllW = TextWidth(pinAllLabel, 18.0f) + 24.0f;
            const float pinAllX = undoX + undoW + 8.0f;
            const Rectangle pinAllR = {pinAllX, row3, pinAllW, 30.0f};
            const bool pinEnabled = !drawing && app.necklace != nullptr;
            const bool pinAllHover = pinEnabled && CheckCollisionPointRec(mouse, pinAllR);
            DrawRectangleRounded(pinAllR, 0.35f, 6,
                                 !pinEnabled ? kDisabledFill
                                             : (pinAllHover ? ColorBrightness(kAccent, -0.15f)
                                                            : kAccent));
            DrawTextS(pinAllLabel, pinAllR.x + 12.0f, pinAllR.y + 5.0f, 18.0f,
                      pinEnabled ? WHITE : kBackground);
            if (pinEnabled && click && CheckCollisionPointRec(mouse, pinAllR)) {
                app.SetAllPinned(true);
            }
            const std::string unpinAllLabel = "Unpin all";
            const float unpinAllW = TextWidth(unpinAllLabel, 18.0f) + 24.0f;
            const float unpinAllX = pinAllX + pinAllW + 8.0f;
            const Rectangle unpinAllR = {unpinAllX, row3, unpinAllW, 30.0f};
            const bool unpinAllHover = pinEnabled && CheckCollisionPointRec(mouse, unpinAllR);
            DrawRectangleRounded(unpinAllR, 0.35f, 6,
                                 !pinEnabled ? kDisabledFill
                                             : (unpinAllHover ? ColorBrightness(kAccent, -0.15f)
                                                              : kAccent));
            DrawTextS(unpinAllLabel, unpinAllR.x + 12.0f, unpinAllR.y + 5.0f, 18.0f,
                      pinEnabled ? WHITE : kBackground);
            if (pinEnabled && click && CheckCollisionPointRec(mouse, unpinAllR)) {
                app.SetAllPinned(false);
            }
            const float reidX = unpinAllX + unpinAllW + 16.0f;
            const float reidW =
                checkbox(reidX, row3, "Highlight Reidemeister regions",
                         &app.highlightReidemeister, !drawing);
            const float afterButtonsX = reidX + reidW + 16.0f;
            if (app.physics) {
                char buf[128];
                std::snprintf(buf, sizeof(buf),
                              "tension: %.2f   (stretch %.2f  +  compression %.2f)",
                              app.physics->TotalTension(), app.physics->StretchTension(),
                              app.physics->CompressionTension());
                DrawTextS(buf, afterButtonsX, row3 + 6.0f, 18.0f, kTextDim);
                if (app.physics->Tangled()) {
                    const float tangleX = afterButtonsX + TextWidth(buf, 18.0f) + 12.0f;
                    DrawTextS("tangle", tangleX, row3 + 6.0f, 18.0f, kErrorRed);
                }
            }
            else {
                DrawTextS("Drag beads to perturb; auto-tighten removes slack.", afterButtonsX,
                          row3 + 6.0f, 18.0f, kTextDim);
            }
            if (!app.debug) {
                app.tuningUI.open = false;
            }
            if (app.debug) {
                const std::string tuneLabel =
                    app.tuningUI.open ? "Hide physics params" : "Physics params";
                const float tuneW = TextWidth(tuneLabel, 18.0f) + 24.0f;
                const float tuneX = static_cast<float>(screenWidth) - tuneW - pad;
                const Rectangle tuneBtn = {tuneX, row3, tuneW, 30.0f};
                const bool tuneHover = CheckCollisionPointRec(mouse, tuneBtn);
                DrawRectangleRounded(tuneBtn, 0.35f, 6,
                                     tuneHover ? ColorBrightness(kAccent, -0.15f) : kAccent);
                DrawTextS(tuneLabel, tuneBtn.x + 12.0f, tuneBtn.y + 5.0f, 18.0f, WHITE);
                if (click && CheckCollisionPointRec(mouse, tuneBtn)) {
                    app.tuningUI.open = !app.tuningUI.open;
                }
            }
        }

        // ---- Row 4: start a match on the open diagram ----
        {
            const std::string playLabel = "Play the Unpinning Game";
            const float playFont = 18.0f;
            const float playW = TextWidthFont(gFontBold, playLabel, playFont) + 24.0f;
            const Rectangle playR = {pad, row4, playW, 30.0f};
            const bool diagram = app.layout != nullptr && app.necklace != nullptr;
            const bool busy = drawing || (app.physics && (app.physics->MoveActive() ||
                                                          app.physics->AutoTightening() ||
                                                          app.physics->Dragging()));
            const bool playEnabled = diagram && !busy;
            const bool playHover = playEnabled && CheckCollisionPointRec(mouse, playR);
            const Color playBase = {46, 158, 78, 255};
            const Color playFill = !playEnabled ? kDisabledFill
                                                : (playHover ? ColorBrightness(playBase, -0.15f)
                                                             : playBase);
            DrawRectangleRounded(playR, 0.35f, 6, playFill);
            DrawTextFont(gFontBold, playLabel, playR.x + 12.0f, playR.y + 5.0f, playFont,
                         playEnabled ? WHITE : kBackground);
            if (playEnabled && click && CheckCollisionPointRec(mouse, playR)) {
                app.game.CloseAllowanceMenu();
                app.BeginUnpinningMatch();
            }

            float ax = pad + playW + 12.0f;
            const std::string lead = "on this multiloop with an allowance of";
            DrawTextS(lead, ax, row4 + 6.0f, 18.0f, kTextDark);
            ax += TextWidth(lead, 18.0f) + 8.0f;
            allowanceBox = {ax, row4, 64.0f, 30.0f};
            const bool menu = app.game.AllowanceMenuOpen();
            if (click && CheckCollisionPointRec(mouse, allowanceBox)) {
                app.dropdownOpen = false;
                app.game.ToggleAllowanceMenu();
            }
            DrawRectangleRec(allowanceBox, WHITE);
            DrawRectangleLinesEx(allowanceBox, 2.0f, menu ? kAccent : kBorder);
            const std::string allowanceShown =
                UnpinningGame::AllowanceLabel(app.game.AllowanceIndex());
            DrawTextS(allowanceShown, allowanceBox.x + 10.0f, allowanceBox.y + 5.0f, 18.0f,
                      kTextDark);
            const float arrowX = allowanceBox.x + allowanceBox.width - 16.0f;
            const float arrowY = allowanceBox.y + (menu ? 18.0f : 12.0f);
            const float arrowDy = menu ? -6.0f : 6.0f;
            DrawTriangle(Vector2{arrowX - 5.0f, arrowY}, Vector2{arrowX + 5.0f, arrowY},
                         Vector2{arrowX, arrowY + arrowDy}, kTextDim);
            ax += allowanceBox.width + 8.0f;
            DrawTextS("Reidemeister moves per turn.", ax, row4 + 6.0f, 18.0f, kTextDark);
        }
        if (app.game.AllowanceMenuOpen()) {
            const Rectangle list = {allowanceBox.x, allowanceBox.y + allowanceBox.height,
                                    allowanceBox.width, 6.0f * 26.0f};
            if (CheckCollisionPointRec(mouse, list)) {
                uiHoldsClick = true;
            }
        }
        if (app.dropdownOpen && !app.faceLabels.empty()) {
            const int faceCount = static_cast<int>(app.faceLabels.size());
            const int visibleFaces = std::min(faceCount, 14);
            const Rectangle faceList = {dropdownBox.x, dropdownBox.y + dropdownBox.height,
                                        dropdownBox.width, visibleFaces * 26.0f};
            if (CheckCollisionPointRec(mouse, faceList)) {
                uiHoldsClick = true;
            }
        }
        }
        else if (!gameOver) {
            const float infoSize = 22.0f;
            const float titleSize = infoSize * 2.0f;
            const float infoY = (gameBarH - infoSize) * 0.5f;
            const std::string title = app.game.TurnTitle();
            const float titleW = TextWidthFont(gFontBold, title, titleSize);
            const float boxPadX = 18.0f;
            const float boxPadY = 8.0f;
            const float boxW = titleW + boxPadX * 2.0f;
            const float boxH = titleSize + boxPadY * 2.0f;
            const float boxX = (static_cast<float>(screenWidth) - boxW) * 0.5f;
            const float boxY = (gameBarH - boxH) * 0.5f;
            const bool playerTwo = app.game.CurrentPlayer() == 2;
            const Color boxFill = playerTwo ? BLACK : WHITE;
            const Color titleColor = playerTwo ? WHITE : BLACK;
            const Rectangle titleBox = {boxX, boxY, boxW, boxH};
            DrawRectangleRec(titleBox, boxFill);
            if (!playerTwo) {
                DrawRectangleLinesEx(titleBox, 2.0f, BLACK);
            }
            const float titleX = boxX + boxPadX;
            const float titleY = boxY + boxPadY;

            const float rulesFont = 18.0f;
            const std::string rulesLabel = "Rules";
            const float rulesW = TextWidth(rulesLabel, rulesFont) + 24.0f;
            const float rulesH = 32.0f;
            const Rectangle rulesR = {static_cast<float>(screenWidth) - pad - rulesW,
                                      (gameBarH - rulesH) * 0.5f, rulesW, rulesH};
            const bool rulesEnabled = !app.game.RulesOpen() && !app.game.ExitPrompt();
            const bool rulesHover = rulesEnabled && CheckCollisionPointRec(mouse, rulesR);
            DrawRectangleRounded(rulesR, 0.35f, 6,
                                 rulesHover ? ColorBrightness(kAccent, -0.15f) : kAccent);
            DrawTextS(rulesLabel, rulesR.x + 12.0f, rulesR.y + 6.0f, rulesFont, WHITE);
            if (rulesHover && click) {
                app.game.OpenRules();
            }

            const std::string scores = app.game.ScoreLine();
            const float scoreW = TextWidthFont(gFontBold, scores, infoSize);
            const float scoreRight = rulesR.x - 16.0f;
            const float scoreX = scoreRight - scoreW;
            const float scoreClipL = boxX + boxW + 12.0f;
            if (scoreX > scoreClipL) {
                DrawTextFont(gFontBold, scores, scoreX, infoY, infoSize, kTextDark);
            }
            else if (scoreRight - scoreClipL > 40.0f) {
                BeginScissorMode(static_cast<int>(scoreClipL), 0,
                                 static_cast<int>(scoreRight - scoreClipL),
                                 static_cast<int>(gameBarH));
                DrawTextFont(gFontBold, scores, scoreClipL, infoY, infoSize, kTextDark);
                EndScissorMode();
            }

            DrawTextFont(gFontBold, title, titleX, titleY, titleSize, titleColor);

            const std::string detail = app.game.StatusDetail(app.AnyPinned(), app.HasPlayableMove());
            const float leftLimit = std::max(40.0f, boxX - pad - 16.0f);
            float detailSize = infoSize;
            const float detailW = TextWidthFont(gFontBold, detail, infoSize);
            if (detailW > leftLimit) {
                detailSize = std::max(14.0f, infoSize * leftLimit / detailW);
            }
            const float detailY = (gameBarH - detailSize) * 0.5f;
            BeginScissorMode(static_cast<int>(pad), 0, static_cast<int>(leftLimit),
                             static_cast<int>(gameBarH));
            DrawTextFont(gFontBold, detail, pad, detailY, detailSize, kTextDark);
            EndScissorMode();
            if (app.game.InTangle() && !app.game.RulesOpen() && !app.game.ExitPrompt()) {
                float bx = pad + TextWidthFont(gFontBold, detail, detailSize) + 16.0f;
                const auto tangleButton = [&](const std::string& label, bool enabled) -> bool {
                    const float w = TextWidth(label, 18.0f) + 24.0f;
                    if (bx + w > boxX - 8.0f) {
                        return false;
                    }
                    const Rectangle r = {bx, (gameBarH - 32.0f) * 0.5f, w, 32.0f};
                    const bool hover = enabled && CheckCollisionPointRec(mouse, r);
                    const Color fill = !enabled ? kDisabledFill
                                                : (hover ? ColorBrightness(kAccent, -0.15f) : kAccent);
                    DrawRectangleRounded(r, 0.35f, 6, fill);
                    DrawTextS(label, r.x + 12.0f, r.y + 6.0f, 18.0f, enabled ? WHITE : kBackground);
                    bx += w + 8.0f;
                    return enabled && click && CheckCollisionPointRec(mouse, r);
                };
                if (tangleButton("Undo", !app.r3Undo.empty())) {
                    app.UndoR3();
                    app.game.OnTangleUndone();
                }
                if (tangleButton("Reset", app.physics != nullptr)) {
                    app.ResolveTangleReset();
                }
            }
        }

        // ---- Canvas ----
        if (drawing) {
            if (click && !uiHoldsClick && CheckCollisionPointRec(mouse, canvas)) {
                const double now = GetTime();
                const bool doubleClick = (now - app.lastCanvasClick) < App::kDoubleClickTime;
                app.lastCanvasClick = now;
                if (doubleClick && !app.drawCurrent.empty()) {
                    app.HandleDoubleClick(mouse);
                }
                else {
                    app.HandleDrawClick(mouse);
                }
            }
            if (IsKeyPressed(KEY_DELETE) || IsKeyPressed(KEY_BACKSPACE)) {
                app.DeleteLastPoint();
            }
            if (IsKeyPressed(KEY_ENTER) || IsKeyPressed(KEY_KP_ENTER)) {
                app.FinishDrawing();
            }
            if (IsKeyPressed(KEY_ESCAPE)) {
                app.CancelDrawing();
            }
            DrawDrawingCanvas(app, canvas, mouse);
            const std::string hint =
                "Click to place vertices. Click the start vertex or double-click to close a loop. "
                "Delete removes the last point. Enter or Done finishes. Escape cancels.";
            const float hw = TextWidth(hint, 18.0f);
            DrawTextS(hint, (screenWidth - hw) / 2.0f, canvas.y + 6.0f, 18.0f, kTextDim);
        }
        else if (app.layout && app.necklace) {
            const View view = MakeView(app, *app.layout, canvas);
            const bool tuningBlocks = PhysicsTuningBlocksInput(app.tuningUI, mouse);
            if (app.physics) {
                const Vector2 gridMouse = view.ToGrid(mouse);
                const bool blockPlay = inGame && (app.game.ExitPrompt() || app.game.Over() ||
                                                   app.game.RulesOpen());
                const bool allowDrag = (!inGame || app.game.AllowsDrag()) && !app.game.RulesOpen();
                // Right-click toggles the innermost region under the cursor.
                if (rightClick && !blockPlay && !app.dropdownOpen && !tuningBlocks &&
                    CheckCollisionPointRec(mouse, canvas)) {
                    const int region = app.physics->PickRegionAt(gridMouse);
                    if (region < 0) {
                        app.exteriorHighlighted = !app.exteriorHighlighted;
                    }
                    else if (region < static_cast<int>(app.necklace->faceBeads.size())) {
                        FaceBead& fb = app.necklace->faceBeads[static_cast<size_t>(region)];
                        fb.highlighted = !fb.highlighted;
                    }
                }
                // Click and drag a bead (unless interacting with the dropdown).
                // A match allows dragging only while a tangle is being resolved.
                if (allowDrag && click && !uiHoldsClick && !blockPlay && !app.dropdownOpen &&
                    !tuningBlocks && CheckCollisionPointRec(mouse, canvas) &&
                    !app.physics->Dragging()) {
                    app.pickPressedOnCanvas = true;
                    app.pickClickPos = mouse;
                    const float pickDist =
                        std::max(app.necklace->radius, 14.0f / view.cell);
                    const int body = app.physics->PickBody(gridMouse, pickDist);
                    if (body >= 0) {
                        app.pickClickBody = body;
                        app.physics->BeginDrag(body, gridMouse);
                    }
                    else {
                        app.pickClickBody = -1;
                    }
                }
                else if (!allowDrag && click && !uiHoldsClick && !blockPlay && !app.dropdownOpen &&
                         !tuningBlocks && CheckCollisionPointRec(mouse, canvas)) {
                    const float pickDist =
                        std::max(app.necklace->radius, 14.0f / view.cell);
                    app.HandleGameClick(gridMouse, pickDist);
                }
                if (app.physics->Dragging() && !allowDrag) {
                    app.pickClickBody = -1;
                    app.pickPressedOnCanvas = false;
                    app.physics->EndDrag();
                    app.physics->RefreshTangle();
                    if (inGame && !app.physics->Tangled()) {
                        app.game.OnDragCleared(app.AnyPinned(), GetTime());
                    }
                }
                else if (app.physics->Dragging()) {
                    if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
                        app.physics->Drag(gridMouse);
                    }
                    else {
                        const int beadCount = static_cast<int>(app.necklace->beads.size());
                        const float dx = mouse.x - app.pickClickPos.x;
                        const float dy = mouse.y - app.pickClickPos.y;
                        const bool isClick = dx * dx + dy * dy < 64.0f;
                        const bool pickedFace =
                            app.pickClickBody >= beadCount &&
                            app.pickClickBody <
                                beadCount + static_cast<int>(app.necklace->faceBeads.size());
                        const int clickedFaceBead =
                            pickedFace ? app.pickClickBody - beadCount : -1;
                        app.pickClickBody = -1;
                        app.pickPressedOnCanvas = false;
                        app.physics->EndDrag();
                        app.physics->RefreshTangle();
                        if (inGame && app.game.InTangle() && !app.physics->Tangled()) {
                            if (app.physics->AutoTightening()) {
                                app.physics->StopAutoTighten();
                            }
                            app.game.OnDragCleared(app.AnyPinned(), GetTime());
                        }
                        if (isClick && clickedFaceBead >= 0 && !inGame) {
                            FaceBead& fb = app.necklace->faceBeads[static_cast<size_t>(
                                clickedFaceBead)];
                            if (app.debug) {
                                app.debugFaceBeadIndex = clickedFaceBead;
                                app.physics->PrintShoelaceDebug(fb.face);
                            }
                            if (fb.pinned) {
                                app.PushR3Undo();
                                fb.pinned = false;
                            }
                            else {
                                app.TryStartReidemeisterAt(gridMouse);
                            }
                        }
                    }
                }
                else if (app.pickPressedOnCanvas && !IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
                    const float dx = mouse.x - app.pickClickPos.x;
                    const float dy = mouse.y - app.pickClickPos.y;
                    const bool isClick = dx * dx + dy * dy < 64.0f;
                    app.pickPressedOnCanvas = false;
                    app.pickClickBody = -1;
                    if (isClick && !inGame && !uiHoldsClick && !app.dropdownOpen && !tuningBlocks) {
                        app.TryStartReidemeisterAt(gridMouse);
                    }
                }
                app.physics->SetGlobalCollisions(app.debug && app.showGlobalCollisions);
                const bool debugTiming = app.debug && app.showDebugTiming && app.physics;
                if (debugTiming) {
                    app.physics->ResetDebugTimings();
                }
                const auto physicsT0 = std::chrono::high_resolution_clock::now();
                const auto noteGameRewrite = [&](bool ok) {
                    if (!app.game.MovePending()) {
                        return;
                    }
                    if (!ok) {
                        app.game.OnMoveAborted();
                        app.UndoR3();
                    }
                    else {
                        app.game.OnRewriteFinished(GetTime());
                    }
                };
                if (app.physics->MoveActive()) {
                    app.physics->MoveStep();
                    if (app.physics->R0Complete()) {
                        const App::MoveFinish fin = app.CompleteR0Move();
                        if (fin == App::MoveFinish::Cleared) {
                            if (app.game.Active()) {
                                app.game.OnBoardCleared(GetTime());
                            }
                        }
                        else {
                            noteGameRewrite(fin == App::MoveFinish::Done);
                        }
                    }
                    else if (app.physics->R1Complete()) {
                        noteGameRewrite(app.CompleteR1Move());
                    }
                    else if (app.physics->R2Complete()) {
                        noteGameRewrite(app.CompleteR2Move());
                    }
                    else if (app.physics->R3Complete()) {
                        noteGameRewrite(app.CompleteR3Move());
                    }
                }
                else if (app.physics->AutoTightening()) {
                    app.physics->AutoStep();
                }
                else if (app.physics) {
                    app.physics->IdleStep();
                }
                const bool autoNow = app.physics && app.physics->AutoTightening();
                const bool r3Now = app.physics && app.physics->MoveActive();
                const bool wentIdle =
                    (app.wasAutoTightening && !autoNow) || (app.wasR3Active && !r3Now);
                const bool settleNow =
                    wentIdle || (app.game.NeedsSettle() && !autoNow && !r3Now);
                if (settleNow) {
                    const bool searchRotation = app.wasAutoTightening && !autoNow;
                    BeginViewFit(app, canvas, searchRotation);
                    if (app.physics) {
                        app.physics->RefreshTangle();
                    }
                    if (app.game.Active() && !app.game.Over()) {
                        const bool tangled = app.physics && app.physics->Tangled();
                        app.game.OnSettled(tangled, app.AnyPinned(), GetTime());
                    }
                }
                app.wasAutoTightening = autoNow;
                app.wasR3Active = r3Now;
                if (debugTiming && app.physics) {
                    const auto physicsT1 = std::chrono::high_resolution_clock::now();
                    const double frameMs =
                        std::chrono::duration<double, std::milli>(physicsT1 - physicsT0).count();
                    app.physics->SetDebugFrameTotalMs(frameMs);
                }
            }
            if (app.necklace && app.layout) {
                DrawLayout(app, view, canvas);
            }

            if (!inGame && app.debug && app.physics && app.necklace) {
                const float fontSize = 15.0f;
                const float lineH = 18.0f;
                const float panelPad = 10.0f;
                const float canvasMargin = 12.0f;
                const float panelGap = 16.0f;
                const Color panelFill = Color{252, 252, 255, 235};
                const float panelX = canvas.x + canvasMargin;

                float reservedTop = canvas.y + canvasMargin;
                if (app.showDebugTiming) {
                    const std::vector<std::string> timingLines =
                        BuildPhysicsTimingLines(app.physics->DebugTimings(),
                                              app.physics->AutoTightening());
                    DrawDebugTextPanel(panelX, reservedTop, timingLines, fontSize, lineH, panelPad,
                                       panelFill);
                    reservedTop += DebugPanelHeight(static_cast<int>(timingLines.size()), lineH,
                                                    panelPad) +
                                   panelGap;
                }

                if (app.showDebugExcessAreas) {
                    std::vector<std::string> lines;
                    lines.push_back("Excess areas:");
                    for (const FaceBead& fb : app.necklace->faceBeads) {
                        char buf[48];
                        std::snprintf(buf, sizeof(buf), "  %d: %.3f",
                                      fb.label, app.physics->FaceExcessArea(fb));
                        lines.push_back(buf);
                    }
                    if (app.debugFaceBeadIndex >= 0 &&
                        app.debugFaceBeadIndex <
                            static_cast<int>(app.necklace->faceBeads.size())) {
                        const FaceBead& sel =
                            app.necklace->faceBeads[static_cast<size_t>(app.debugFaceBeadIndex)];
                        lines.push_back("phi: " + CycleToString(sel.phiCycle));
                    }
                    const float maxHeight =
                        canvas.y + canvas.height - canvasMargin - reservedTop;
                    const float fullHeight =
                        DebugPanelHeight(static_cast<int>(lines.size()), lineH, panelPad);
                    const float boxH = (maxHeight > lineH + 2.0f * panelPad)
                                           ? std::min(fullHeight, maxHeight)
                                           : fullHeight;
                    const float by = canvas.y + canvas.height - canvasMargin - boxH;
                    DrawDebugTextPanel(panelX, by, lines, fontSize, lineH, panelPad, panelFill,
                                       maxHeight);
                }
            }
        }
        else if (!app.game.Over()) {
            const std::string hint = app.loop
                ? "Choose the infinite face from the dropdown above."
                : "Draw a multiloop, generate one, or enter a vertex permutation and press Enter.";
            const float hw = TextWidth(hint, 20.0f);
            DrawTextS(hint, (screenWidth - hw) / 2.0f, canvas.y + canvas.height / 2.0f - 10.0f,
                      20.0f, kTextDim);
        }

        // ---- Dropdown list overlay (drawn last, on top of the canvas) ----
        if (!inGame && app.dropdownOpen && !app.faceLabels.empty()) {
            const float itemHeight = 26.0f;
            const int count = static_cast<int>(app.faceLabels.size());
            const int visible = std::min(count, 14);
            const Rectangle list = {dropdownBox.x, dropdownBox.y + dropdownBox.height,
                                    dropdownBox.width, visible * itemHeight};

            app.dropdownScroll -= GetMouseWheelMove() * 2.0f;
            app.dropdownScroll = Clampf(app.dropdownScroll, 0.0f, static_cast<float>(count - visible));
            const int first = static_cast<int>(app.dropdownScroll);

            DrawRectangleRec(list, WHITE);
            DrawRectangleLinesEx(list, 2.0f, kAccent);
            BeginScissorMode(static_cast<int>(list.x), static_cast<int>(list.y),
                             static_cast<int>(list.width), static_cast<int>(list.height));
            for (int row = 0; row < visible; ++row) {
                const int i = first + row;
                if (i >= count) {
                    break;
                }
                const Rectangle item = {list.x, list.y + row * itemHeight, list.width, itemHeight};
                const bool hover = CheckCollisionPointRec(mouse, item);
                if (i == app.selectedFace) {
                    DrawRectangleRec(item, kAccent);
                }
                else if (hover) {
                    DrawRectangleRec(item, kHover);
                }
                DrawTextS(app.faceLabels[i], item.x + 8.0f, item.y + 3.0f, 16.0f,
                          i == app.selectedFace ? WHITE : kTextDark);
                if (hover && click) {
                    app.SelectFace(i);
                    app.dropdownOpen = false;
                }
            }
            EndScissorMode();
            if (count > visible) {
                // Scrollbar hint.
                const float barHeight = list.height * static_cast<float>(visible) / count;
                const float barY = list.y + (list.height - barHeight) *
                                              app.dropdownScroll / static_cast<float>(count - visible);
                DrawRectangle(static_cast<int>(list.x + list.width - 6.0f), static_cast<int>(barY),
                              4, static_cast<int>(barHeight), kBorder);
            }
            if (click && !CheckCollisionPointRec(mouse, list) &&
                !CheckCollisionPointRec(mouse, dropdownBox)) {
                app.dropdownOpen = false;
            }
        }

        // ---- Drawing warning overlay near the cursor ----
        if (app.drawing && !app.drawWarning.empty() && GetTime() < app.drawWarningUntil) {
            const float fontSize = 16.0f;
            const float w = TextWidth(app.drawWarning, fontSize) + 16.0f;
            const float h = fontSize + 12.0f;
            float bx = app.drawWarningPos.x + 14.0f;
            float by = app.drawWarningPos.y + 14.0f;
            bx = Clampf(bx, 0.0f, screenWidth - w);
            by = Clampf(by, 0.0f, screenHeight - h);
            const Rectangle box = {bx, by, w, h};
            DrawRectangleRec(box, Color{60, 30, 32, 235});
            DrawRectangleLinesEx(box, 1.5f, kErrorRed);
            DrawTextS(app.drawWarning, bx + 8.0f, by + 6.0f, fontSize, Color{255, 225, 225, 255});
        }

        if (!inGame && app.game.AllowanceMenuOpen()) {
            const float itemHeight = 26.0f;
            const Rectangle list = {allowanceBox.x, allowanceBox.y + allowanceBox.height,
                                    allowanceBox.width,
                                    UnpinningGame::kAllowanceCount * itemHeight};
            DrawRectangleRec(list, WHITE);
            DrawRectangleLinesEx(list, 2.0f, kAccent);
            for (int i = 0; i < UnpinningGame::kAllowanceCount; ++i) {
                const Rectangle item = {list.x, list.y + i * itemHeight, list.width, itemHeight};
                const bool hover = CheckCollisionPointRec(mouse, item);
                const bool selected = i == app.game.AllowanceIndex();
                if (selected) {
                    DrawRectangleRec(item, kAccent);
                }
                else if (hover) {
                    DrawRectangleRec(item, kHover);
                }
                const std::string label = UnpinningGame::AllowanceLabel(i);
                DrawTextS(label, item.x + 10.0f, item.y + 4.0f, 16.0f,
                          selected ? WHITE : kTextDark);
                if (hover && click) {
                    app.game.SetAllowanceIndex(i);
                    app.game.CloseAllowanceMenu();
                }
            }
            if (click && !CheckCollisionPointRec(mouse, list) &&
                !CheckCollisionPointRec(mouse, allowanceBox)) {
                app.game.CloseAllowanceMenu();
            }
        }

        if (app.game.RulesOpen()) {
            DrawRectangle(0, 0, screenWidth, screenHeight, Color{40, 42, 48, 90});
            const float boxW = 760.0f;
            const float inset = 32.0f;
            const float titleSize = 32.0f;
            const float bodySize = 20.0f;
            const float lineH = 28.0f;
            const float textW = boxW - 2.0f * inset;
            const char* paragraphs[] = {
                "The Unpinning Game is a 2-player game played on a loop in the plane (drawn in red) which intersects itself, with pins (drawn in green) in numerically-labeled complementary regions.",
                "Your turn consists of two phases:",
                "1) Click on the numerically-labeled regions without pins, attempting to decrease the number of self-intersections. You are awarded one point for each self-intersection that you remove. You can perform moves only on regions with zero, one, two, or three adjacent intersections (Reidemeister moves).",
                "2) Remove a pin, if there is at least one remaining.",
                "The winner is the player with more points at the end.",
            };
            std::vector<std::string> bodyLines;
            for (const char* paragraph : paragraphs) {
                if (!bodyLines.empty()) {
                    bodyLines.push_back("");
                }
                const std::vector<std::string> wrapped =
                    WrapText(gFont, paragraph, bodySize, textW);
                bodyLines.insert(bodyLines.end(), wrapped.begin(), wrapped.end());
            }
            const float btnH = 36.0f;
            const float boxH = inset + titleSize + 18.0f +
                               static_cast<float>(bodyLines.size()) * lineH + 22.0f + btnH + inset;
            const Rectangle box = {(static_cast<float>(screenWidth) - boxW) * 0.5f,
                                   (static_cast<float>(screenHeight) - boxH) * 0.5f, boxW, boxH};
            DrawRectangleRounded(box, 0.04f, 8, Color{252, 252, 255, 255});
            DrawRectangleRoundedLinesEx(box, 0.04f, 8, 2.0f, kAccent);
            const std::string rulesTitle = "The Unpinning Game";
            const float rulesTitleW = TextWidthFont(gFontBold, rulesTitle, titleSize);
            float y = box.y + inset;
            DrawTextFont(gFontBold, rulesTitle, box.x + (boxW - rulesTitleW) * 0.5f, y, titleSize,
                         kTextDark);
            y += titleSize + 18.0f;
            for (const std::string& line : bodyLines) {
                if (!line.empty()) {
                    DrawTextS(line, box.x + inset, y, bodySize, kTextDark);
                }
                y += lineH;
            }
            y += 8.0f;
            const std::string backLabel = "Back";
            const float backW = TextWidth(backLabel, 18.0f) + 36.0f;
            const Rectangle backR = {box.x + (boxW - backW) * 0.5f, y, backW, btnH};
            const bool backHover = CheckCollisionPointRec(mouse, backR);
            DrawRectangleRounded(backR, 0.35f, 6,
                                 backHover ? ColorBrightness(kAccent, -0.15f) : kAccent);
            const float backTextW = TextWidth(backLabel, 18.0f);
            DrawTextS(backLabel, backR.x + (backR.width - backTextW) * 0.5f, backR.y + 8.0f, 18.0f,
                      WHITE);
            if (click && backHover) {
                app.game.CloseRules();
            }
        }

        if (app.game.ExitPrompt()) {
            DrawRectangle(0, 0, screenWidth, screenHeight, Color{40, 42, 48, 90});
            const float boxW = 420.0f;
            const float boxH = 156.0f;
            const Rectangle box = {(static_cast<float>(screenWidth) - boxW) * 0.5f,
                                   (static_cast<float>(screenHeight) - boxH) * 0.5f, boxW, boxH};
            DrawRectangleRounded(box, 0.08f, 8, Color{252, 252, 255, 255});
            DrawRectangleRoundedLinesEx(box, 0.08f, 8, 2.0f, kAccent);
            const std::string prompt = "Exit to sandbox?";
            const float promptSize = 28.0f;
            const float promptW = TextWidthFont(gFontBold, prompt, promptSize);
            DrawTextFont(gFontBold, prompt, box.x + (boxW - promptW) * 0.5f, box.y + 32.0f,
                         promptSize, kTextDark);
            const float btnW = 110.0f;
            const float btnH = 34.0f;
            const Rectangle yesR = {box.x + 64.0f, box.y + 92.0f, btnW, btnH};
            const Rectangle noR = {box.x + boxW - 64.0f - btnW, box.y + 92.0f, btnW, btnH};
            const auto promptButton = [&](Rectangle r, const std::string& label, Color base) {
                const bool hover = CheckCollisionPointRec(mouse, r);
                DrawRectangleRounded(r, 0.35f, 6, hover ? ColorBrightness(base, -0.15f) : base);
                const float tw = TextWidth(label, 18.0f);
                DrawTextS(label, r.x + (r.width - tw) * 0.5f, r.y + 7.0f, 18.0f, WHITE);
                return click && hover;
            };
            if (promptButton(yesR, "Yes", kAccent)) {
                app.ConfirmExitToSandbox();
            }
            if (promptButton(noR, "No", kErrorRed)) {
                app.game.DismissExit();
            }
        }

        if (app.game.Over()) {
            DrawRectangle(0, 0, screenWidth, screenHeight, kBackground);
            const auto pointsWord = [](int n) {
                return std::to_string(n) + (n == 1 ? std::string(" point") : std::string(" points"));
            };
            const float titleSize = 42.0f;
            const std::string title = "Game over!";
            float y = static_cast<float>(screenHeight) * 0.5f - 130.0f;
            const float titleW = TextWidthFont(gFontBold, title, titleSize);
            DrawTextFont(gFontBold, title,
                         (static_cast<float>(screenWidth) - titleW) * 0.5f, y, titleSize, kTextDark);
            y += 68.0f;
            const float lineSize = 26.0f;
            const std::string line1 = "Player 1: " + pointsWord(app.game.Score(1));
            const std::string line2 = "Player 2: " + pointsWord(app.game.Score(2));
            const float line1W = TextWidthFont(gFontBold, line1, lineSize);
            const float line2W = TextWidthFont(gFontBold, line2, lineSize);
            DrawTextFont(gFontBold, line1, (static_cast<float>(screenWidth) - line1W) * 0.5f, y,
                         lineSize, kTextDark);
            y += 36.0f;
            DrawTextFont(gFontBold, line2, (static_cast<float>(screenWidth) - line2W) * 0.5f, y,
                         lineSize, kTextDark);
            y += 40.0f;
            const std::string outcome = app.game.Outcome();
            const float outcomeW = TextWidthFont(gFontBold, outcome, lineSize);
            DrawTextFont(gFontBold, outcome, (static_cast<float>(screenWidth) - outcomeW) * 0.5f, y,
                         lineSize, kTextDark);
            y += 56.0f;
            const float btnFont = 18.0f;
            const std::string again = "Play again (same loop)";
            const std::string fresh = "Play again (new loop)";
            const std::string back = "Sandbox";
            const float againW = TextWidthFont(gFontBold, again, btnFont) + 28.0f;
            const float freshW = TextWidthFont(gFontBold, fresh, btnFont) + 28.0f;
            const float backW = TextWidthFont(gFontBold, back, btnFont) + 28.0f;
            const float gap = 16.0f;
            float bx = (static_cast<float>(screenWidth) - againW - freshW - backW - 2.0f * gap) * 0.5f;
            const Rectangle againR = {bx, y, againW, 36.0f};
            const Rectangle freshR = {bx + againW + gap, y, freshW, 36.0f};
            const Rectangle backR = {freshR.x + freshW + gap, y, backW, 36.0f};
            const Color playBase = {46, 158, 78, 255};
            const auto endButton = [&](Rectangle r, const std::string& label, Color base, bool bold) {
                const bool hover = CheckCollisionPointRec(mouse, r);
                DrawRectangleRounded(r, 0.35f, 6, hover ? ColorBrightness(base, -0.15f) : base);
                const float tw = bold ? TextWidthFont(gFontBold, label, btnFont) : TextWidth(label, btnFont);
                const float tx = r.x + (r.width - tw) * 0.5f;
                if (bold) {
                    DrawTextFont(gFontBold, label, tx, r.y + 8.0f, btnFont, WHITE);
                }
                else {
                    DrawTextS(label, tx, r.y + 8.0f, btnFont, WHITE);
                }
                return click && hover;
            };
            if (endButton(againR, again, playBase, true)) {
                app.ReplayUnpinningMatch();
            }
            if (endButton(freshR, fresh, playBase, true)) {
                app.StartRandomOddMatch();
            }
            if (endButton(backR, back, kAccent, true)) {
                app.BackToSandboxFromGame();
            }
        }

        if (!inGame && app.debug) {
            FramePhysicsTuning(app.tuningUI, app.physicsParams, gFont, mouse, click,
                               screenWidth, screenHeight);
        }

        EndDrawing();
    }
}

void RunApp()
{
    // Render at the monitor scale and multisample curves, on macOS and Windows.
    SetConfigFlags(FLAG_WINDOW_RESIZABLE | FLAG_WINDOW_HIGHDPI | FLAG_MSAA_4X_HINT);
    InitWindow(1460, 850, "The Unpinning Game");
    OpenWindowMaximized();
    SetExitKey(KEY_NULL); // ESC is used to unfocus the textbox
    SetTargetFPS(60);
    LoadUiFont();

    App app;
    while (!WindowShouldClose()) {
        Frame(app);
    }

    UnloadUiFont();
    CloseWindow();
}
