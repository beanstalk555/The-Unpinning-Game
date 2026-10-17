#pragma once

#include <string>

// Two-player match on the diagram that is already open in the sandbox.
// App owns one of these, draws it, and routes clicks. Physics stays in App.
class UnpinningGame
{
public:
    enum class MoveKind
    {
        R0,
        R1,
        R2,
        R3
    };

    static constexpr int kAllowanceCount = 6; // 1, 2, 3, 4, 5, unlimited

    bool Active() const { return phase_ != Phase::Off; }
    bool Over() const { return phase_ == Phase::Over; }
    bool ExitPrompt() const { return phase_ == Phase::ExitPrompt; }
    // Dragging is allowed only while resolving a tangle, not under the exit prompt.
    bool AllowsDrag() const { return phase_ == Phase::Tangle; }
    bool InTangle() const { return phase_ == Phase::Tangle; }

    bool AcceptsMove() const;
    bool AcceptsPinRemoval() const;
    bool MovePending() const { return pending_; }
    // A rewrite has finished, or a tangle reset is waiting for autotighten to stop.
    bool NeedsSettle() const;

    int AllowanceIndex() const { return allowanceIndex_; }
    void SetAllowanceIndex(int index);
    bool AllowanceMenuOpen() const { return allowanceOpen_; }
    void ToggleAllowanceMenu();
    void CloseAllowanceMenu();
    static std::string AllowanceLabel(int index);

    // Opens a new log and starts player 1 on the current diagram.
    void Begin(const std::string& sigma, const std::string& infiniteRegion,
               const std::string& pins);
    void Deactivate();

    void Tick(double now);

    void OnMoveStarted(MoveKind kind, const std::string& region);
    void OnRewriteFinished(double now);
    void OnMoveAborted();
    // The last circle was removed. Logs the move and ends the match.
    void OnBoardCleared(double now);
    // A move or a tangle reset has gone idle. `tangled` is the fresh check.
    void OnSettled(bool tangled, bool anyPins, double now);
    // The player dragged the diagram out of a tangle.
    void OnDragCleared(bool anyPins, double now);
    void OnTangleUndone();
    void OnResetPressed();
    void OnPinRemoved(const std::string& region);

    void AskExit();
    void DismissExit();
    // Writes the early-exit line and leaves the match. The live diagram stays.
    void MarkEarlyExit();

    bool RulesOpen() const { return rulesOpen_; }
    void OpenRules() { rulesOpen_ = true; }
    void CloseRules() { rulesOpen_ = false; }

    std::string TurnTitle() const;
    int CurrentPlayer() const { return player_; }
    // Move instruction, score flash, or tangle notice. Both flags come from the live necklace.
    std::string StatusDetail(bool anyPins, bool movesAvailable) const;
    std::string ScoreLine() const;
    int Score(int player) const;
    std::string Outcome() const;

private:
    enum class Phase
    {
        Off,
        Play,
        Tangle,
        ExitPrompt,
        Over
    };

    Phase Live() const;
    Phase& Live();
    void ApplyScore(double now);
    void CommitPending(bool anyPins, double now);
    void EndTurn();
    void Append(const std::string& line);
    void Flush() const;
    void OpenLog();
    // Moves the current log into tangled_games, completed_games, or other_games.
    void PlaceLog(const char* folder);
    std::string MoveLine(bool includePoints) const;
    static const char* MoveName(MoveKind kind);
    static int PointsFor(MoveKind kind);
    static std::string PointsPhrase(int points);

    Phase phase_ = Phase::Off;
    Phase returnPhase_ = Phase::Play;
    int allowanceIndex_ = 5;
    bool allowanceOpen_ = false;
    bool unlimited_ = true;
    int allowance_ = 0;
    int player_ = 1;
    int scores_[2] = {0, 0};
    int movesUsed_ = 0;

    bool pending_ = false;
    bool rewritten_ = false;
    bool pointsApplied_ = false;
    bool resetAwaiting_ = false;
    MoveKind pendingKind_ = MoveKind::R0;
    int pendingPoints_ = 0;
    std::string pendingRegion_;

    std::string flash_;
    double flashUntil_ = 0.0;
    bool rulesOpen_ = false;

    std::string log_;
    std::string logPath_;
    bool sawTangle_ = false;
};
