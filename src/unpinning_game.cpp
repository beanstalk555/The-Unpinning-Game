#include "unpinning_game.h"

#include <ctime>
#include <fstream>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

namespace
{
    const char* kInfinityUtf8 = "\xE2\x88\x9E";

    void EnsureLogDirectory()
    {
#ifdef _WIN32
        _mkdir("unpinning_logs");
#else
        mkdir("unpinning_logs", 0755);
#endif
    }

    bool FileExists(const std::string& path)
    {
        std::ifstream in(path.c_str());
        return static_cast<bool>(in);
    }
}

std::string UnpinningGame::AllowanceLabel(int index)
{
    if (index >= 5 || index < 0) {
        return kInfinityUtf8;
    }
    return std::to_string(index + 1);
}

void UnpinningGame::SetAllowanceIndex(int index)
{
    if (index < 0) {
        index = 0;
    }
    if (index >= kAllowanceCount) {
        index = kAllowanceCount - 1;
    }
    allowanceIndex_ = index;
}

void UnpinningGame::ToggleAllowanceMenu()
{
    allowanceOpen_ = !allowanceOpen_;
}

void UnpinningGame::CloseAllowanceMenu()
{
    allowanceOpen_ = false;
}

UnpinningGame::Phase UnpinningGame::Live() const
{
    return phase_ == Phase::ExitPrompt ? returnPhase_ : phase_;
}

UnpinningGame::Phase& UnpinningGame::Live()
{
    return phase_ == Phase::ExitPrompt ? returnPhase_ : phase_;
}

bool UnpinningGame::AcceptsMove() const
{
    if (Live() != Phase::Play || pending_) {
        return false;
    }
    if (unlimited_) {
        return true;
    }
    return movesUsed_ < allowance_;
}

bool UnpinningGame::AcceptsPinRemoval() const
{
    return Live() == Phase::Play && !pending_;
}

bool UnpinningGame::NeedsSettle() const
{
    if (phase_ == Phase::Off || phase_ == Phase::Over) {
        return false;
    }
    const Phase live = Live();
    if (live == Phase::Tangle) {
        return resetAwaiting_;
    }
    return pending_ && rewritten_;
}

const char* UnpinningGame::MoveName(MoveKind kind)
{
    switch (kind) {
    case MoveKind::R0:
        return "R0";
    case MoveKind::R1:
        return "R1";
    case MoveKind::R2:
        return "R2";
    case MoveKind::R3:
        return "R3";
    }
    return "R0";
}

int UnpinningGame::PointsFor(MoveKind kind)
{
    if (kind == MoveKind::R1) {
        return 1;
    }
    if (kind == MoveKind::R2) {
        return 2;
    }
    return 0;
}

std::string UnpinningGame::PointsPhrase(int points)
{
    if (points == 1) {
        return "1 point";
    }
    return std::to_string(points) + " points";
}

void UnpinningGame::Append(const std::string& line)
{
    log_ += line;
    log_ += '\n';
    Flush();
}

void UnpinningGame::Flush() const
{
    if (logPath_.empty()) {
        return;
    }
    std::ofstream out(logPath_.c_str(), std::ios::binary | std::ios::trunc);
    if (!out) {
        return;
    }
    out << log_;
}

void UnpinningGame::OpenLog()
{
    EnsureLogDirectory();
    const std::time_t now = std::time(nullptr);
    const std::tm local = *std::localtime(&now);
    char stamp[32];
    std::strftime(stamp, sizeof(stamp), "%Y%m%d-%H%M%S", &local);
    std::string path = std::string("unpinning_logs/game-") + stamp + ".txt";
    int suffix = 2;
    while (FileExists(path)) {
        path = std::string("unpinning_logs/game-") + stamp + "-" + std::to_string(suffix) + ".txt";
        ++suffix;
    }
    logPath_ = path;
    log_.clear();
}

void UnpinningGame::Begin(const std::string& sigma, const std::string& infiniteRegion,
                          const std::string& pins)
{
    unlimited_ = allowanceIndex_ >= 5;
    allowance_ = unlimited_ ? 0 : allowanceIndex_ + 1;
    phase_ = Phase::Play;
    returnPhase_ = Phase::Play;
    player_ = 1;
    scores_[0] = 0;
    scores_[1] = 0;
    movesUsed_ = 0;
    pending_ = false;
    rewritten_ = false;
    pointsApplied_ = false;
    resetAwaiting_ = false;
    flash_.clear();
    flashUntil_ = 0.0;
    rulesOpen_ = false;
    allowanceOpen_ = false;
    OpenLog();
    Append(std::string("Sigma: ") + sigma);
    Append(std::string("Infinite region: ") + infiniteRegion);
    Append(std::string("Moves per turn: ") + (unlimited_ ? "infinity" : std::to_string(allowance_)));
    Append(std::string("Starting pins: ") + (pins.empty() ? "none" : pins));
    Append("");
    Append("Player 1's turn");
}

void UnpinningGame::Deactivate()
{
    phase_ = Phase::Off;
    returnPhase_ = Phase::Play;
    pending_ = false;
    rewritten_ = false;
    resetAwaiting_ = false;
    flash_.clear();
    allowanceOpen_ = false;
    rulesOpen_ = false;
}

void UnpinningGame::Tick(double now)
{
    if (!flash_.empty() && now >= flashUntil_) {
        flash_.clear();
    }
}

void UnpinningGame::OnMoveStarted(MoveKind kind, const std::string& region)
{
    pending_ = true;
    rewritten_ = false;
    pointsApplied_ = false;
    pendingKind_ = kind;
    pendingPoints_ = PointsFor(kind);
    pendingRegion_ = region;
}

void UnpinningGame::ApplyScore(double now)
{
    if (!pending_ || pointsApplied_) {
        return;
    }
    scores_[player_ - 1] += pendingPoints_;
    pointsApplied_ = true;
    if (pendingPoints_ <= 0) {
        return;
    }
    const int points = pendingPoints_;
    const std::string noun = (points == 1) ? "point" : "points";
    const std::string intersections = (points == 1) ? "intersection" : "intersections";
    flash_ = "Player " + std::to_string(player_) + " scores " + std::to_string(points) + " " + noun +
             " by removing " + std::to_string(points) + " " + intersections + "!";
    flashUntil_ = now + 3.0;
}

void UnpinningGame::OnRewriteFinished(double now)
{
    if (!pending_) {
        return;
    }
    rewritten_ = true;
    ApplyScore(now);
}

void UnpinningGame::OnMoveAborted()
{
    pending_ = false;
    rewritten_ = false;
    pointsApplied_ = false;
}

std::string UnpinningGame::MoveLine(bool includePoints) const
{
    std::string line = std::string(MoveName(pendingKind_)) + " on " + pendingRegion_;
    if (includePoints) {
        line += ", " + PointsPhrase(pendingPoints_);
    }
    return line;
}

void UnpinningGame::CommitPending(bool anyPins, double now)
{
    (void)now;
    if (!pending_) {
        return;
    }
    ++movesUsed_;
    Append(MoveLine(true));
    pending_ = false;
    rewritten_ = false;
    pointsApplied_ = false;
    const bool exhausted = !unlimited_ && movesUsed_ >= allowance_;
    if (exhausted && !anyPins) {
        EndTurn();
    }
}

void UnpinningGame::EndTurn()
{
    player_ = (player_ == 1) ? 2 : 1;
    movesUsed_ = 0;
    Append("");
    Append("Player " + std::to_string(player_) + "'s turn");
}

void UnpinningGame::OnBoardCleared(double now)
{
    if (pending_) {
        ApplyScore(now);
        Append(MoveLine(true));
        pending_ = false;
        rewritten_ = false;
        pointsApplied_ = false;
        flash_.clear();
    }
    phase_ = Phase::Over;
    returnPhase_ = Phase::Over;
    rulesOpen_ = false;
    Append("");
    Append("Game over.");
    Append(ScoreLine());
    Append(Outcome());
}

void UnpinningGame::OnSettled(bool tangled, bool anyPins, double now)
{
    if (phase_ == Phase::Off || phase_ == Phase::Over) {
        return;
    }
    Phase& live = Live();
    if (live == Phase::Tangle) {
        if (tangled) {
            resetAwaiting_ = false;
            return;
        }
        Append(resetAwaiting_ ? "Cleared by reset and autotighten." : "Cleared by dragging.");
        resetAwaiting_ = false;
        live = Phase::Play;
        if (pending_) {
            CommitPending(anyPins, now);
        }
        return;
    }
    if (!pending_ || !rewritten_) {
        return;
    }
    if (tangled) {
        Append(MoveLine(false) + " tangled.");
        Append("Tangle detected.");
        live = Phase::Tangle;
        return;
    }
    CommitPending(anyPins, now);
}

void UnpinningGame::OnDragCleared(bool anyPins, double now)
{
    if (phase_ == Phase::Off || phase_ == Phase::Over) {
        return;
    }
    Phase& live = Live();
    if (live != Phase::Tangle) {
        return;
    }
    Append("Cleared by dragging.");
    resetAwaiting_ = false;
    live = Phase::Play;
    if (pending_) {
        CommitPending(anyPins, now);
    }
}

void UnpinningGame::OnTangleUndone()
{
    if (phase_ == Phase::Off || phase_ == Phase::Over) {
        return;
    }
    Phase& live = Live();
    if (live != Phase::Tangle) {
        return;
    }
    Append("Cleared by undo.");
    if (pointsApplied_) {
        scores_[player_ - 1] -= pendingPoints_;
        if (scores_[player_ - 1] < 0) {
            scores_[player_ - 1] = 0;
        }
    }
    pending_ = false;
    rewritten_ = false;
    pointsApplied_ = false;
    resetAwaiting_ = false;
    flash_.clear();
    live = Phase::Play;
}

void UnpinningGame::OnResetPressed()
{
    resetAwaiting_ = true;
}

void UnpinningGame::OnPinRemoved(const std::string& region)
{
    if (!AcceptsPinRemoval()) {
        return;
    }
    Append("Removed pin on " + region + ".");
    EndTurn();
}

void UnpinningGame::AskExit()
{
    if (phase_ == Phase::Play || phase_ == Phase::Tangle) {
        returnPhase_ = phase_;
        phase_ = Phase::ExitPrompt;
    }
}

void UnpinningGame::DismissExit()
{
    if (phase_ == Phase::ExitPrompt) {
        phase_ = returnPhase_;
    }
}

void UnpinningGame::MarkEarlyExit()
{
    if (phase_ == Phase::Off) {
        return;
    }
    Append("");
    Append("Early exit to sandbox.");
    Deactivate();
}

std::string UnpinningGame::TurnTitle() const
{
    return "Player " + std::to_string(player_) + "'s turn.";
}

std::string UnpinningGame::StatusDetail(bool anyPins) const
{
    if (Live() == Phase::Tangle) {
        return "Tangle detected.";
    }
    if (!flash_.empty()) {
        return flash_;
    }
    const bool noMoves = !unlimited_ && movesUsed_ >= allowance_;
    if (noMoves && anyPins) {
        return "No moves left: Remove a pin.";
    }
    std::string sentence;
    if (unlimited_) {
        sentence = "Make as many moves as you want";
    }
    else {
        const int remaining = allowance_ - movesUsed_;
        sentence = "Make at most " + std::to_string(remaining) + (remaining == 1 ? " move" : " moves");
    }
    if (anyPins) {
        sentence += " and then remove a pin.";
    }
    else {
        sentence += ".";
    }
    return sentence;
}

std::string UnpinningGame::ScoreLine() const
{
    return "Player 1: " + PointsPhrase(scores_[0]) + "     Player 2: " + PointsPhrase(scores_[1]);
}

int UnpinningGame::Score(int player) const
{
    if (player == 2) {
        return scores_[1];
    }
    return scores_[0];
}

std::string UnpinningGame::Outcome() const
{
    if (scores_[0] > scores_[1]) {
        return "Player 1 wins!";
    }
    if (scores_[1] > scores_[0]) {
        return "Player 2 wins!";
    }
    return "The game is a draw.";
}
