# The Unpinning Game

A desktop implementation of the unpinning game.

## Credits and acknowledgements

This project grew out of a collaboration between Araf Israk, Christopher-Lloyd Simon, and Ben Stucky; all have contributed ideas to this work. In particular the implementation of the dynamic Reidemeister moves shares ideas from an earlier implementation of a "loop sandbox" coded by Araf in Summer 2026. We thank Lucas Burson for his help as well.

## Rules

The unpinning game is a 2-player game played on a loop in the plane which intersects itself, with pins in numerically-labeled complementary regions.

Your turn consists of two phases:

1. Click on the numerically-labeled regions without pins, attempting to decrease the number of self-intersections. You are awarded one point for each self-intersection that you remove. You can perform moves only on regions with zero, one, two, or three adjacent intersections.

2. Remove a pin, if there is at least one remaining.

The winner is the player with more points at the end.

### Technical notes

This version of the unpinning game is a dynamic variant of the "static" variant presented in this [YouTube lecture](https://www.youtube.com/watch?v=X9h2pMM0wrM) and [on this poster](https://christopherlloyd.github.io/LooPindex/unpinning_game/poster.svg). 

The moves that players can perform on their turn are the intersection-non-increasing Reidemeister moves (shortening R1, shortening R2, R3) which will be familiar to knot-theorists (though we don't care about over-under crossing information).

The game always ends with the loop becoming completely untangled. There is a very nice theorem of Hass and Scott that loops in surfaces in general position may be tightened using only these shortening moves. This is in stark contrast to the situation in knot theory where untangling diagrammatically may require introducing many extra crossings before a minimum number of crossings is achieved.

## App features

### No draws

One feature worth noting is that if the unpinning game is played on a loop with an odd number of crossings, it will never result in a draw.

### Implementation features

The implementation features:

- A custom physics engine for tightening and settling
- Dynamic performance of Reidemeister moves. 
- The permutation representation as the primary data structure for storing and manipulating multiloops

It also features a sandbox with random multiloop generation and drawing, with dedicated options for one-strand and simple cases. Pins can be added and removed freely, the multiloop can be dragged with left click, and regions can be highlighted with right click. From here the game can be played with "unconventional" settings including on multiloops with more than one strand, or with pins in only some regions.

## Installation

The easiest way to install the game is to use a pre-built binary. The [latest release](https://github.com/beanstalk555/The-Unpinning-Game/releases) contains binaries for Windows, macOS, and Linux under Assets.

## Building from source

The instructions below are for the Windows platform. They should be fairly straightforward to adapt for any platform, but I have not tested on other platforms.

### Requirements

- [raylib](https://www.raylib.com/) with the w64devkit toolchain
- Default install layout: `C:\raylib\raylib` and `C:\raylib\w64devkit`

### Download

```sh
git clone https://github.com/beanstalk555/The-Unpinning-Game/tree/main
cd The-Unpinning-Game
```

Or download and unzip the repository.

### Build

Open a terminal in the project folder and run:

```sh
make
```

If `make` is not on your PATH, use:

```sh
C:\raylib\w64devkit\bin\mingw32-make.exe
```

If raylib is installed somewhere else:

```sh
make RAYLIB_PATH=C:/path/to/raylib
```

### Run

```sh
.\game.exe
```

Close any running copy of the app before rebuilding.

## AI Disclosure

The codebase for this program was mostly AI-generated using Cursor (Pro plan). The models used include:

- Composer 2
- Composer 2.5
- Cursor Grok 4.6
- Cursor Grok 4.7
- Claude Fable 5
- Claude Fable 5.1
