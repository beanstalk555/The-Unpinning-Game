# The Unpinning Game

A raylib prototype for the unpinning game. Features random multiloop generation and drawing, a custom physics engine for tightening and settling, the ability to add and remove pins, and dynamic performance of Reidemeister moves. The combinatorics of the multiloop are stored and manipulated using the permutation representation.

## Windows

### Requirements

- [raylib](https://www.raylib.com/) with the w64devkit toolchain
- Default install layout: `C:\raylib\raylib` and `C:\raylib\w64devkit`

### Download

```sh
git clone <your-repo-url>
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

### AI Disclosure

The codebase for this prototype was mostly AI-generated using Cursor (Pro plan). The models used include:

- Composer 2
- Composer 2.5
- Cursor Grok 4.6
- Claude Fable 5
- Claude Fable 5.1