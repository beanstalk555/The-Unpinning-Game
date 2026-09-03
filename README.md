# The Unpinning Game

A raylib prototype for drawing and tightening multiloop diagrams.

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
