# Testing TMPEST on macOS

Thank you for testing. This takes about 15 minutes, most of it waiting for installs. Run the
lines one at a time, in order, in the Terminal app (zsh). Each step says what you should see. If a step
fails, stop there and send back the full output of that step.

TMPEST is a Python library with a small compiled core written in Rust. To install it from
source you need three things on your machine: git, Python 3.10 or newer, and the Rust compiler.
Nothing here needs administrator rights except where a step says so, and everything is
installed inside one folder that you can delete afterwards.

> Status of this guide: written on 2026-11-15 and **not yet run on macOS** by the author.
> It was written from the Linux build, which works. Finding where it breaks is the point of the test.

Works on both Apple Silicon (M1 and later) and Intel Macs.

## 1. Check what you already have

```bash
git --version
python3 --version
```

You need git, and Python 3.10 or newer. If `git` is missing, macOS offers to install the
Command Line Tools; accept, wait for it to finish, and run the line again. If it does not
offer, run:

```bash
xcode-select --install
```

If Python is older than 3.10 or missing, install it from https://www.python.org/downloads/ (or
`brew install python` if you use Homebrew), then open a new Terminal window.

## 2. Install the Rust compiler

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
```

```bash
source "$HOME/.cargo/env"
```

```bash
rustc --version
```

If the first line stops with `could not amend shell profile ... Permission denied`, the
installer could not edit one of your shell start-up files. Run it again with
`--no-modify-path` added, which skips that edit:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --no-modify-path
```

Then carry on with the `source` line. In that case you need to run the `source` line again
in every new terminal window before building.

Expected: `rustc 1.83` or newer. This installs into `~/.cargo` and `~/.rustup` in your home
folder only.

## 3. Get the code

The repository is private, so git will
ask you to sign in to GitHub.

```bash
git clone https://github.com/Astrophysicscolin/tmpest.git TMPEST
```

```bash
cd TMPEST
```

## 4. Create and activate a virtual environment

A virtual environment is a private Python install inside the project folder, so nothing here
touches the rest of your system.

```bash
python3 -m venv .venv
```

```bash
source .venv/bin/activate
```

Your prompt should now start with `(.venv)`. If you close the Terminal, come back with
`cd TMPEST` and the `source` line again.

## 5. Build and install TMPEST

```bash
python -m pip install --upgrade pip
```

```bash
python -m pip install -e ".[dev,examples]"
```

The second line compiles the Rust core; the first time takes one to three minutes. Expected
last line: `Successfully installed tmpest-1.4.0` (other packages may be listed too).

## 6. Check that the compiled core is in use

```bash
python examples/system_report.py
```

Expected: `rust extension built: True` and `active backend: rust`, with `threads` equal to
your number of CPU cores.

## 7. Run the tests with the Rust core

```bash
python -m pytest -q
```

Expected: `43 passed` in a few seconds.

## 8. Run the tests again with the pure-numpy fallback

```bash
TMPEST_BACKEND=numpy python -m pytest -q
```

Expected: `43 passed` again. The two backends are required to give identical numbers.

## 9. Run the demo

```bash
python examples/demo.py
```

Expected: a table, then four lines starting `[OK]`. On the last three, the two numbers on each
line should agree within about 30 percent.

```bash
python examples/gallery.py
```

Expected: three lines starting `[SAVE]`, and three PNG files in a new `figures` folder. Open
one and check that it shows plots, not blank panels.

## 10. What to send back

Copy the output of these two commands, plus the `passed` / `failed` lines from steps 7 and 8:

```bash
python examples/system_report.py
python examples/benchmark.py
```

The benchmark takes one to three minutes and prints a table comparing the plain-numpy code
with the Rust code on your machine. Every row should show a speed-up of 1.0x or more; the
"running median" and "1000000 x 9" rows should be tens of times faster.

## Cleaning up

```bash
deactivate
```

Delete the `TMPEST` folder to remove everything from this test. To remove Rust as well:
`rustup self uninstall`.
