# Testing TMPEST on Windows

Thank you for testing. This takes about 15 minutes, most of it waiting for installs. Run the
lines one at a time, in order, in **PowerShell** (not Command Prompt, and not Git Bash). Each step says what you should see. If a step
fails, stop there and send back the full output of that step.

TMPEST is a Python library with a small compiled core written in Rust. To install it from
source you need three things on your machine: git, Python 3.10 or newer, and the Rust compiler.
Nothing here needs administrator rights except where a step says so, and everything is
installed inside one folder that you can delete afterwards.

> Status of this guide: written on 2026-11-15 and **not yet run on Windows** by the author.
> It was written from the Linux build, which works; Windows is the platform most likely to need a fix. Finding where it breaks is the point of the test.

These are PowerShell lines, not bash: Windows does not ship bash. If you use WSL, follow
`LINUX.md` inside WSL instead. Windows 10 or 11, 64-bit.

## 1. Install the tools

`winget` is the installer built into Windows 10 and 11. Each line may ask for permission.

```powershell
winget install --id Git.Git -e
```

```powershell
winget install --id Python.Python.3.12 -e
```

```powershell
winget install --id Rustlang.Rustup -e
```

Rust on Windows also needs Microsoft's C++ build tools (a large download, several GB):

```powershell
winget install --id Microsoft.VisualStudio.2022.BuildTools -e --override "--wait --passive --add Microsoft.VisualStudio.Workload.VCTools --includeRecommended"
```

**Close PowerShell and open a new window** so the new programs are found, then check:

```powershell
git --version
```

```powershell
python --version
```

```powershell
rustc --version
```

Expected: git prints a version, Python is 3.10 or newer, `rustc` is 1.83 or newer. If
`python` opens the Microsoft Store instead, use `py` in place of `python` in step 3 only.

## 2. Get the code

The repository is private, so a
GitHub sign-in window will open.

```powershell
git clone https://github.com/Astrophysicscolin/tmpest.git TMPEST
```

```powershell
cd TMPEST
```

## 3. Create and activate a virtual environment

A virtual environment is a private Python install inside the project folder.

```powershell
python -m venv .venv
```

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
```

```powershell
.venv\Scripts\Activate.ps1
```

Your prompt should now start with `(.venv)`. The `Set-ExecutionPolicy` line only affects this
PowerShell window; it lets the activation script run.

## 4. Build and install TMPEST

```powershell
python -m pip install --upgrade pip
```

```powershell
python -m pip install -e ".[dev,examples]"
```

The second line compiles the Rust core; the first time takes a few minutes. Expected last
line: `Successfully installed tmpest-1.4.0` (other packages may be listed too). An error
mentioning `link.exe` means the C++ build tools from step 1 are missing or PowerShell was not
reopened.

## 5. Check that the compiled core is in use

```powershell
python examples\system_report.py
```

Expected: `rust extension built: True` and `active backend: rust`.

## 6. Run the tests with the Rust core

```powershell
python -m pytest -q
```

Expected: `43 passed` in a few seconds.

## 7. Run the tests again with the pure-numpy fallback

```powershell
$env:TMPEST_BACKEND = "numpy"
```

```powershell
python -m pytest -q
```

```powershell
Remove-Item Env:TMPEST_BACKEND
```

Expected: `43 passed` again. The two backends are required to give identical numbers.

## 8. Run the demo

```powershell
python examples\demo.py
```

Expected: a table, then four lines starting `[OK]`. On the last three, the two numbers on each
line should agree within about 30 percent.

```powershell
python examples\gallery.py
```

Expected: three lines starting `[SAVE]`, and three PNG files in a new `figures` folder.

## 9. What to send back

Copy the output of these two commands, plus the `passed` / `failed` lines from steps 6 and 7:

```powershell
python examples\system_report.py
python examples\benchmark.py
```

The benchmark takes one to three minutes and prints a table comparing the plain-numpy code
with the Rust code on your machine. Every row should show a speed-up of 1.0x or more; the
"running median" and "1000000 x 9" rows should be tens of times faster.

## Cleaning up

```powershell
deactivate
```

Delete the `TMPEST` folder to remove everything from this test. Rust, Python, git and the
build tools can be removed from Settings, Apps.
