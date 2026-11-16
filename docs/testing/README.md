# Native test guides

One guide per operating system. Each is a list of single lines to run in order, with the
expected result after each step, ending with what to send back.

| System | Guide | Shell | Run by the author? |
|---|---|---|---|
| macOS (Apple Silicon and Intel) | `MACOS.md` | Terminal (zsh) | Run on Apple Silicon 2026-11-15, Python 3.12: 43 passed on each backend, demo and benchmark run; Intel not yet |
| Linux (Debian/Ubuntu; notes for Fedora, Arch, WSL) | `LINUX.md` | bash | Same steps passed in a Linux container, x86-64, Python 3.13 |
| Windows 10/11 | `WINDOWS.md` | PowerShell | Not yet |

Before sending these out: add each tester as a collaborator on the private repository so `git clone` works for them.

Send results back to the author.

The package is `tmpest`, without an E after the T. `pip install tempest` installs an
unrelated package.
