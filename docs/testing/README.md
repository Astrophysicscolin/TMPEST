# Native test guides

One guide per operating system. Each is a list of single lines to run in order, with the
expected result after each step, ending with what to send back.

| System | Guide | Shell | Run by the author? |
|---|---|---|---|
| macOS (Apple Silicon and Intel) | `MACOS.md` | Terminal (zsh) | Run on Apple Silicon 2026-11-15, Python 3.12: 43 passed on each backend, demo and benchmark run; Intel not yet |
| Linux (Debian/Ubuntu; notes for Fedora, Arch, WSL) | `LINUX.md` | bash | Guide not yet followed by a tester. The automated tests pass on GitHub's Ubuntu runner, Python 3.10 and 3.14, both backends (2026-11-16) |
| Windows 10/11 | `WINDOWS.md` | PowerShell | Guide not yet followed by a tester. The automated tests pass on GitHub's Windows runner, Python 3.10 and 3.14, both backends (2026-11-16) |

Testers get the code as a source package, `tmpest-1.5.0.tar.gz`, built with `maturin sdist`; it holds the library, the tests, the example scripts these guides use and the small public CSVs in `data/`. Step 3 (step 2 on Windows) unpacks it.

Send results back to the author.

The package is `tmpest`, without an E after the T. `pip install tempest` installs an
unrelated package.
