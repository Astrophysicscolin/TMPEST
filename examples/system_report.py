import os
import platform
import sys

import numpy as np

import tmpest as tp


def main():
    print("TMPEST SYSTEM REPORT")
    print(f"os: {platform.system()} {platform.release()} ({platform.platform()})")
    print(f"machine: {platform.machine()}, logical cpus: {os.cpu_count()}")
    print(f"python: {sys.version.split()[0]} ({platform.python_implementation()})")
    print(f"numpy: {np.__version__}")
    print(f"tmpest: {tp.__version__}")
    print(f"rust extension built: {tp.rust_available()}")
    print(f"active backend: {tp.backend()}, threads: {tp.num_threads()}")


if __name__ == "__main__":
    main()
