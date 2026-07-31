#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import importlib
import shutil
import sys


PYTHON_PACKAGES = ("fastapi", "multipart", "uvicorn")
COMMANDS = ("node", "npm", "samtools", "bcftools", "minimap2")


def main():
    problems = []
    if sys.version_info[:2] != (3, 11):
        problems.append(
            f"Python 3.11 is required; found {sys.version.split()[0]}."
        )

    for package in PYTHON_PACKAGES:
        try:
            importlib.import_module(package)
        except ImportError:
            problems.append(f"Missing Python package: {package}")

    for command in COMMANDS:
        if not shutil.which(command):
            problems.append(f"Missing command: {command}")

    if problems:
        print("Environment check failed:")
        for problem in problems:
            print(f"  - {problem}")
        return 1

    print("Environment check passed.")
    print(f"Python: {sys.version.split()[0]}")
    for command in COMMANDS:
        print(f"{command}: {shutil.which(command)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
