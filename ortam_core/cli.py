"""ORTAM command-line interface.

Usage:
  ortam list                               List available modules
  ortam run <module> <input.yaml>          Run a module standalone
  ortam run <module> <input.yaml> --json   Write the result as JSON
"""

import argparse
import sys
from pathlib import Path

import yaml

from modules.registry import MODULES, get_module


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="ortam", description="ORTAM environmental modelling platform")
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("list", help="List available modules")
    r = sub.add_parser("run", help="Run a module standalone")
    r.add_argument("module")
    r.add_argument("input", type=Path)
    r.add_argument("--json", action="store_true", help="Write the result as JSON")
    a = p.parse_args(argv)

    if a.command == "list":
        for name in MODULES:
            print(f"{name:30s} {get_module(name).title}")
        return 0

    module = get_module(a.module)
    inputs = module.validate(yaml.safe_load(a.input.read_text(encoding="utf-8")))
    output = module.run(inputs)
    print(output.model_dump_json(indent=2) if a.json else module.summary(output))
    return 0


if __name__ == "__main__":
    sys.exit(main())
