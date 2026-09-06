#!/usr/bin/env python3
"""Compile, and optionally run, complete Rust listings in design.tex."""

import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--offline", action="store_true", help="use cached dependencies")
    parser.add_argument("--run", action="store_true", help="also run each bounded example")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    source = (root / "docs/design.tex").read_text()
    examples = re.findall(
        r"^% example: ([a-z_]+)\n\\begin\{lstlisting\}[^\n]*\n(.*?)"
        r"\\end\{lstlisting\}",
        source,
        re.MULTILINE | re.DOTALL,
    )
    declared = re.findall(r"^% example: (.+)$", source, re.MULTILINE)
    names = [name for name, _ in examples]
    if not examples or names != declared or len(set(names)) != len(names):
        raise SystemExit("Missing, malformed, or duplicate example markers")

    with tempfile.TemporaryDirectory(prefix="csta-design-examples-") as temporary:
        project = Path(temporary)
        (project / "src/bin").mkdir(parents=True)
        manifest = [
            '[package]\nname = "csta-design-examples"\nversion = "0.0.0"',
            'edition = "2024"\nrust-version = "1.98.0"\n[workspace]',
            '[features]\ncheckpoint = ["csta/checkpoint"]',
            '[dependencies]\nrand = "=0.10.2"',
            f'csta = {{ path = {json.dumps(str(root / "csta"))} }}',
        ]
        for name, code in examples:
            (project / f"src/bin/{name}.rs").write_text(code)
            if name == "checkpoint":
                manifest.append(
                    '[[bin]]\nname = "checkpoint"\npath = "src/bin/checkpoint.rs"\n'
                    'required-features = ["checkpoint"]'
                )
        (project / "Cargo.toml").write_text("\n".join(manifest) + "\n")
        shutil.copyfile(root / "Cargo.lock", project / "Cargo.lock")
        common = [
            "--manifest-path", str(project / "Cargo.toml"),
            "--target-dir", str(root / "target/design-examples"),
        ]
        if args.offline:
            common.append("--offline")
        for features in ([], ["--features", "checkpoint"]):
            subprocess.run(
                ["cargo", "check", "--all-targets", *common, *features],
                check=True, timeout=300,
            )
        if args.run:
            for name, _ in examples:
                print(f"Running {name}", flush=True)
                features = ["--features", "checkpoint"] if name == "checkpoint" else []
                subprocess.run(
                    ["cargo", "run", "--quiet", *common, "--bin", name, *features],
                    check=True, timeout=120,
                )
    action = "compiled and ran" if args.run else "compiled"
    print(f"Successfully {action} {len(examples)} document examples.")


if __name__ == "__main__":
    main()
