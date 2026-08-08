#!/usr/bin/env python3
"""Launcher for the exact TESS radius-valley analysis source.

The original analysis file is preserved verbatim in ordered chunks under
analysis_parts/ so the complete source can be published through the connected
GitHub interface without changing its execution semantics. Each chunk is
compiled with its own filename and executed into the same global namespace.
The final chunk contains the original main() guard.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PARTS = sorted((ROOT / "analysis_parts").glob("part_*.py.inc"))
if not PARTS:
    raise SystemExit("No analysis source parts found under analysis_parts/")

for part in PARTS:
    source = part.read_text(encoding="utf-8")
    exec(compile(source, str(part), "exec"), globals(), globals())
