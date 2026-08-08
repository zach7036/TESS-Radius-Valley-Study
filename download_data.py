#!/usr/bin/env python3
"""Download the public MAST TOI catalog used by the TESS radius-valley study."""

from __future__ import annotations

import hashlib
from pathlib import Path
from urllib.request import Request, urlopen

URL = "https://archive.stsci.edu/missions/tess/catalogs/toi/tois.csv"
OUTPUT = Path(__file__).resolve().parent / "tois_mast.csv"
STUDY_SHA256 = "6797d7ea49dedd95cf5a5711fe4c93553f218268d3ce68f9fbde569d8ff06a7e"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    req = Request(URL, headers={"User-Agent": "TESS-Radius-Valley-Study/1.0"})
    with urlopen(req, timeout=60) as response:
        data = response.read()

    digest = sha256_bytes(data)
    OUTPUT.write_bytes(data)
    print(f"Downloaded {len(data):,} bytes to {OUTPUT.name}")
    print(f"SHA-256: {digest}")

    if digest == STUDY_SHA256:
        print("Checksum matches the exact catalog snapshot used in the August 8, 2026 study.")
    else:
        print("WARNING: checksum differs from the study snapshot.")
        print("The live MAST catalog has changed; this run should be treated as an updated replication.")


if __name__ == "__main__":
    main()
