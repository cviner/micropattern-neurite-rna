#!/usr/bin/env python3
"""Download the Middleton et al. (2019) supplementary tables used for overlap tests."""

import hashlib
import urllib.request

from common import EXTERNAL

BASE = "https://static-content.springer.com/esm/art%3A10.1186%2Fs12915-019-0630-z/MediaObjects/"
FILES = {
    "Middleton2019_MOESM1.xlsx": ("12915_2019_630_MOESM1_ESM.xlsx",
                                  "b074b4199e0f6b6a8a485d84f099375c2f04d4327f5bbc549f5784759ece648b"),
    "Middleton2019_MOESM5.xlsx": ("12915_2019_630_MOESM5_ESM.xlsx",
                                  "c8a05c603473866e6e051c4d6dc621a46176704ee5faa344f4e5f8887716c6d2"),
    "Middleton2019_MOESM6.xlsx": ("12915_2019_630_MOESM6_ESM.xlsx",
                                  "bdf5e0422f9a2a59836596d5c3dd98e88a47a7eb8529a0b5e4f2ae27fc975e70"),
}


def main() -> None:
    EXTERNAL.mkdir(parents=True, exist_ok=True)
    for local, (remote, sha) in FILES.items():
        path = EXTERNAL / local
        if not path.exists():
            urllib.request.urlretrieve(BASE + remote, path)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != sha:
            raise SystemExit(f"{local}: checksum {digest} does not match {sha}")
        print(f"{local}\tok")


if __name__ == "__main__":
    main()
