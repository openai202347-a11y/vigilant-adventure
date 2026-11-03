"""Verify CDC input files against the archived SHA-256 manifest."""
import hashlib
import json
from pathlib import Path
import sys
root = Path(__file__).resolve().parents[1]
manifest = json.loads((root / "03_data/file_hashes.json").read_text(encoding="utf-8-sig"))
failed = False
for item in manifest:
    path = root / item["Path"]
    if not path.is_file():
        print(f"MISSING {item['Path']}")
        failed = True
        continue
    with path.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256").hexdigest()
    ok = digest.lower() == item["Hash"].lower()
    print(f"{'OK' if ok else 'MISMATCH'} {item['Path']}")
    failed |= not ok
sys.exit(1 if failed else 0)
