#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


NOTES = {"INPUT_MANIFEST.tsv", "SHA256SUMS.txt", "GEO_SOURCE.tsv"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_manifest(path: Path) -> list[dict[str, str]]:
    rows = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("short_name"):
                continue
            fields = line.split("\t")
            if len(fields) == 3:
                name, size, digest = fields
            elif len(fields) == 4:
                name, _geo_filename, size, digest = fields
            else:
                raise ValueError(f"{path} has a row with {len(fields)} columns: {line}")
            rows.append({"short_name": name, "bytes": size, "sha256": digest})
    if not rows:
        raise ValueError(f"{path} has no files")
    return rows


def check(root: Path, samples_csv: Path) -> str:
    root = root.resolve()
    inputs = root / "inputs"
    manifest_path = inputs / "GEO_SOURCE.tsv"
    if not manifest_path.is_file():
        raise ValueError(
            f"missing {manifest_path}. Offline mode checks this table of the "
            "matrices the workflow reads. GSE274934_suppl/manifest.tsv is not required."
        )
    listed = read_manifest(manifest_path)
    extra_manifest = inputs / "INPUT_MANIFEST.tsv"
    if extra_manifest.is_file():
        extra = {row["short_name"]: row for row in read_manifest(extra_manifest)}
        for row in listed:
            other = extra.get(row["short_name"])
            if other is None or other["bytes"] != row["bytes"] or other["sha256"] != row["sha256"]:
                raise ValueError(f"INPUT_MANIFEST.tsv disagrees with GEO_SOURCE.tsv for {row['short_name']}")
    by_name = {row["short_name"]: row for row in listed}
    if len(by_name) != len(listed):
        raise ValueError(f"duplicate short_name in {manifest_path}")

    with samples_csv.open() as handle:
        samples = [row["sample_id"] for row in csv.DictReader(handle)]
    if len(samples) != len(set(samples)):
        raise ValueError(f"duplicate sample_id in {samples_csv}")

    required = []
    for sample in samples:
        required.extend([
            f"{sample}_filtered.h5",
            f"{sample}_raw.h5",
            f"{sample}_tcr.csv.gz",
        ])
        peaks = f"{sample}_atac_peaks.h5"
        single = f"{sample}_atac_singlecell.csv.gz"
        have = peaks in by_name, single in by_name
        if have == (True, False) or have == (False, True):
            raise ValueError(f"{sample} has only one of the two ATAC files in {manifest_path}")
    missing = [name for name in required if name not in by_name]
    if missing:
        raise ValueError("cohort sample sheet files absent from GEO_SOURCE.tsv: " + ", ".join(missing))

    filtered = {name[: -len("_filtered.h5")] for name in by_name if name.endswith("_filtered.h5")}
    if filtered != set(samples):
        raise ValueError(
            "filtered matrices and the sample sheet name different donors: "
            f"manifest {sorted(filtered)} sheet {samples}"
        )

    present = {path.name for path in inputs.iterdir() if path.name not in NOTES}
    if present != set(by_name):
        extra = sorted(present - set(by_name))
        absent = sorted(set(by_name) - present)
        raise ValueError(f"inputs/ does not match GEO_SOURCE.tsv. extra={extra} absent={absent}")

    for name, row in sorted(by_name.items()):
        path = inputs / name
        size = path.stat().st_size
        if size != int(row["bytes"]):
            raise ValueError(f"{name} is {size} bytes, manifest says {row['bytes']}")
        digest = sha256(path)
        if digest != row["sha256"]:
            raise ValueError(f"{name} sha256 {digest} does not match GEO_SOURCE.tsv")

    sums = inputs / "SHA256SUMS.txt"
    if sums.is_file():
        for line in sums.read_text().splitlines():
            if not line.strip():
                continue
            digest, name = line.split()
            row = by_name.get(name)
            if row is None or row["sha256"] != digest:
                raise ValueError(f"SHA256SUMS.txt disagrees with GEO_SOURCE.tsv for {name}")

    supplement = root / "GSE274934_suppl" / "manifest.tsv"
    supplement_note = "GSE274934_suppl/manifest.tsv absent; full GEO download was not checked."
    if supplement.is_file():
        bad = []
        n = 0
        for line in supplement.read_text().splitlines():
            if not line.strip():
                continue
            name, size = line.split("\t")
            path = root / "GSE274934_suppl" / name
            got = path.stat().st_size if path.exists() else -1
            n += 1
            if got != int(size):
                bad.append(f"{name} {got} != {size}")
        if bad:
            raise ValueError("GEO size check failed:\n" + "\n".join(bad))
        supplement_note = f"GEO supplement size check OK ({n} files)."
    return f"Workflow input checksums OK ({len(by_name)} files). {supplement_note}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--samples", required=True, type=Path)
    args = parser.parse_args()
    try:
        print(">> " + check(args.root, args.samples))
    except ValueError as exc:
        raise SystemExit(f"ERROR: {exc}")


if __name__ == "__main__":
    main()
