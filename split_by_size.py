#!/usr/bin/env python3
"""
split_by_size.py
----------------
Split a du-format file list into batches of up to MAX_TB terabytes each.

Input format (one entry per line, tab-separated):
    <size_in_1K_blocks>  <file_path>

Output:
    For an input file named  bam.list.cu01.txt
    produces               bam.list.cu01.batch01.txt
                           bam.list.cu01.batch02.txt  ...

Usage:
    python3 split_by_size.py bam.list.cu01.txt [--max-tb 10] [--output-dir .]

    # Process multiple files at once:
    python3 split_by_size.py data_paths/bam.list.*.txt --max-tb 10
"""

import argparse
import os
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
BYTES_PER_KB   = 1_024                          # du reports 1 K-blocks
TB_IN_BYTES    = 1_000_000_000_000              # 1 TB  (SI, storage standard)
KB_PER_TB      = TB_IN_BYTES / BYTES_PER_KB     # 976_562_500  KB per TB


def parse_entries(input_path: Path) -> tuple[list[tuple[int, str]], int]:
    """Read and parse a du-format list file.

    Returns (entries, skipped_count) where each entry is (size_kb, original_line).
    """
    entries: list[tuple[int, str]] = []
    skipped = 0
    with open(input_path) as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            if len(parts) != 2:
                print(f"  WARNING line {lineno}: cannot parse '{line}' — skipped",
                      file=sys.stderr)
                skipped += 1
                continue
            try:
                size_kb = int(parts[0])
            except ValueError:
                print(f"  WARNING line {lineno}: non-integer size '{parts[0]}' — skipped",
                      file=sys.stderr)
                skipped += 1
                continue
            entries.append((size_kb, line))
    return entries, skipped


def group_pairs(entries: list[tuple[int, str]]) -> list[list[tuple[int, str]]]:
    """Group entries so that a .bam and its .bam.bai are one atomic unit.

    Strategy: build a lookup of bam_path → entry, then for each .bam.bai entry
    attach it to the preceding .bam.  Any unpaired files form their own unit.
    """
    # Map bam_path → (size_kb, line) for quick lookup
    bam_index: dict[str, int] = {}   # bam_path -> position in units list
    units: list[list[tuple[int, str]]] = []

    for size_kb, line in entries:
        path = line.split(None, 1)[1].strip()

        if path.endswith(".bam.bai"):
            bam_path = path[: -len(".bai")]   # strip the .bai suffix
            if bam_path in bam_index:
                # Attach BAI to its BAM unit
                units[bam_index[bam_path]].append((size_kb, line))
                continue
            else:
                print(
                    f"  WARNING: BAI has no matching BAM in this file, treating standalone:\n"
                    f"    {path}",
                    file=sys.stderr,
                )

        # New unit (BAM or anything else)
        idx = len(units)
        units.append([(size_kb, line)])
        if path.endswith(".bam"):
            bam_index[path] = idx

    return units


def split_file(input_path: Path, max_tb: float, output_dir: Path) -> None:
    """Split *input_path* into batch files, each ≤ max_tb TB."""

    max_kb = max_tb * KB_PER_TB

    # ---- read & parse -------------------------------------------------------
    entries, skipped = parse_entries(input_path)

    if not entries:
        print(f"  No valid entries found in {input_path.name} — nothing written.")
        return

    # ---- group BAM + BAI pairs into atomic units ----------------------------
    units = group_pairs(entries)   # list of lists; each inner list goes together

    # ---- bin units into batches ---------------------------------------------
    batches: list[list[str]] = []
    current_batch: list[str] = []
    current_kb = 0

    for unit in units:
        unit_kb = sum(s for s, _ in unit)
        unit_lines = [line for _, line in unit]

        # Warn if a single unit exceeds the batch limit
        if unit_kb > max_kb:
            print(
                f"  WARNING: unit larger than {max_tb} TB "
                f"({unit_kb / KB_PER_TB:.2f} TB) — placed in its own batch:\n"
                + "\n".join(f"    {l.split(None,1)[1]}" for l in unit_lines),
                file=sys.stderr,
            )
            # Flush whatever is in the current batch first, then give the
            # oversized unit its own batch immediately so it never combines
            # with any other unit and inflates a batch beyond the limit.
            if current_batch:
                batches.append(current_batch)
            batches.append(unit_lines)
            current_batch = []
            current_kb = 0
            continue

        # Start a new batch when adding this unit would exceed the limit
        if current_batch and (current_kb + unit_kb) > max_kb:
            batches.append(current_batch)
            current_batch = []
            current_kb = 0

        current_batch.extend(unit_lines)
        current_kb += unit_kb

    if current_batch:
        batches.append(current_batch)

    # ---- write output files -------------------------------------------------
    stem = input_path.stem          # e.g. "bam.list.cu01"
    suffix = input_path.suffix      # e.g. ".txt"
    output_dir.mkdir(parents=True, exist_ok=True)

    pad = len(str(len(batches)))    # zero-padding width

    print(f"\n{input_path.name}")
    print(f"  Total entries : {len(entries):,}  ({skipped} skipped)")
    print(f"  Atomic units  : {len(units):,}  (BAM+BAI pairs kept together)")
    total_tb = sum(s for s, _ in entries) / KB_PER_TB
    print(f"  Total size    : {total_tb:.2f} TB")
    print(f"  Batches       : {len(batches)}  (limit {max_tb} TB each)")

    for i, batch in enumerate(batches, 1):
        batch_kb = sum(int(line.split(None, 1)[0]) for line in batch)
        out_name = f"{stem}.batch{i:0{pad}d}{suffix}"
        out_path = output_dir / out_name
        with open(out_path, "w") as fh:
            fh.write("\n".join(batch) + "\n")
        print(f"  [{i:0{pad}d}] {out_name}  —  {len(batch):,} files,  "
              f"{batch_kb / KB_PER_TB:.2f} TB  →  {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Split du-format file lists into ≤N TB batches for EGA upload."
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        metavar="LIST_FILE",
        help="One or more du-format list files to split.",
    )
    parser.add_argument(
        "--max-tb",
        type=float,
        default=10.0,
        metavar="TB",
        help="Maximum size per batch in terabytes (default: 10).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        metavar="DIR",
        help="Directory for output files (default: same directory as each input).",
    )
    args = parser.parse_args()

    for raw in args.inputs:
        path = Path(raw)
        if not path.is_file():
            print(f"ERROR: {raw} is not a file — skipped", file=sys.stderr)
            continue
        out_dir = args.output_dir if args.output_dir else path.parent
        split_file(path, args.max_tb, out_dir)


if __name__ == "__main__":
    main()
