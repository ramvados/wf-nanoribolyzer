#!/usr/bin/env python3
"""Retain raw reads absent from Barbell's trimmed output."""
import argparse
import csv
import gzip
from pathlib import Path

def records(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if (not header.startswith("@") or not plus.startswith("+")
                    or not quality
                    or len(sequence.rstrip("\r\n")) != len(quality.rstrip("\r\n"))):
                raise SystemExit("Invalid FASTQ record in {}".format(path))
            yield header[1:].split()[0], header + sequence + plus + quality

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--raw", type=Path, required=True)
parser.add_argument("--trimmed", type=Path, required=True)
parser.add_argument("--filtered", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--counts", type=Path, required=True)
parser.add_argument("--trim-context", action="store_true")
parser.add_argument("--barcode-fasta", type=Path)
args = parser.parse_args()

prefixes = {}
if args.trim_context:
    if args.barcode_fasta is None:
        parser.error("--trim-context requires --barcode-fasta")
    queries = {}
    label = None
    with args.barcode_fasta.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                label = line[1:].split()[0]
                if label in queries:
                    raise SystemExit("Duplicate barcode label: " + label)
                queries[label] = ""
            elif label is None:
                raise SystemExit("Invalid barcode FASTA.")
            else:
                queries[label] += line.upper()
    if len(queries) < 2:
        raise SystemExit("At least two barcode queries are required.")
    for label, sequence in queries.items():
        if not sequence.startswith("GGCTTCTT") or not sequence.endswith("CCAATTTGTGGGTTCGT"):
            raise SystemExit("Context trimming requires DRB004 flanks: " + label)
        core = sequence[8:-17]
        if len(core) < 18 or any(base not in "ACGT" for base in core):
            raise SystemExit("Unsupported barcode core: " + label)
        prefixes.setdefault(core[:18], set()).add(label)

def context_trim(record):
    if not args.trim_context:
        return record, 0
    lines = record.splitlines()
    sequence = lines[1].upper()
    hits = []
    position = sequence.find("GGCTTCTT")
    while position >= 0:
        prefix = sequence[position + 8:position + 26]
        for label in prefixes.get(prefix, ()):
            hits.append(position)
        position = sequence.find("GGCTTCTT", position + 1)
    if len(hits) != 1:
        return record, 0
    start = hits[0]
    if (start < 10 or len(sequence) - start > 250
            or sequence[start - 10:start] != "AAAAAAAAAA"):
        return record, 0
    removed = len(sequence) - start
    lines[1] = lines[1][:start]
    lines[3] = lines[3][:start]
    return "\n".join(lines) + "\n", removed

with args.filtered.open() as handle:
    accepted = {row["read_id"] for row in csv.DictReader(handle, delimiter="\t")}

trimmed = set()
for path in sorted(list(args.trimmed.glob("*.fastq")) + list(args.trimmed.glob("*.fastq.gz"))):
    for read_id, record in records(path):
        if read_id in trimmed:
            raise SystemExit("Duplicate trimmed read: " + read_id)
        trimmed.add(read_id)

if not trimmed <= accepted:
    raise SystemExit("Trimmed reads absent from filter table.")

if args.output.exists() or args.counts.exists():
    raise SystemExit("Output already exists; nothing overwritten.")

seen = set()
fallback = 0
context_trimmed = 0
try:
    with gzip.open(args.output, "xt") as output:
        for read_id, record in records(args.raw):
            if read_id in seen:
                raise SystemExit("Duplicate raw read: " + read_id)
            seen.add(read_id)
            if read_id not in trimmed:
                record, removed = context_trim(record)
                output.write(record)
                fallback += 1
                context_trimmed += int(removed > 0)
    if not accepted <= seen:
        raise SystemExit("Filter IDs absent from raw input.")
    with args.counts.open("x") as output:
        output.write("category\treads\n")
        for label, number in (
            ("raw_input", len(seen)),
            ("barbell_trimmed", len(trimmed)),
            ("accepted_without_trimmed_output", len(accepted - trimmed)),
            ("not_retained_by_filter", len(seen - accepted)),
            ("fallback_raw", fallback),
            ("fallback_context_trimmed", context_trimmed),
            ("fallback_unmodified", fallback - context_trimmed),
        ):
            output.write("{}\t{}\n".format(label, number))
except BaseException:
    args.output.unlink(missing_ok=True)
    raise
