import gzip
import os
import tempfile
import time
import numpy as np
from prody import (
    applyMutinfoCorr,
    buildMutinfoMatrix,
    buildSeqidMatrix,
    calcShannonEntropy,
    parseMSA,
)

FASTA = "FASTA"
STOCKHOLM = "Stockholm"
EXPECTED = {
    "PF08171": (208, 68),
    "PF00089": (20152, 221),
    "PF16957": (154, 548),
}


def print_perf(name, timing_ms):
    print("ProDy," + name + "," + str(timing_ms))


def assert_msa(msa, label, expected_key):
    exp_seq, exp_len = EXPECTED[expected_key]
    assert msa.numSequences() == exp_seq, f"{label}: expected {exp_seq} sequences, got {msa.numSequences()}"
    assert msa.numResidues() == exp_len, f"{label}: expected {exp_len} columns, got {msa.numResidues()}"


def trim_inserts(msa):
    arr = msa.getArray().astype("U1")
    lower_or_gap = np.char.islower(arr) | (arr == "-") | (arr == ".")
    keep = ~np.all(lower_or_gap, axis=0)
    trimmed = msa[:, keep]
    return trimmed


def parse_stockholm(path, expected_key):
    if path.endswith(".gz"):
        # ProDy's parser can choke on blank lines in compressed inputs; filter them.
        with gzip.open(path, "rt") as handle, tempfile.NamedTemporaryFile(
            mode="wb", suffix=".sth", delete=False
        ) as tmp:
            for line in handle:
                if line.strip() == "":
                    continue
                tmp.write(line.encode())
        try:
            msa = parseMSA(tmp.name, format=STOCKHOLM)
        finally:
            os.remove(tmp.name)
    else:
        msa = parseMSA(path, format=STOCKHOLM)
    trimmed = trim_inserts(msa)
    assert_msa(trimmed, f"{expected_key} Stockholm", expected_key)
    return trimmed


def parse_fasta(path, expected_key, **kwargs):
    msa = parseMSA(path, **kwargs)
    trimmed = trim_inserts(msa)
    assert_msa(trimmed, f"{expected_key} FASTA", expected_key)
    return trimmed


def bench_read(path, label, expected_key, **kwargs):
    start = time.perf_counter()
    if kwargs.get("format") == STOCKHOLM:
        _ = parse_stockholm(path, expected_key)
    else:
        _ = parse_fasta(path, expected_key, **kwargs)
    elapsed = time.perf_counter() - start
    print(f"[BENCH] ProDy Read {label}: {elapsed}")


msa = parse_fasta("../data/PF08171.fasta", "PF08171", format=FASTA)

tmin = float("inf")
for _ in range(5):
    start = time.perf_counter()
    parse_stockholm("../data/PF08171.sth", "PF08171")
    elapsed = time.perf_counter() - start
    if elapsed < tmin:
        tmin = elapsed

print_perf("Read Pfam Stockholm MSA", 1000 * tmin)


def mip(msa):
    mi = buildMutinfoMatrix(msa)
    return applyMutinfoCorr(mi, corr="prod")


tmin = float("inf")
for _ in range(5):
    start = time.perf_counter()
    mip(msa)
    elapsed = time.perf_counter() - start
    if elapsed < tmin:
        tmin = elapsed

print_perf("Mutual Information APC", 1000 * tmin)

tmin = float("inf")
for _ in range(5):
    start = time.perf_counter()
    buildSeqidMatrix(msa)
    elapsed = time.perf_counter() - start
    if elapsed < tmin:
        tmin = elapsed

print_perf("Percent Identity Matrix", 1000 * tmin)

msa_long = parse_fasta("../data/PF00089_aligned.fasta", "PF00089", format=FASTA)
msa_wide = parse_fasta("../data/PF16957_aligned.fasta", "PF16957", format=FASTA)

bench_read("../data/PF00089.stockholm.gz", "compressed Stockholm PF00089", "PF00089", compressed=True, format=STOCKHOLM)
bench_read("../data/PF00089.sth", "uncompressed Stockholm PF00089", "PF00089", format=STOCKHOLM)
bench_read("../data/PF00089.fasta.gz", "compressed FASTA PF00089", "PF00089", compressed=True)
bench_read("../data/PF00089.fasta", "uncompressed FASTA PF00089", "PF00089")

bench_read("../data/PF16957.stockholm.gz", "compressed Stockholm PF16957", "PF16957", compressed=True, format=STOCKHOLM)
bench_read("../data/PF16957.sth", "uncompressed Stockholm PF16957", "PF16957", format=STOCKHOLM)
bench_read("../data/PF16957.fasta.gz", "compressed FASTA PF16957", "PF16957", compressed=True)
bench_read("../data/PF16957.fasta", "uncompressed FASTA PF16957", "PF16957")

start = time.perf_counter()
mip(msa_long)
elapsed = time.perf_counter() - start
print("[BENCH] MIp PF00089:", elapsed)

start = time.perf_counter()
mip(msa_wide)
elapsed = time.perf_counter() - start
print("[BENCH] MIp PF16957:", elapsed)

start = time.perf_counter()
calcShannonEntropy(msa_long)
elapsed = time.perf_counter() - start
print("[BENCH] Shannon entropy PF00089:", elapsed)

start = time.perf_counter()
calcShannonEntropy(msa_wide)
elapsed = time.perf_counter() - start
print("[BENCH] Shannon entropy PF16957:", elapsed)

start = time.perf_counter()
buildSeqidMatrix(msa_wide)
elapsed = time.perf_counter() - start
print("[BENCH] Percent Identity PF16957:", elapsed)
