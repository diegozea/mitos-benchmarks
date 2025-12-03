import gzip
import time
from collections import Counter

import numpy as np
from Bio import AlignIO

EXPECTED = {
    "PF08171": (208, 68),
    "PF00089": (20152, 221),
    "PF16957": (154, 548),
}

def print_perf(name, time):
    print("BioPython," + name + "," + str(time))

tmin = float('inf')
for i in range(5):
    t = time.time()
    AlignIO.read("../data/PF08171.sth", "stockholm")
    t = time.time()-t
    if t < tmin: tmin = t

print_perf("Read Pfam Stockholm MSA", 1000*tmin)


def read_stockholm(path):
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as handle:
            return list(AlignIO.parse(handle, "stockholm"))
    with open(path) as handle:
        return list(AlignIO.parse(handle, "stockholm"))


def assert_shape(arr, expected_key, label):
    exp_seq, exp_len = EXPECTED[expected_key]
    assert arr.shape[0] == exp_seq, f"{label}: expected {exp_seq} sequences, got {arr.shape[0]}"
    assert arr.shape[1] == exp_len, f"{label}: expected {exp_len} columns, got {arr.shape[1]}"


def match_mask(arr):
    lower_or_gap = np.char.islower(arr) | (arr == "-") | (arr == ".")
    return ~np.all(lower_or_gap, axis=0)


def trim_alignment_array(aln, expected_key):
    arr = alignment_array(aln)
    keep = match_mask(arr)
    trimmed = arr[:, keep]
    assert_shape(trimmed, expected_key, f"{expected_key} trimmed")
    return trimmed


def bench_stockholm(path, label, expected_key):
    start = time.perf_counter()
    alignments = read_stockholm(path)
    if not alignments:
        raise RuntimeError(f"No alignments parsed from {path}")
    arr = trim_alignment_array(alignments[0], expected_key)
    elapsed = time.perf_counter() - start
    print(f"[BENCH] Biopython Read Stockholm {label}:", elapsed)
    return arr


def bench_fasta(path, label):
    start = time.perf_counter()
    aln = AlignIO.read(path, "fasta")
    arr = alignment_array(aln)
    expected_key = "PF08171" if "08171" in path else "PF00089" if "00089" in path else "PF16957"
    assert_shape(arr, expected_key, f"{label} FASTA")
    elapsed = time.perf_counter() - start
    print(f"[BENCH] Biopython Read FASTA {label}:", elapsed)
    return arr


GAP_TOKENS = ("-", ".")


def alignment_array(aln):
    return np.array([list(str(record.seq)) for record in aln], dtype="U1")


def shannon_entropy(arr):
    ent = np.zeros(arr.shape[1], dtype=float)
    for idx, col in enumerate(arr.T):
        mask = ~np.isin(col, GAP_TOKENS)
        if not mask.any():
            continue
        symbols, counts = np.unique(col[mask], return_counts=True)
        probs = counts / counts.sum()
        ent[idx] = float(-(probs * np.log2(probs)).sum())
    return ent


def mutual_information(arr):
    arr = np.asarray(arr)
    ncol = arr.shape[1]
    mi = np.zeros((ncol, ncol), dtype=float)
    for i in range(ncol - 1):
        col_i = arr[:, i]
        mask_i = ~np.isin(col_i, GAP_TOKENS)
        for j in range(i + 1, ncol):
            col_j = arr[:, j]
            mask = mask_i & ~np.isin(col_j, GAP_TOKENS)
            if not mask.any():
                continue
            pairs = Counter(zip(col_i[mask], col_j[mask]))
            total = float(sum(pairs.values()))
            if total == 0:
                continue
            pi = Counter(col_i[mask])
            pj = Counter(col_j[mask])
            mi_val = 0.0
            for (ai, aj), count in pairs.items():
                p_ij = count / total
                p_i = pi[ai] / total
                p_j = pj[aj] / total
                mi_val += p_ij * np.log2(p_ij / (p_i * p_j))
            mi[i, j] = mi[j, i] = mi_val
    return mi


def _to_array(aln):
    if isinstance(aln, np.ndarray):
        return aln
    return np.array([list(str(rec.seq)) for rec in aln], dtype="U1")


def pid_matrix(aln, block=128):
    arr = _to_array(aln)
    n, _ = arr.shape
    mask = ~np.isin(arr, GAP_TOKENS)
    result = np.zeros((n, n), dtype=np.float32)

    for i in range(0, n, block):
        bi = arr[i : i + block]
        mi = mask[i : i + block]
        for j in range(i, n, block):
            bj = arr[j : j + block]
            mj = mask[j : j + block]
            both_mask = mi[:, None, :] & mj[None, :, :]
            matches = np.sum((bi[:, None, :] == bj[None, :, :]) & both_mask, axis=2, dtype=np.int32)
            aligned = np.sum(both_mask, axis=2, dtype=np.int32)
            pid_block = np.divide(
                matches,
                aligned,
                out=np.zeros_like(matches, dtype=np.float32),
                where=aligned > 0,
            )
            pid_block *= 100.0
            result[i : i + pid_block.shape[0], j : j + pid_block.shape[1]] = pid_block
            if j != i:
                result[j : j + pid_block.shape[1], i : i + pid_block.shape[0]] = pid_block.T
    return result


def bench_entropy(path, label):
    aln = AlignIO.read(path, "fasta")
    arr = alignment_array(aln)
    expected_key = "PF00089" if "00089" in path else "PF16957"
    assert_shape(arr, expected_key, f"{label} FASTA")
    start = time.perf_counter()
    shannon_entropy(arr)
    elapsed = time.perf_counter() - start
    print(f"[BENCH] Biopython Shannon entropy {label}:", elapsed)


def bench_mi(path, label):
    aln = AlignIO.read(path, "stockholm")
    arr = trim_alignment_array(aln, "PF08171")
    start = time.perf_counter()
    mutual_information(arr)
    elapsed = time.perf_counter() - start
    print(f"[BENCH] Biopython Mutual Information {label}:", elapsed)


def bench_pid(path, label):
    aln = AlignIO.read(path, "fasta")
    arr = _to_array(aln)
    expected_key = "PF08171" if "08171" in path else "PF16957" if "16957" in path else "PF00089"
    assert_shape(arr, expected_key, f"{label} FASTA")
    start = time.perf_counter()
    pid_matrix(arr)
    elapsed = time.perf_counter() - start
    print(f"[BENCH] Biopython Percent Identity {label}:", elapsed)


def bench_annotations(path, label):
    alignments = read_stockholm(path)
    if not alignments:
        raise RuntimeError("No alignments parsed from Stockholm file")
    alignment = alignments[0]
    assert len(alignment) == EXPECTED["PF08171"][0], "PF08171 Stockholm sequence count changed"
    assert alignment.get_alignment_length() >= EXPECTED["PF08171"][1], "PF08171 Stockholm length changed"
    ann = alignment.annotations
    col_ann = alignment.column_annotations
    seqs = list(alignment)

    ss_cons_key = "consensus secondary structure"
    if ss_cons_key not in col_ann and "secondary_structure" in col_ann:
        ss_cons_key = "secondary_structure"
    residue_ss_key = "secondary structure"
    if residue_ss_key not in seqs[0].letter_annotations and "secondary_structure" in seqs[0].letter_annotations:
        residue_ss_key = "secondary_structure"

    start = time.perf_counter()
    _ = ann.get("accession")
    _ = ann.get("identifier")
    _ = col_ann.get(ss_cons_key)
    _ = [record.annotations.get("accession") for record in seqs]
    _ = [record.letter_annotations.get(residue_ss_key) for record in seqs]
    access_elapsed = time.perf_counter() - start
    print(f"[BENCH] Biopython Stockholm {label} annotation access:", access_elapsed)


bench_fasta("../data/PF00089_aligned.fasta", "PF00089")
bench_fasta("../data/PF16957_aligned.fasta", "PF16957")

bench_stockholm("../data/PF08171.sth", "PF08171", "PF08171")
bench_stockholm("../data/PF00089.stockholm.gz", "PF00089.gz", "PF00089")
bench_stockholm("../data/PF16957.sth", "PF16957", "PF16957")

bench_annotations("../data/PF08171.sth", "PF08171")

bench_entropy("../data/PF00089_aligned.fasta", "PF00089")
bench_entropy("../data/PF16957_aligned.fasta", "PF16957")
bench_mi("../data/PF08171.sth", "PF08171")

bench_pid("../data/PF08171.fasta", "PF08171")
bench_pid("../data/PF16957_aligned.fasta", "PF16957")
