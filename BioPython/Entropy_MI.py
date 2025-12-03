import time
from collections import Counter

import numpy as np
from Bio import AlignIO

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


def bench_entropy(path, label):
    aln = AlignIO.read(path, "fasta")
    arr = alignment_array(aln)
    start = time.perf_counter()
    shannon_entropy(arr)
    elapsed = time.perf_counter() - start
    print(f"[BENCH] Biopython Shannon entropy {label}:", elapsed)


def bench_mi(path, label):
    aln = AlignIO.read(path, "stockholm")
    arr = alignment_array(aln)
    start = time.perf_counter()
    mutual_information(arr)
    elapsed = time.perf_counter() - start
    print(f"[BENCH] Biopython Mutual Information {label}:", elapsed)


if __name__ == "__main__":
    bench_entropy("../data/PF00089_aligned.fasta", "PF00089")
    bench_entropy("../data/PF16957_aligned.fasta", "PF16957")
    bench_mi("../data/PF08171.sth", "PF08171")
