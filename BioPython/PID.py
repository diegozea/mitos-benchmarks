import numpy as np
import time
from Bio import AlignIO

GAP_TOKENS = ("-", ".")


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


def bench_pid(path, label):
    aln = AlignIO.read(path, "fasta")
    arr = _to_array(aln)
    start = time.perf_counter()
    pid_matrix(arr)
    elapsed = time.perf_counter() - start
    print(f"[BENCH] Biopython Percent Identity {label}:", elapsed)


if __name__ == "__main__":
    bench_pid("../data/PF16957_aligned.fasta", "PF16957")
