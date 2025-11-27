import time
import numpy as np
from Bio import AlignIO


def pid_matrix(arr, block=128):
    n, _ = arr.shape
    mask = arr != "-"
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
    arr = np.array([list(str(rec.seq)) for rec in aln], dtype="U1")

    start = time.perf_counter()
    pid_matrix(arr)
    elapsed = time.perf_counter() - start
    print(f"[BENCH] Biopython Percent Identity {label}:", elapsed)


if __name__ == "__main__":
    bench_pid("../data/PF16957_aligned.fasta", "PF16957")
