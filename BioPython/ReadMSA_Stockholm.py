import gzip
import time
from Bio import AlignIO


def read_stockholm(path):
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as handle:
            return list(AlignIO.parse(handle, "stockholm"))
    with open(path) as handle:
        return list(AlignIO.parse(handle, "stockholm"))


def bench_stockholm(path, label):
    start = time.perf_counter()
    alignments = read_stockholm(path)
    elapsed = time.perf_counter() - start
    print(f"[BENCH] Biopython Read Stockholm {label}:", elapsed)
    return alignments


if __name__ == "__main__":
    bench_stockholm("../data/PF08171.sth", "PF08171")
    bench_stockholm("../data/PF00089.stockholm.gz", "PF00089.gz")
