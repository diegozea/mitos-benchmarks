import time
from Bio import AlignIO


def bench_fasta(path, label):
    start = time.perf_counter()
    alignment = AlignIO.read(path, "fasta")
    elapsed = time.perf_counter() - start
    print(f"[BENCH] Biopython Read FASTA {label}:", elapsed)
    return alignment


if __name__ == "__main__":
    bench_fasta("../data/PF00089_aligned.fasta", "PF00089")
    bench_fasta("../data/PF16957_aligned.fasta", "PF16957")
