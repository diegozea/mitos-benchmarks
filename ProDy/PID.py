import time
from prody import buildSeqidMatrix, parseMSA

msa_pid = parseMSA("../data/PF16957_aligned.fasta")

start = time.perf_counter()
buildSeqidMatrix(msa_pid)
elapsed = time.perf_counter() - start
print("[BENCH] Percent Identity PF16957:", elapsed)
