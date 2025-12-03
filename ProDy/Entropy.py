import time
from prody import calcShannonEntropy, parseMSA

msa_long = parseMSA("../data/PF00089_aligned.fasta")
msa_wide = parseMSA("../data/PF16957_aligned.fasta")

start = time.perf_counter()
calcShannonEntropy(msa_long)
elapsed = time.perf_counter() - start
print("[BENCH] Shannon entropy PF00089:", elapsed)

start = time.perf_counter()
calcShannonEntropy(msa_wide)
elapsed = time.perf_counter() - start
print("[BENCH] Shannon entropy PF16957:", elapsed)
