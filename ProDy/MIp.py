import time
from prody import applyMutinfoCorr, buildMutinfoMatrix, parseMSA

msa_long = parseMSA("../data/PF00089_aligned.fasta")
msa_wide = parseMSA("../data/PF16957_aligned.fasta")


def mip(msa):
    mi = buildMutinfoMatrix(msa)
    return applyMutinfoCorr(mi, corr="prod")


start = time.perf_counter()
mip(msa_long)
elapsed = time.perf_counter() - start
print("[BENCH] MIp PF00089:", elapsed)

start = time.perf_counter()
mip(msa_wide)
elapsed = time.perf_counter() - start
print("[BENCH] MIp PF16957:", elapsed)
