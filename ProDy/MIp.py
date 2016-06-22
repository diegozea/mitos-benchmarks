import time
from prody import *

msa_long = parseMSA("../data/PF00089_aligned.fasta")
msa_wide = parseMSA("../data/PF16957_aligned.fasta")

def mip(msa):
    MI = buildMutinfoMatrix(msa)
    MIp = applyMutinfoCorr(MI)
    return MIp

start = time.time()
mip(msa_long)
elapsed = time.time() - start

print "[BENCH] MIp PF00089: ", elapsed

start = time.time()
mip(msa_wide)
elapsed = time.time() - start

print "[BENCH] MIp PF16957: ", elapsed
