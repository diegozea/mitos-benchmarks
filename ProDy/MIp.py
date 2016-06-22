import time
from prody import *

msa = parseMSA("../data/PF00089_aligned.fasta")

def mip(msa):
    MI = buildMutinfoMatrix(msa)
    MIp = applyMutinfoCorr(MI)
    return MIp

start = time.time()
mip(msa)
elapsed = time.time() - start

print "[BENCH] MIp PF00089: ", elapsed
