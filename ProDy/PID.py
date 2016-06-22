import time
from prody import *

msa_long = parseMSA("../data/PF00089_aligned.fasta")
msa_wide = parseMSA("../data/PF16957_aligned.fasta")

start = time.time()
buildSeqidMatrix(msa_long)
elapsed = time.time() - start

print "[BENCH] Percent Identity PF00089: ", elapsed

start = time.time()
buildSeqidMatrix(msa_wide)
elapsed = time.time() - start

print "[BENCH] Percent Identity PF16957: ", elapsed
