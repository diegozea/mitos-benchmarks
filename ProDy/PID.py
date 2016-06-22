import time
from prody import *

msa = parseMSA("../data/PF00089_aligned.fasta")

start = time.time()
buildSeqidMatrix(msa)
elapsed = time.time() - start

print "[BENCH] Percent Identity PF00089: ", elapsed
