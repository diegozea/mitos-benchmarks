import time
from prody import *

msa = parseMSA("../data/PF00089_aligned.fasta")

start = time.time()
calcShannonEntropy(msa)
elapsed = time.time() - start

print "[BENCH] Shannon entropy PF00089: ", elapsed
