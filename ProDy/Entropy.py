import time
from prody import *

msa_long = parseMSA("../data/PF00089_aligned.fasta")
msa_wide = parseMSA("../data/PF16957_aligned.fasta")

start = time.time()
calcShannonEntropy(msa_long)
elapsed = time.time() - start

print "[BENCH] Shannon entropy PF00089: ", elapsed


start = time.time()
calcShannonEntropy(msa_wide)
elapsed = time.time() - start

print "[BENCH] Shannon entropy PF16957: ", elapsed
