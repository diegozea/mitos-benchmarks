import time
from prody import *

start = time.time()
parseMSA("../data/PF00089.stockholm.gz", compressed=True, format="Stockholm")
elapsed = time.time() - start

print "[BENCH] Read compressed Stockholm PF00089: ", elapsed

start = time.time()
parseMSA("../data/PF00089.sth",) # .stockholm causes an error
elapsed = time.time() - start

print "[BENCH] Read uncompressed Stockholm PF00089: ", elapsed

start = time.time()
parseMSA("../data/PF00089.fasta.gz", compressed=True)
elapsed = time.time() - start

print "[BENCH] Read compressed FASTA PF00089: ", elapsed

start = time.time()
parseMSA("../data/PF00089.fasta")
elapsed = time.time() - start

print "[BENCH] Read uncompressed FASTA PF00089: ", elapsed

start = time.time()
parseMSA("../data/PF16957.stockholm.gz", compressed=True, format="Stockholm")
elapsed = time.time() - start

print "[BENCH] Read compressed Stockholm PF16957: ", elapsed

start = time.time()
parseMSA("../data/PF16957.sth",) # .stockholm causes an error
elapsed = time.time() - start

print "[BENCH] Read uncompressed Stockholm PF16957: ", elapsed

start = time.time()
parseMSA("../data/PF16957.fasta.gz", compressed=True)
elapsed = time.time() - start

print "[BENCH] Read compressed FASTA PF16957: ", elapsed

start = time.time()
parseMSA("../data/PF16957.fasta")
elapsed = time.time() - start

print "[BENCH] Read uncompressed FASTA PF16957: ", elapsed
