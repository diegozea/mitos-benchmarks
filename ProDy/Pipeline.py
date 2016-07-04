import time
from prody import *

def print_perf(name, time):
    print("ProDy," + name + "," + str(time))

msa = parseMSA("../data/PF08171.fasta", format="fasta")

tmin = float('inf')
for i in range(5):
    t = time.time()
    parseMSA("../data/PF08171.sth", format="Stockholm")
    t = time.time()-t
    if t < tmin: tmin = t

print_perf("Read Pfam Stockholm MSA", 1000*tmin)

def mip(msa):
    MI = buildMutinfoMatrix(msa)
    MIp = applyMutinfoCorr(MI)
    return MIp

tmin = float('inf')
for i in range(5):
    t = time.time()
    mip(msa)
    t = time.time()-t
    if t < tmin: tmin = t

print_perf("Mutual Information APC", 1000*tmin)

tmin = float('inf')
for i in range(5):
    t = time.time()
    buildSeqidMatrix(msa)
    t = time.time()-t
    if t < tmin: tmin = t

print_perf("Percent Identity Matrix", 1000*tmin)
