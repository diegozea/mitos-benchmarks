import time
from prody import applyMutinfoCorr, buildMutinfoMatrix, buildSeqidMatrix, parseMSA


def print_perf(name, timing_ms):
    print("ProDy," + name + "," + str(timing_ms))


msa = parseMSA("../data/PF08171.fasta", format="fasta")

tmin = float("inf")
for _ in range(5):
    start = time.perf_counter()
    parseMSA("../data/PF08171.sth", format="Stockholm")
    elapsed = time.perf_counter() - start
    if elapsed < tmin:
        tmin = elapsed

print_perf("Read Pfam Stockholm MSA", 1000 * tmin)


def mip(msa):
    mi = buildMutinfoMatrix(msa)
    return applyMutinfoCorr(mi, corr="prod")


tmin = float("inf")
for _ in range(5):
    start = time.perf_counter()
    mip(msa)
    elapsed = time.perf_counter() - start
    if elapsed < tmin:
        tmin = elapsed

print_perf("Mutual Information APC", 1000 * tmin)

tmin = float("inf")
for _ in range(5):
    start = time.perf_counter()
    buildSeqidMatrix(msa)
    elapsed = time.perf_counter() - start
    if elapsed < tmin:
        tmin = elapsed

print_perf("Percent Identity Matrix", 1000 * tmin)
