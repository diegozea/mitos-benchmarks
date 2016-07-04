import time
from Bio import AlignIO

def print_perf(name, time):
    print("BioPython," + name + "," + str(time))

tmin = float('inf')
for i in range(5):
    t = time.time()
    AlignIO.read("../data/PF08171.sth", "stockholm")
    t = time.time()-t
    if t < tmin: tmin = t

print_perf("Read Pfam Stockholm MSA", 1000*tmin)

