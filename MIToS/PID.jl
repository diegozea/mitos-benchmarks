using MIToS.MSA
using MIToS.Information

const msa = read("../data/PF00089_aligned.fasta", FASTA)

percentidentity(msa)

println("[BENCH] Percent Identity PF00089: ", @elapsed percentidentity(msa))
