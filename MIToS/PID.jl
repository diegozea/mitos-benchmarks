using MIToS.MSA
using MIToS.Information

const msa_long = read("../data/PF00089_aligned.fasta", FASTA)
const msa_wide = read("../data/PF16957_aligned.fasta", FASTA)

percentidentity(msa_long)
percentidentity(msa_wide)

println("[BENCH] Percent Identity PF00089: ", @elapsed percentidentity(msa_long))
println("[BENCH] Percent Identity PF16957: ", @elapsed percentidentity(msa_wide))
