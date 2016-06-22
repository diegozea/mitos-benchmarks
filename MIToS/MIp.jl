using MIToS.MSA
using MIToS.Information

const msa_long = read("../data/PF00089_aligned.fasta", FASTA)
const msa_wide = read("../data/PF16957_aligned.fasta", FASTA)

mip(msa) = APC!(estimateincolumns(msa, ResidueCount{Int,2,true}, MutualInformation{Float64}()))
mip(msa_long)
mip(msa_wide)

println("[BENCH] MIp PF00089: ", @elapsed mip(msa_long))
println("[BENCH] MIp PF16957: ", @elapsed mip(msa_wide))
