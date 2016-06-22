using MIToS.MSA
using MIToS.Information

const msa = read("../data/PF00089_aligned.fasta", FASTA)

mip(msa) = APC!(estimateincolumns(msa, ResidueCount{Int,2,true}, MutualInformation{Float64}()))

println("[BENCH] MIp PF00089: ", @elapsed mip(msa))
