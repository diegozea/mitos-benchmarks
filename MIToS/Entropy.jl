using MIToS.MSA
using MIToS.Information

const msa_long = read("../data/PF00089_aligned.fasta", FASTA)
const msa_wide = read("../data/PF16957_aligned.fasta", FASTA)

estimateincolumns(msa_long, Int, ResidueProbability{Float64,1,false}, Entropy{Float64}())
estimateincolumns(msa_wide, Int, ResidueProbability{Float64,1,false}, Entropy{Float64}())

println("[BENCH] Shannon entropy PF00089: ", @elapsed estimateincolumns(msa_long, Int, ResidueProbability{Float64,1,false}, Entropy{Float64}()))
println("[BENCH] Shannon entropy PF16957: ", @elapsed estimateincolumns(msa_wide, Int, ResidueProbability{Float64,1,false}, Entropy{Float64}()))
