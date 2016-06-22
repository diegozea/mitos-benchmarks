using MIToS.MSA
using MIToS.Information

const msa = read("../data/PF00089_aligned.fasta", FASTA)

estimateincolumns(msa, Int, ResidueProbability{Float64,1,false}, Entropy{Float64}())

println("[BENCH] Shannon entropy PF00089: ", @elapsed estimateincolumns(msa, Int, ResidueProbability{Float64,1,false}, Entropy{Float64}()))
