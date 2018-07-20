using MIToS.MSA
using MIToS.Information

const msa_long = read("../data/PF00089_aligned.fasta", FASTA)
const msa_wide = read("../data/PF16957_aligned.fasta", FASTA)

mapcolfreq!(entropy, msa_long, Probabilities{Float64,1,UngappedAlphabet}(ContingencyTable(Float64,Val{1},UngappedAlphabet())))
mapcolfreq!(entropy, msa_wide, Probabilities{Float64,1,UngappedAlphabet}(ContingencyTable(Float64,Val{1},UngappedAlphabet())))

println("[BENCH] Shannon entropy PF00089: ", @elapsed mapcolfreq!(entropy, msa_long, Probabilities{Float64,1,UngappedAlphabet}(ContingencyTable(Float64,Val{1},UngappedAlphabet()))))
println("[BENCH] Shannon entropy PF16957: ", @elapsed mapcolfreq!(entropy, msa_wide, Probabilities{Float64,1,UngappedAlphabet}(ContingencyTable(Float64,Val{1},UngappedAlphabet()))))
