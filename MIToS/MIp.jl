using MIToS.MSA
using MIToS.Information

const msa_long = read("../data/PF00089_aligned.fasta", FASTA)
const msa_wide = read("../data/PF16957_aligned.fasta", FASTA)

mip(msa) = APC!(mapcolpairfreq!(mutual_information, msa, Counts{Float64,2,GappedAlphabet}(ContingencyTable(Float64,Val{2},GappedAlphabet()))))
mip(msa_long)
mip(msa_wide)

println("[BENCH] MIp PF00089: ", @elapsed mip(msa_long))
println("[BENCH] MIp PF16957: ", @elapsed mip(msa_wide))
