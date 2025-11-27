using MIToS.MSA
using MIToS.Information

const msa_long = read_file("../data/PF00089_aligned.fasta", FASTA)
const msa_wide = read_file("../data/PF16957_aligned.fasta", FASTA)

mip(msa) = last(
    buslje09(
        msa;
        samples = 0,           # no random shuffling (for comparability)
        clustering = false,    # no Hobohm clustering
        lambda = 0.0,          # no pseudocount
        maxgap = 1.0,          # keep all columns (match ProDy defaults)
        alphabet = UngappedAlphabet(),
    ),
)
_ = mip(msa_long) # warm-up
_ = mip(msa_wide)

println("[BENCH] MIp PF00089: ", @elapsed mip(msa_long))
println("[BENCH] MIp PF16957: ", @elapsed mip(msa_wide))
