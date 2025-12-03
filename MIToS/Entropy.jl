using MIToS.MSA
using MIToS.Information

const data_dir = normpath(joinpath(@__DIR__, "..", "data"))

const msa_long = read_file(joinpath(data_dir, "PF00089_aligned.fasta"), FASTA)
const msa_wide = read_file(joinpath(data_dir, "PF16957_aligned.fasta"), FASTA)

table = Frequencies(ContingencyTable(Float64, Val{1}, UngappedAlphabet()))

_ = mapcolfreq!(shannon_entropy, msa_long, table)
_ = mapcolfreq!(shannon_entropy, msa_wide, table)

println(
    "[BENCH] Shannon entropy PF00089: ",
    @elapsed mapcolfreq!(shannon_entropy, msa_long, table)
)
println(
    "[BENCH] Shannon entropy PF16957: ",
    @elapsed mapcolfreq!(shannon_entropy, msa_wide, table)
)
