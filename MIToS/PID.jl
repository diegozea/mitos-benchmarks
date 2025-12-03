using MIToS.MSA

const data_dir = normpath(joinpath(@__DIR__, "..", "data"))

const msa_long = read_file(joinpath(data_dir, "PF00089_aligned.fasta"), FASTA)
const msa_wide = read_file(joinpath(data_dir, "PF16957_aligned.fasta"), FASTA)

_ = percentidentity(msa_long)
_ = percentidentity(msa_wide)

println("[BENCH] Percent Identity PF00089: ", @elapsed percentidentity(msa_long))
println("[BENCH] Percent Identity PF16957: ", @elapsed percentidentity(msa_wide))
