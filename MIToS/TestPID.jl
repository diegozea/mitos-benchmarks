using MIToS.MSA
using LinearAlgebra
using Statistics

root = normpath(joinpath(@__DIR__, ".."))
fasta_path = joinpath(root, "tests", "data", "test_alignment.fasta")

msa = read_file(fasta_path, FASTA)
pid = percentidentity(msa)

expected = [100.0 75.0 75.0; 75.0 100.0 50.0; 75.0 50.0 100.0]

scale = median(diag(pid))
@assert isapprox(scale, 100.0; atol = 1e-6) "Diagonal of PID matrix deviates from 100"

@assert all(isapprox.(pid, expected; atol=1e-6))
println("[TEST] MIToS percent identity toy alignment ok")
