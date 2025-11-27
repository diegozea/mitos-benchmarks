using MIToS.MSA
using MIToS.Information

const msa_pid = read_file("../data/PF16957_aligned.fasta", FASTA)

percentidentity(msa_pid)

println("[BENCH] Percent Identity PF16957: ", @elapsed percentidentity(msa_pid))
