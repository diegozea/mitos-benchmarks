using MIToS.MSA

read("../data/PF00089.stockholm.gz", Stockholm, deletefullgaps=false)
read("../data/PF00089.sth", Stockholm, deletefullgaps=false)
read("../data/PF00089.fasta.gz", FASTA, deletefullgaps=false)
read("../data/PF00089.fasta", FASTA, deletefullgaps=false)

println("[BENCH] Read compressed Stockholm PF00089:", @elapsed read("../data/PF00089.stockholm.gz", Stockholm, deletefullgaps=false))
println("[BENCH] Read uncompressed Stockholm PF00089:", @elapsed read("../data/PF00089.sth", Stockholm, deletefullgaps=false))
println("[BENCH] Read compressed FASTA PF00089:", @elapsed read("../data/PF00089.fasta.gz", FASTA, deletefullgaps=false))
println("[BENCH] Read uncompressed FASTA PF00089:", @elapsed read("../data/PF00089.fasta", FASTA, deletefullgaps=false))
