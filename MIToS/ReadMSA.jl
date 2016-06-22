using MIToS.MSA

read("../data/PF00089.stockholm.gz", Stockholm, MultipleSequenceAlignment, deletefullgaps=false)
read("../data/PF00089.sth", Stockholm, MultipleSequenceAlignment, deletefullgaps=false)
read("../data/PF00089.fasta.gz", FASTA, MultipleSequenceAlignment, deletefullgaps=false)
read("../data/PF00089.fasta", FASTA, MultipleSequenceAlignment, deletefullgaps=false)

println("[BENCH] Read compressed Stockholm PF00089: ", @elapsed read("../data/PF00089.stockholm.gz", Stockholm, MultipleSequenceAlignment, deletefullgaps=false))
println("[BENCH] Read uncompressed Stockholm PF00089: ", @elapsed read("../data/PF00089.sth", Stockholm, MultipleSequenceAlignment, deletefullgaps=false))
println("[BENCH] Read compressed FASTA PF00089: ", @elapsed read("../data/PF00089.fasta.gz", FASTA, MultipleSequenceAlignment, deletefullgaps=false))
println("[BENCH] Read uncompressed FASTA PF00089: ", @elapsed read("../data/PF00089.fasta", FASTA, MultipleSequenceAlignment, deletefullgaps=false))

read("../data/PF16957.stockholm.gz", Stockholm, MultipleSequenceAlignment, deletefullgaps=false)
read("../data/PF16957.sth", Stockholm, MultipleSequenceAlignment, deletefullgaps=false)
read("../data/PF16957.fasta.gz", FASTA, MultipleSequenceAlignment, deletefullgaps=false)
read("../data/PF16957.fasta", FASTA, MultipleSequenceAlignment, deletefullgaps=false)

println("[BENCH] Read compressed Stockholm PF16957: ", @elapsed read("../data/PF16957.stockholm.gz", Stockholm, MultipleSequenceAlignment, deletefullgaps=false))
println("[BENCH] Read uncompressed Stockholm PF16957: ", @elapsed read("../data/PF16957.sth", Stockholm, MultipleSequenceAlignment, deletefullgaps=false))
println("[BENCH] Read compressed FASTA PF16957: ", @elapsed read("../data/PF16957.fasta.gz", FASTA, MultipleSequenceAlignment, deletefullgaps=false))
println("[BENCH] Read uncompressed FASTA PF16957: ", @elapsed read("../data/PF16957.fasta", FASTA, MultipleSequenceAlignment, deletefullgaps=false))
