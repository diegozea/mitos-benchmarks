using MIToS.MSA

const data_dir = normpath(joinpath(@__DIR__, "..", "data"))

bench_read(path, format, label) = println(label, @elapsed read_file(path, format))

bench_read(joinpath(data_dir, "PF00089.stockholm.gz"), Stockholm, "[BENCH] Read compressed Stockholm PF00089: ")
bench_read(joinpath(data_dir, "PF00089.sth"), Stockholm, "[BENCH] Read uncompressed Stockholm PF00089: ")
bench_read(joinpath(data_dir, "PF00089.fasta.gz"), FASTA, "[BENCH] Read compressed FASTA PF00089: ")
bench_read(joinpath(data_dir, "PF00089.fasta"), FASTA, "[BENCH] Read uncompressed FASTA PF00089: ")

bench_read(joinpath(data_dir, "PF16957.stockholm.gz"), Stockholm, "[BENCH] Read compressed Stockholm PF16957: ")
bench_read(joinpath(data_dir, "PF16957.sth"), Stockholm, "[BENCH] Read uncompressed Stockholm PF16957: ")
bench_read(joinpath(data_dir, "PF16957.fasta.gz"), FASTA, "[BENCH] Read compressed FASTA PF16957: ")
bench_read(joinpath(data_dir, "PF16957.fasta"), FASTA, "[BENCH] Read uncompressed FASTA PF16957: ")
