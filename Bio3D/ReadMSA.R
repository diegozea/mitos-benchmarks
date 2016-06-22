library(bio3d)
library(microbenchmark)

readfastalong <- function() {
    read.fasta("../data/PF00089.fasta")
}

readfastawide <- function() {
    read.fasta("../data/PF16957.fasta")
}

bench <- microbenchmark(readfastalong(), times=1)

cat("[BENCH] Read uncompressed FASTA PF00089:", bench$time / 10^9, "\n", sep="")

bench <- microbenchmark(readfastawide(), times=1)

cat("[BENCH] Read uncompressed FASTA PF16957:", bench$time / 10^9, "\n", sep="")
