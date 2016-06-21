library(bio3d)
library(microbenchmark)

readfasta <- function() {
    read.fasta("../data/PF00089.fasta")
}

bench <- microbenchmark(readfasta(), times=1)

cat("[BENCH] Read uncompressed FASTA PF00089:", bench$time / 10^9, "\n", sep="")
