library(bio3d)
library(microbenchmark)

msa <- read.fasta("../data/PF00089_aligned.fasta")

# It's more feature than ProDy and MIToS implementations
h <- function() {
    entropy(msa)
}

bench <- microbenchmark(h(), times=1)

cat("[BENCH] Shannon entropy PF00089:", bench$time / 10^9, "\n", sep="")
