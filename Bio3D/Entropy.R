library(bio3d)
library(microbenchmark)

msa_long <- read.fasta("../data/PF00089_aligned.fasta")
msa_wide <- read.fasta("../data/PF16957_aligned.fasta")

# It's more feature than ProDy and MIToS implementations
long <- function() {
    entropy(msa_long)
}
wide <- function() {
    entropy(msa_wide)
}

bench <- microbenchmark(long(), times=1)

cat("[BENCH] Shannon entropy PF00089:", bench$time / 10^9, "\n", sep="")

bench <- microbenchmark(wide(), times=1)

cat("[BENCH] Shannon entropy PF16957:", bench$time / 10^9, "\n", sep="")
