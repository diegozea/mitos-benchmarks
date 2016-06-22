library(bio3d)
library(microbenchmark)

msa <- read.fasta("../data/PF00089_aligned.fasta")

# It's more feature than ProDy and MIToS implementations
pid <- function() {
    seqidentity(msa, normalize=FALSE)
}

bench <- microbenchmark(pid(), times=1)

cat("[BENCH] Percent Identity PF00089:", bench$time / 10^9, "\n", sep="")
