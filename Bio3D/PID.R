library(bio3d)
library(microbenchmark)

# msa_long <- read.fasta("../data/PF00089_aligned.fasta") # It uses too much RAM
msa_wide <- read.fasta("../data/PF16957_aligned.fasta")

# It's more feature than ProDy and MIToS implementations
pid <- function() {
    seqidentity(msa_wide, normalize=FALSE)
}

bench <- microbenchmark(pid(), times=1)

cat("[BENCH] Percent Identity PF16957:", bench$time / 10^9, "\n", sep="")
