library(bio3d)
library(microbenchmark)

# PF00089 (20k seqs) is too large for seqidentity on typical machines; report NA
cat("[BENCH] Percent Identity PF00089: NA (dataset too large for seqidentity here)\n")

msa_wide <- read.fasta("../data/PF16957_aligned.fasta")
pid_wide <- function() {
    seqidentity(msa_wide, normalize=FALSE)
}

bench <- microbenchmark(pid_wide(), times=1)
cat("[BENCH] Percent Identity PF16957:", bench$time / 10^9, "\n", sep="")
