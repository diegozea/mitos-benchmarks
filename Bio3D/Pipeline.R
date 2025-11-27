library(bio3d)
library(microbenchmark)
library(compiler)

if(abs(system.time(Sys.sleep(1))["elapsed"] - 1) < 0.01){ mult <<- 1000 }else{cat("Unit ERROR\n")}

check_msa <- function(msa, expected_seq, expected_len, label) {
    if (nrow(msa$ali) != expected_seq || ncol(msa$ali) != expected_len) {
        stop(paste(label, "dimensions changed:", nrow(msa$ali), "x", ncol(msa$ali)))
    }
}

timeit = function(name, f, ..., times=5) {
    tmin = Inf
    f = cmpfun(f)
    for (t in 1:times) {
        t = system.time(f(...))["elapsed"]
        if (t < tmin) tmin = t
    }
    cat(sprintf("r,%s,%.8f\n", name, mult * tmin))
}

pdb <- read.pdb("../data/4BL0.pdb")
cm <- function(){
    cmap(pdb$xyz, grpby=pdb$atom$resno, dcut=6.03, scut=0)
    }

timeit("Protein Contact Map", cm)

msa <- read.fasta("../data/PF08171.fasta")
check_msa(msa, 208, 68, "PF08171")
pid <- function() {
    seqidentity(msa, normalize=FALSE)
}

timeit("Percent Identity Matrix", pid)

bench_read <- function(path, label) {
    bench <- microbenchmark(read.fasta(path), times=1)
    cat("[BENCH] Read uncompressed FASTA ", label, ":", bench$time / 10^9, "\n", sep="")
}

bench_read("../data/PF00089_aligned.fasta", "PF00089")
bench_read("../data/PF16957_aligned.fasta", "PF16957")

msa_long <- read.fasta("../data/PF00089_aligned.fasta")
msa_wide <- read.fasta("../data/PF16957_aligned.fasta")
check_msa(msa_long, 20152, 221, "PF00089 aligned")
check_msa(msa_wide, 154, 548, "PF16957 aligned")

entropy_long <- function() {
    entropy(msa_long)
}
entropy_wide <- function() {
    entropy(msa_wide)
}

bench <- microbenchmark(entropy_long(), times=1)
cat("[BENCH] Shannon entropy PF00089:", bench$time / 10^9, "\n", sep="")

bench <- microbenchmark(entropy_wide(), times=1)
cat("[BENCH] Shannon entropy PF16957:", bench$time / 10^9, "\n", sep="")

cat("[BENCH] Percent Identity PF00089: NA (dataset too large for seqidentity here)\n")

pid_wide <- function() {
    seqidentity(msa_wide, normalize=FALSE)
}

bench <- microbenchmark(pid_wide(), times=1)
cat("[BENCH] Percent Identity PF16957:", bench$time / 10^9, "\n", sep="")
