library(bio3d)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("--file=", args, value = TRUE)
script_dir <- if (length(script_arg) == 0) {
  getwd()
} else {
  dirname(normalizePath(sub("--file=", "", script_arg)))
}
root <- normalizePath(file.path(script_dir, ".."))
fasta_path <- file.path(root, "tests", "data", "test_alignment.fasta")

msa <- read.fasta(fasta_path)
pid <- seqidentity(msa, normalize = FALSE)

expected <- matrix(
  c(100, 75, 75,
    75, 100, 50,
    75, 50, 100),
  nrow = 3,
  byrow = TRUE
)

scale <- median(diag(pid))
normalized <- pid / scale * 100

if (!all(abs(normalized - expected) < 1e-6)) {
  stop("Bio3D percent identity does not match expected values")
}

cat("[TEST] Bio3D percent identity toy alignment ok\n")
