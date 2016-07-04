library(bio3d)
library(compiler)

if(abs(system.time(Sys.sleep(1))["elapsed"] - 1) < 0.01){ mult <<- 1000 }else{cat("Unit ERROR\n")}

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
pid <- function() {
    seqidentity(msa, normalize=FALSE)
}

timeit("Percent Identity Matrix", pid)
