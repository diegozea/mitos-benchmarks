library(bio3d)
library(microbenchmark)

pdb_path <- "../data/4BL0.pdb"
if (!file.exists(pdb_path)) {
  stop("PDB file not found: ", pdb_path)
}

pdb <- read.pdb(pdb_path)
ca_inds <- atom.select(pdb, "calpha")
ca_xyz <- pdb$xyz[ca_inds$xyz]

bench_maps <- microbenchmark(
  contact = cmap(ca_xyz, dcut = 8.0, scut = 0, mask.lower = FALSE, mask.upper = FALSE),
  distance = dm.xyz(ca_xyz),
  times = 1
)

exprs <- as.character(bench_maps$expr)
contact_time <- bench_maps$time[exprs == "contact"] / 1e9
distance_time <- bench_maps$time[exprs == "distance"] / 1e9

contact_map <- cmap(ca_xyz, dcut = 8.0, scut = 0, mask.lower = FALSE, mask.upper = FALSE)
contact_mat <- if (is.list(contact_map) && !is.null(contact_map$map)) contact_map$map else contact_map
distance_mat <- dm.xyz(ca_xyz)

cat("[BENCH] Contact Map 4BL0 (Bio3D CA, 8A):", contact_time, "\n")
cat("[BENCH] Distance Matrix 4BL0 (Bio3D CA):", distance_time, "\n")

bench_dccm <- microbenchmark(
  dccm_nma = dccm.nma(pdb),
  times = 1
)
dccm_time <- bench_dccm$time / 1e9
cij <- dccm.nma(pdb)
cat("[BENCH] DCCM (NMA) 4BL0:", dccm_time, "\n")

if (requireNamespace("igraph", quietly = TRUE)) {
  net_time <- system.time({
    net <- cna(cij, cm = contact_mat)
  })["elapsed"]
  cat("[BENCH] CNA network (contact-filtered) 4BL0:", net_time, "\n")
  cat("[INFO] CNA nodes:", nrow(net$acc), "edges:", nrow(net$edge), "\n")
} else {
  cat("[SKIP] igraph not installed; skipping CNA network benchmark\n")
}
