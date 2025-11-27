try
    import BioStructures
catch e
    println("[SKIP] BioStructures benchmarks: BioStructures.jl not available ($e)")
    exit()
end

graphs_available = true
try
    import Graphs
catch e
    graphs_available = false
    println("[SKIP] Graphs.jl not available; graph benchmark will be skipped ($e)")
end

meta_available = true
try
    import MetaGraphs
catch e
    meta_available = false
    println("[SKIP] MetaGraphs.jl not available; contact graph benchmark will be skipped ($e)")
end

using BioStructures: calphaselector, collectatoms, ContactMap, DistanceMap

data_dir = normpath(joinpath(@__DIR__, "..", "data"))
pdb_path = joinpath(data_dir, "4BL0.pdb")

if !isfile(pdb_path)
    println("[SKIP] BioStructures benchmarks: missing file $pdb_path")
    exit()
end

struc = BioStructures.read(pdb_path, BioStructures.PDB)
chain_id = "B"

chain = try
    struc[chain_id]
catch e
    println("[SKIP] BioStructures benchmarks: chain $chain_id not present ($e)")
    exit()
end

calphas = collectatoms(chain, calphaselector)

contact_time = @elapsed cmap = ContactMap(calphas; cutoff = 8.0)
distance_time = @elapsed dmap = DistanceMap(calphas)

println("[BENCH] BioStructures ContactMap 4BL0 chain $chain_id CA (8A): ", contact_time)
println("[BENCH] BioStructures DistanceMap 4BL0 chain $chain_id CA: ", distance_time)

if graphs_available && meta_available
    try
        g = nothing
        graph_time = @elapsed g = MetaGraphs.MetaGraph(cmap)
        println("[BENCH] BioStructures ContactGraph 4BL0 chain $chain_id: ", graph_time)
        println("[INFO] BioStructures graph nodes=", Graphs.nv(g), " edges=", Graphs.ne(g))
    catch e
        println("[SKIP] Contact graph benchmark failed: $e")
    end
end
