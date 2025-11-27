using MIToS.MSA
using MIToS.Pfam
using MIToS.PDB
using MIToS.Information
using Printf
using ROCAnalysis

function bench_ms(label::AbstractString, f::Function)
    # warm up once, then take minimum of 5 timings (in ms)
    f()
    times = Float64[]
    for _ in 1:5
        push!(times, 1000 * @elapsed f())
    end
    @printf "MIToS,%s,%f\n" label minimum(times)
end

aln_mapping = read_file("../data/PF08171.sth", Stockholm; generatemapping = true, useidcoordinates = true)
col2res = msacolumn2pdbresidue(aln_mapping, "BUB1_YEAST/291-355", "4BL0", "B", "PF08171", "../data/4bl0.xml.gz")
pdb_residues = read_file("../data/4BL0.xml", PDBML)
resdict = residuesdict(pdb_residues; model = "1", chain = "B", group = "ATOM", residue = All)
cmap = msacontacts(aln_mapping, resdict, col2res)
_, MIp = buslje09(aln_mapping; samples = 0, clustering = false, lambda = 0.0, alphabet = UngappedAlphabet(), maxgap = 1.0)

bench_ms("Read Pfam Stockholm MSA", () -> read_file("../data/PF08171.sth", Stockholm, Matrix{Residue}))

bench_ms("Read MSA and annotations", () -> read_file("../data/PF08171.sth", Stockholm))

bench_ms(
    "Read MSA and annotations, generate coordinates",
    () -> read_file("../data/PF08171.sth", Stockholm; generatemapping = true, useidcoordinates = true),
)

bench_ms("Percent Identity Matrix", () -> percentidentity(aln_mapping))

bench_ms(
    "SIFTS residue level mapping",
    () -> msacolumn2pdbresidue(aln_mapping, "BUB1_YEAST/291-355", "4BL0", "B", "PF08171", "../data/4bl0.xml.gz"),
)

bench_ms("Read PDBML", () -> read_file("../data/4BL0.xml", PDBML))

bench_ms("Protein Contact Map", () -> msacontacts(aln_mapping, resdict, col2res))

bench_ms(
    "Mutual Information APC",
    () -> APC!(mapcolpairfreq!(mutual_information, aln_mapping, Counts{Float64,2,GappedAlphabet}(ContingencyTable(Float64, Val{2}, GappedAlphabet())))),
)

bench_ms("AUC (ROC) for contact prediction, MIp", () -> AUC(MIp, cmap))
