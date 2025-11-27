using MIToS.MSA
using MIToS.Pfam
using MIToS.PDB
using MIToS.Information
using Printf
using ROCAnalysis

const data_dir = normpath(joinpath(@__DIR__, "..", "data"))
const EXPECT_SEQ_PF08171 = 208
const EXPECT_LEN_PF08171 = 68
const EXPECT_SEQ_PF00089 = 20152
const EXPECT_LEN_PF00089 = 221
const EXPECT_SEQ_PF16957 = 154
const EXPECT_LEN_PF16957 = 548

function bench_ms(label::AbstractString, f::Function)
    # warm up once, then take minimum of 5 timings (in ms)
    f()
    times = Float64[]
    for _ in 1:5
        push!(times, 1000 * @elapsed f())
    end
    @printf "MIToS,%s,%f\n" label minimum(times)
end

bench_read(path, format, label) = println(label, @elapsed read_file(path, format))

aln_mapping = read_file(joinpath(data_dir, "PF08171.sth"), Stockholm; generatemapping = true, useidcoordinates = true)
@assert size(aln_mapping) == (EXPECT_SEQ_PF08171, EXPECT_LEN_PF08171) "PF08171 mapping dimensions changed"
col2res = msacolumn2pdbresidue(aln_mapping, "BUB1_YEAST/291-355", "4BL0", "B", "PF08171", joinpath(data_dir, "4bl0.xml.gz"))
pdb_residues = read_file(joinpath(data_dir, "4BL0.xml"), PDBML)
resdict = residuesdict(pdb_residues; model = "1", chain = "B", group = "ATOM", residue = All)
cmap = msacontacts(aln_mapping, resdict, col2res)
_, MIp = buslje09(aln_mapping; samples = 0, clustering = false, lambda = 0.0, alphabet = UngappedAlphabet(), maxgap = 1.0, apc = true)

bench_ms("Read Pfam Stockholm MSA", () -> read_file(joinpath(data_dir, "PF08171.sth"), Stockholm, Matrix{Residue}))

bench_ms("Read MSA and annotations", () -> read_file(joinpath(data_dir, "PF08171.sth"), Stockholm))

bench_ms(
    "Read MSA and annotations, generate coordinates",
    () -> read_file(joinpath(data_dir, "PF08171.sth"), Stockholm; generatemapping = true, useidcoordinates = true),
)

bench_ms("Percent Identity Matrix", () -> percentidentity(aln_mapping))

bench_ms(
    "SIFTS residue level mapping",
    () -> msacolumn2pdbresidue(aln_mapping, "BUB1_YEAST/291-355", "4BL0", "B", "PF08171", joinpath(data_dir, "4bl0.xml.gz")),
)

bench_ms("Read PDBML", () -> read_file(joinpath(data_dir, "4BL0.xml"), PDBML))

bench_ms("Protein Contact Map", () -> msacontacts(aln_mapping, resdict, col2res))

bench_ms(
    "Mutual Information APC",
    () -> APC!(mapcolpairfreq!(mutual_information, aln_mapping, Frequencies(ContingencyTable(Float64, Val{2}, GappedAlphabet())))),
)

bench_ms("AUC (ROC) for contact prediction, MIp", () -> AUC(MIp, cmap))

# Additional MIToS benchmarks from standalone scripts

bench_read(joinpath(data_dir, "PF00089.stockholm.gz"), Stockholm, "[BENCH] Read compressed Stockholm PF00089: ")
bench_read(joinpath(data_dir, "PF00089.sth"), Stockholm, "[BENCH] Read uncompressed Stockholm PF00089: ")
bench_read(joinpath(data_dir, "PF00089.fasta.gz"), FASTA, "[BENCH] Read compressed FASTA PF00089: ")
bench_read(joinpath(data_dir, "PF00089.fasta"), FASTA, "[BENCH] Read uncompressed FASTA PF00089: ")

bench_read(joinpath(data_dir, "PF16957.stockholm.gz"), Stockholm, "[BENCH] Read compressed Stockholm PF16957: ")
bench_read(joinpath(data_dir, "PF16957.sth"), Stockholm, "[BENCH] Read uncompressed Stockholm PF16957: ")
bench_read(joinpath(data_dir, "PF16957.fasta.gz"), FASTA, "[BENCH] Read compressed FASTA PF16957: ")
bench_read(joinpath(data_dir, "PF16957.fasta"), FASTA, "[BENCH] Read uncompressed FASTA PF16957: ")

msa_long = read_file(joinpath(data_dir, "PF00089_aligned.fasta"), FASTA)
msa_wide = read_file(joinpath(data_dir, "PF16957_aligned.fasta"), FASTA)
@assert size(msa_long) == (EXPECT_SEQ_PF00089, EXPECT_LEN_PF00089) "PF00089 aligned dimensions changed"
@assert size(msa_wide) == (EXPECT_SEQ_PF16957, EXPECT_LEN_PF16957) "PF16957 aligned dimensions changed"

freq_table = Frequencies(ContingencyTable(Float64, Val{1}, UngappedAlphabet()))
_ = mapcolfreq!(shannon_entropy, msa_long, freq_table)
_ = mapcolfreq!(shannon_entropy, msa_wide, freq_table)

println("[BENCH] Shannon entropy PF00089: ", @elapsed mapcolfreq!(shannon_entropy, msa_long, freq_table))
println("[BENCH] Shannon entropy PF16957: ", @elapsed mapcolfreq!(shannon_entropy, msa_wide, freq_table))

_ = percentidentity(msa_long)
_ = percentidentity(msa_wide)

println("[BENCH] Percent Identity PF00089: ", @elapsed percentidentity(msa_long))
println("[BENCH] Percent Identity PF16957: ", @elapsed percentidentity(msa_wide))

buslje_mip(msa) = last(
    buslje09(
        msa;
        samples = 0,           # no random shuffling (for comparability)
        clustering = false,    # no Hobohm clustering
        lambda = 0.0,          # no pseudocount
        maxgap = 1.0,          # keep all columns (match ProDy defaults)
        alphabet = UngappedAlphabet(),
        apc = true,            # report APC-corrected MI (MIp)
    ),
)
_ = buslje_mip(msa_long)
_ = buslje_mip(msa_wide)

println("[BENCH] MIp PF00089: ", @elapsed buslje_mip(msa_long))
println("[BENCH] MIp PF16957: ", @elapsed buslje_mip(msa_wide))
