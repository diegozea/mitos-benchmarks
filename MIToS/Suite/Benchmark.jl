using BenchmarkTools
using JLD
using Test
using MIToS
using MIToS.MSA
using MIToS.Information
using MIToS.PDB
using MIToS.SIFTS
using MIToS.Pfam
using MIToS.Utils

# --------------------------------------------------------------------------- #
# MSA module benchmarks
const msa = BenchmarkGroup()

## Set up
### Files
const msafile_long_sth_gz = "../../data/PF00089.stockholm.gz"
const msafile_long_sth    = "../../data/PF00089.sth"
const msafile_long_fas_gz = "../../data/PF00089.fasta.gz"
const msafile_long_fas    = "../../data/PF00089.fasta"
const msafile_wide_sth_gz = "../../data/PF16957.stockholm.gz"
const msafile_wide_sth    = "../../data/PF16957.sth"
const msafile_wide_fas_gz = "../../data/PF16957.fasta.gz"
const msafile_wide_fas    = "../../data/PF16957.fasta"
### MSAs
const msa_long    = read(msafile_long_sth   , Stockholm)
const msa_wide    = read(msafile_wide_sth   , Stockholm)

#### Parse benchmarks

##### input
msa["input"] = BenchmarkGroup()

for (file,gzipped,shape,format) in [(msafile_long_sth_gz, "gzipped", "long", MIToS.MSA.Stockholm),
                                    (msafile_long_sth,  "ungzipped", "long", MIToS.MSA.Stockholm),
                                    (msafile_long_fas_gz, "gzipped", "long", MIToS.MSA.FASTA),
                                    (msafile_long_fas,  "ungzipped", "long", MIToS.MSA.FASTA),
                                    (msafile_wide_sth_gz, "gzipped", "wide", MIToS.MSA.Stockholm),
                                    (msafile_wide_sth,  "ungzipped", "wide", MIToS.MSA.Stockholm),
                                    (msafile_wide_fas_gz, "gzipped", "wide", MIToS.MSA.FASTA),
                                    (msafile_wide_fas,  "ungzipped", "wide", MIToS.MSA.FASTA)]

    # Default parser
    msa["input"][string(shape,"_",gzipped,"_",format)] = @benchmarkable read($file, $format)::MIToS.MSA.AnnotatedMultipleSequenceAlignment
    # With mapping
    msa["input"][string(shape,"_",gzipped,"_",format,"_mapping")] = @benchmarkable read($file, $format, generatemapping=true, useidcoordinates=true)::MIToS.MSA.AnnotatedMultipleSequenceAlignment
end

##### output
msa["output"] = BenchmarkGroup()

for (file,gzipped,shape,format) in [(msafile_long_sth_gz, "gzipped", "long", MIToS.MSA.Stockholm),
                                    (msafile_long_sth,  "ungzipped", "long", MIToS.MSA.Stockholm),
                                    (msafile_long_fas_gz, "gzipped", "long", MIToS.MSA.FASTA),
                                    (msafile_long_fas,  "ungzipped", "long", MIToS.MSA.FASTA),
                                    (msafile_wide_sth_gz, "gzipped", "wide", MIToS.MSA.Stockholm),
                                    (msafile_wide_sth,  "ungzipped", "wide", MIToS.MSA.Stockholm),
                                    (msafile_wide_fas_gz, "gzipped", "wide", MIToS.MSA.FASTA),
                                    (msafile_wide_fas,  "ungzipped", "wide", MIToS.MSA.FASTA)]
    outfile = string("./tmp/",split(file,"/")[end])
    aln = read(file, format)
    msa["output"][string(shape,"_",format,"_",gzipped)] = @benchmarkable write($outfile, $aln, $format)
end

##### Identity
msa["identity"] = BenchmarkGroup()

for (aln,label) in ((msa_long,"long"), (msa_wide,"wide"))
    for t in (Float16,Float32)
        msa["identity"][string("matrix_",t,"_",label)] = @benchmarkable percentidentity($aln, $t)
    end
    if label == "wide"
        msa["identity"][string("matrix_Float64_wide")] = @benchmarkable percentidentity($aln, Float64)
    end
    msa["identity"][string("mean_",label)] = @benchmarkable meanpercentidentity($aln)
end

##### Clustering
msa["hobohmI"] = BenchmarkGroup()

for (aln,label) in ((msa_long,"long"), (msa_wide,"wide"))
    for pid in 10:10:90
        msa["hobohmI"][string(pid,"_",label)] = @benchmarkable hobohmI($aln, $pid)
    end
end

# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# Information module benchmarks
const information = BenchmarkGroup()

##### estimateincolumns
information["estimateincolumns"] = BenchmarkGroup()

for (aln,label) in ((msa_long,"long"), (msa_wide,"wide"))
    residues = getresidues(aln)
    information["estimateincolumns"][string("Entropy_Probability_",label)] = @benchmarkable estimateincolumns($residues, Int, ResidueProbability{Float64,1,false}, Entropy{Float64}())
    information["estimateincolumns"][string("Entropy_Count_",label)] = @benchmarkable estimateincolumns($residues, ResidueCount{Int,2,false}, Entropy{Float64}())
    information["estimateincolumns"][string("MI_Probability_",label)] = @benchmarkable estimateincolumns($residues, Int, ResidueProbability{Float64,2,false}, MutualInformation{Float64}())
    information["estimateincolumns"][string("MI_Count_",label)] = @benchmarkable estimateincolumns($residues, ResidueCount{Int,2,false}, MutualInformation{Float64}())
    information["estimateincolumns"][string("Entropy_Probability_Gapped_",label)] = @benchmarkable estimateincolumns($residues, Int, ResidueProbability{Float64,1,true}, Entropy{Float64}())
    information["estimateincolumns"][string("Entropy_Count_Gapped_",label)] = @benchmarkable estimateincolumns($residues, ResidueCount{Int,2,true}, Entropy{Float64}())
    information["estimateincolumns"][string("MI_Probability_Gapped_",label)] = @benchmarkable estimateincolumns($residues, Int, ResidueProbability{Float64,2,true}, MutualInformation{Float64}())
    information["estimateincolumns"][string("MI_Count_Gapped_",label)] = @benchmarkable estimateincolumns($residues, ResidueCount{Int,2,true}, MutualInformation{Float64}())
end

##### high level
information["highlevel"] = BenchmarkGroup()

for (aln,label) in ((msa_long,"long"), (msa_wide,"wide"))
    # default
    information["highlevel"][string("Buslje09_",label)] = @benchmarkable buslje09($aln)
    information["highlevel"][string("ZBLMIp_",label)] = @benchmarkable BLMI($aln)
end

##### low level
information["lowlevel"] = BenchmarkGroup()

const residues = getresidues(msa_long)
const column_i = residues[:,10]
const column_j = residues[:,20]
const column_k = residues[:,30]
const Pij = probabilities(column_i, column_j)
const Gij = ContingencyTable{Float64,  2, UngappedAlphabet}(UngappedAlphabet())
const nseq_msa_long = nsequences(msa_long)
const clusters_long = hobohmI(msa_long, 62)

information["lowlevel"]["count_col"] = @benchmarkable count($column_i)
information["lowlevel"]["count_col_col"] = @benchmarkable count($column_i, $column_j)
information["lowlevel"]["count_col_col_col"] = @benchmarkable count($column_i, $column_j, $column_k)
information["lowlevel"]["count_col_clusters"] = @benchmarkable count($column_i, weight=$clusters_long)
information["lowlevel"]["count_col_col_clusters"] = @benchmarkable count($column_i, $column_j, weight=$clusters_long)
information["lowlevel"]["count_col_col_col_clusters"] = @benchmarkable count($column_i, $column_j, $column_k, weight=$clusters_long)
information["lowlevel"]["probabilities_col"] = @benchmarkable probabilities($column_i)
information["lowlevel"]["probabilities_col_col"] = @benchmarkable probabilities($column_i, $column_j)
information["lowlevel"]["probabilities_col_col_col"] = @benchmarkable probabilities($column_i, $column_j, $column_k)

information["lowlevel"]["blosum_pseudofrequencies"] = @benchmarkable blosum_pseudofrequencies!($Gij, $Pij)
information["lowlevel"]["probabilities_blosum"] = @benchmarkable probabilities($column_i, $column_j, pseudofrequencies=BLOSUM_Pseudofrequencies(nsequences(msa), 8.512))

# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# PDB module benchmarks
const mitos_pdb = BenchmarkGroup()

## Set up
### Files
const pdb_pdb_gz = "../../data/2XWB.pdb.gz"
const pdb_xml_gz = "../../data/2XWB.xml.gz"
const pdb_pdb    = "../../data/2XWB.pdb"
const pdb_xml    = "../../data/2XWB.xml"
### Residues
const pdb_residues = read(pdb_xml_gz, PDBML)

#### Parse benchmarks

##### input
mitos_pdb["input"] = BenchmarkGroup()

for (file,gzipped,label,format) in [(pdb_pdb_gz, "gzipped", "pdb", PDBFile),
                                    (pdb_xml_gz, "gzipped", "xml", PDBML),
                                    (pdb_pdb,  "ungzipped", "pdb", PDBFile),
                                    (pdb_xml,  "ungzipped", "xml", PDBML)]
    # Default parser
    mitos_pdb["input"][string(label,"_",format,"_",gzipped)] = @benchmarkable read($file, $format)
end

##### output
mitos_pdb["output"] = BenchmarkGroup()

for (file,gzipped,label,format) in [(pdb_pdb_gz, "gzipped", "pdb", PDBFile),
                                    (pdb_xml_gz, "gzipped", "xml", PDBML),
                                    (pdb_pdb,  "ungzipped", "pdb", PDBFile),
                                    (pdb_xml,  "ungzipped", "xml", PDBML)]
    outfile = string("./tmp/",split(file,"/")[end])
    mitos_pdb["output"][string(label,"_",format,"_",gzipped)] = @benchmarkable write($outfile, $pdb_residues, $format)
end

# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# Pfam module (pipeline) benchmarks
const pfam = BenchmarkGroup()

# Set up
const aln = read(msafile_long_sth_gz, Stockholm, generatemapping=true, useidcoordinates=true)
const col2res = msacolumn2pdbresidue(aln, "CFAB_HUMAN/481-752", "2XWB", "F", "PF00089","../../data/2xwb.xml")
const pdb = read("../../data/2XWB.xml.gz", PDBML)
const resdict = @residuesdict pdb model "1" chain "F" group "ATOM" residue All
const cmap = msacontacts(aln, resdict, col2res)
const ZMIp, MIp = buslje09(aln)

pfam["read_pfam_gzipped"] = @benchmarkable read($msafile_long_sth_gz, Stockholm, generatemapping=true, useidcoordinates=true)
pfam["getseq2pdb"] = @benchmarkable getseq2pdb($aln)
pfam["msacolumn2pdbresidue_sifts"] = @benchmarkable msacolumn2pdbresidue($aln, "CFAB_HUMAN/481-752", "2XWB", "F", "PF00089","../../data/2xwb.xml")
pfam["msacolumn2pdbresidue_sifts_gzipped"] = @benchmarkable msacolumn2pdbresidue($aln, "CFAB_HUMAN/481-752", "2XWB", "F", "PF00089","../../data/2xwb.xml.gz")
pfam["read_PDBML_gzipped"] = @benchmarkable read("../../data/2XWB.xml.gz", PDBML)
pfam["residue_list_to_dict"] = @benchmarkable residuesdict($pdb,"1","F","ATOM",All)
pfam["msaresidues"] = @benchmarkable msaresidues($aln, $resdict, $col2res)
pfam["hasresidues"] = @benchmarkable hasresidues($aln, $col2res)
pfam["contact_map"] = @benchmarkable msacontacts($aln, $resdict, $col2res)
pfam["buslje09"] = @benchmarkable buslje09($aln)
pfam["AUC"] = @benchmarkable AUC($ZMIp, $cmap)

# --------------------------------------------------------------------------- #

function SetUp!(;module_msa::Bool=true,module_information::Bool=true,module_pdb::Bool=true,module_pfam::Bool=true)
    if module_msa
        tune!(msa)
        JLD.save("msa.jld", "msa", params(msa))
    end
    if module_information
        tune!(information)
        JLD.save("information.jld", "information", params(information))
    end
    if module_pdb
        tune!(mitos_pdb)
        JLD.save("mitos_pdb.jld", "mitos_pdb", params(mitos_pdb))
    end
    if module_pfam
        tune!(pfam)
        JLD.save("pfam.jld", "pfam", params(pfam))
    end
end

function Run!(;module_msa::Bool=true,module_information::Bool=true,module_pdb::Bool=true,module_pfam::Bool=true)
    loadparams!(msa, JLD.load("msa.jld", "msa"), :evals, :samples);
    loadparams!(information, JLD.load("information.jld", "information"), :evals, :samples);
    loadparams!(mitos_pdb, JLD.load("mitos_pdb.jld", "mitos_pdb"), :evals, :samples);
    loadparams!(pfam, JLD.load("pfam.jld", "pfam"), :evals, :samples);
    bench = Dict{ASCIIString,BenchmarkGroup}()
    if module_msa
        bench["msa"] = run(msa)
    end
    if module_information
        bench["information"] = run(information)
    end
    if module_pdb
        bench["mitos_pdb"] = run(mitos_pdb)
    end
    if module_pfam
        bench["pfam"] = run(pfam)
    end
    bench
end

# @elapsed include("Benchmark.jl") # 486 s
# @elapsed SetUp!() # 5907 s
# @elapsed result = Run!() # 5642 s
# JLD.save("result_Benchmark.jld", "result", result);
