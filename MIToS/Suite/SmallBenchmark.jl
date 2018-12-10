using BenchmarkTools
using ROCAnalysis
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
const small_bench_msa = BenchmarkGroup()

## Set up
### Files
const msafile_sth_gz = "../../data/PF08171.stockholm.gz"
const msafile_sth    = "../../data/PF08171.sth"
const msafile_fas_gz = "../../data/PF08171.fasta.gz"
const msafile_fas    = "../../data/PF08171.fasta"
### MSAs

#### Parse benchmarks

##### input
small_bench_msa["input"] = BenchmarkGroup()

for (file,gzipped,format) in [(msafile_sth_gz, "gzipped", MIToS.MSA.Stockholm),
                              (msafile_sth,  "ungzipped", MIToS.MSA.Stockholm),
                              (msafile_fas_gz, "gzipped", MIToS.MSA.FASTA),
                              (msafile_fas,  "ungzipped", MIToS.MSA.FASTA)]

    # Default parser
    small_bench_msa["input"][string(format,"_",gzipped)] = @benchmarkable read($file, $format)::MIToS.MSA.AnnotatedMultipleSequenceAlignment
    if format != FASTA
    	# With mapping
    	small_bench_msa["input"][string(format,"_",gzipped,"_mapping")] = @benchmarkable read($file, $format, generatemapping=true, useidcoordinates=true)::MIToS.MSA.AnnotatedMultipleSequenceAlignment
    end
end

##### output
small_bench_msa["output"] = BenchmarkGroup()

for (file,gzipped,format) in [(msafile_sth_gz, "gzipped", MIToS.MSA.Stockholm),
                              (msafile_sth,  "ungzipped", MIToS.MSA.Stockholm),
                              (msafile_fas_gz, "gzipped", MIToS.MSA.FASTA),
                              (msafile_fas,  "ungzipped", MIToS.MSA.FASTA)]

    outfile = string("./tmp/",split(file,"/")[end])
    msa_to_save = read(file, format)
    small_bench_msa["output"][string(format,"_",gzipped)] = @benchmarkable write($outfile, $msa_to_save, $format)
end

##### Identity
small_bench_msa["identity"] = BenchmarkGroup()

const aln = read(msafile_sth_gz, Stockholm);

for t in (Float16,Float32,Float64,BigFloat)
    small_bench_msa["identity"][string("matrix_",t)] = @benchmarkable percentidentity($aln, $t)
end

small_bench_msa["identity"][string("mean")] = @benchmarkable meanpercentidentity($aln)

##### Clustering
small_bench_msa["hobohmI"] = BenchmarkGroup()

for pid in 10:10:90
    small_bench_msa["hobohmI"][string(pid)] = @benchmarkable hobohmI($aln, $pid)
end

# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# Information module benchmarks
const small_bench_information = BenchmarkGroup()

##### estimateincolumns
small_bench_information["mapcolfreq!"] = BenchmarkGroup()

const residues = getresidues(aln);

small_bench_information["mapcolfreq!"][string("Entropy_Probability")] = @benchmarkable mapcolfreq!(entropy, $residues, Probabilities(ContingencyTable(Float64,Val{1},UngappedAlphabet())))
small_bench_information["mapcolfreq!"][string("Entropy_Count")] = @benchmarkable mapcolfreq!(entropy, $residues, Counts(ContingencyTable(Float64,Val{1},UngappedAlphabet())))
small_bench_information["mapcolfreq!"][string("MI_Probability")] = @benchmarkable mapcolpairfreq!(mutual_information, $residues, Probabilities(ContingencyTable(Float64,Val{2},UngappedAlphabet())))
small_bench_information["mapcolfreq!"][string("MI_Count")] = @benchmarkable mapcolpairfreq!(mutual_information, $residues, Counts(ContingencyTable(Float64,Val{2},UngappedAlphabet())))

small_bench_information["mapcolfreq!"][string("Entropy_Probability_Gapped")] = @benchmarkable mapcolfreq!(entropy, $residues, Probabilities(ContingencyTable(Float64,Val{1},GappedAlphabet())))
small_bench_information["mapcolfreq!"][string("Entropy_Count_Gapped")] = @benchmarkable mapcolfreq!(entropy, $residues, Counts(ContingencyTable(Float64,Val{1},GappedAlphabet())))
small_bench_information["mapcolfreq!"][string("MI_Probability_Gapped")] = @benchmarkable mapcolpairfreq!(mutual_information, $residues, Probabilities(ContingencyTable(Float64,Val{2},GappedAlphabet())))
small_bench_information["mapcolfreq!"][string("MI_Count_Gapped")] = @benchmarkable mapcolpairfreq!(mutual_information, $residues, Counts(ContingencyTable(Float64,Val{2},GappedAlphabet())))

##### high level
small_bench_information["highlevel"] = BenchmarkGroup()

# default
small_bench_information["highlevel"][string("Buslje09_")] = @benchmarkable buslje09($aln)
small_bench_information["highlevel"][string("ZBLMIp_")] = @benchmarkable BLMI($aln)

##### low level
small_bench_information["lowlevel"] = BenchmarkGroup()

const column_i = residues[:,10];
const column_j = residues[:,20];
const column_k = residues[:,30];
const Pij = probabilities(column_i, column_j);
const Gij = ContingencyTable{Float64,  2, UngappedAlphabet}(UngappedAlphabet());
const nseq_msa = nsequences(aln);
const clusters = hobohmI(aln, 62);

small_bench_information["lowlevel"]["count_col"] = @benchmarkable count($column_i)
small_bench_information["lowlevel"]["count_col_col"] = @benchmarkable count($column_i, $column_j)
small_bench_information["lowlevel"]["count_col_col_col"] = @benchmarkable count($column_i, $column_j, $column_k)
small_bench_information["lowlevel"]["count_col_clusters"] = @benchmarkable count($column_i, weights=$clusters)
small_bench_information["lowlevel"]["count_col_col_clusters"] = @benchmarkable count($column_i, $column_j, weights=$clusters)
small_bench_information["lowlevel"]["count_col_col_col_clusters"] = @benchmarkable count($column_i, $column_j, $column_k, weights=$clusters)
small_bench_information["lowlevel"]["probabilities_col"] = @benchmarkable probabilities($column_i)
small_bench_information["lowlevel"]["probabilities_col_col"] = @benchmarkable probabilities($column_i, $column_j)
small_bench_information["lowlevel"]["probabilities_col_col_col"] = @benchmarkable probabilities($column_i, $column_j, $column_k)

small_bench_information["lowlevel"]["probabilities_blosum"] = @benchmarkable probabilities($column_i, $column_j, pseudofrequencies=BLOSUM_Pseudofrequencies(nsequences(aln), 8.512))

# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# PDB module benchmarks
const small_bench_mitos_pdb = BenchmarkGroup()

## Set up
### Files
const pdb_pdb_gz = "../../data/4BL0.pdb.gz"
const pdb_xml_gz = "../../data/4BL0.xml.gz"
const pdb_pdb    = "../../data/4BL0.pdb"
const pdb_xml    = "../../data/4BL0.xml"
### Residues
const pdb_residues = read(pdb_xml_gz, PDBML);

#### Parse benchmarks

##### input
small_bench_mitos_pdb["input"] = BenchmarkGroup()

for (file,gzipped,label,format) in [(pdb_pdb_gz, "gzipped", "pdb", PDBFile),
                                    (pdb_xml_gz, "gzipped", "xml", PDBML),
                                    (pdb_pdb,  "ungzipped", "pdb", PDBFile),
                                    (pdb_xml,  "ungzipped", "xml", PDBML)]
    # Default parser
    small_bench_mitos_pdb["input"][string(label,"_",format,"_",gzipped)] = @benchmarkable read($file, $format)
end

##### output
small_bench_mitos_pdb["output"] = BenchmarkGroup()

for (file,gzipped,label,format) in [(pdb_pdb_gz, "gzipped", "pdb", PDBFile),
                                    (pdb_xml_gz, "gzipped", "xml", PDBML),
                                    (pdb_pdb,  "ungzipped", "pdb", PDBFile),
                                    (pdb_xml,  "ungzipped", "xml", PDBML)]
    outfile = string("./tmp/",split(file,"/")[end])
    small_bench_mitos_pdb["output"][string(label,"_",format,"_",gzipped)] = @benchmarkable write($outfile, $pdb_residues, $format)
end

# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# Pfam module (pipeline) benchmarks
const small_bench_pfam = BenchmarkGroup()

# Set up
const aln_mapping = read(msafile_sth_gz, Stockholm, generatemapping=true, useidcoordinates=true);
const col2res = msacolumn2pdbresidue(aln_mapping, "BUB1_YEAST/291-355", "4BL0", "B", "PF08171","../../data/4bl0.xml.gz");
const resdict = @residuesdict pdb_residues model "1" chain "B" group "ATOM" residue All;
const cmap = msacontacts(aln_mapping, resdict, col2res);
const ZMIp, MIp = buslje09(aln_mapping);

small_bench_pfam["read_pfam_gzipped"] = @benchmarkable read($msafile_sth_gz, Stockholm, generatemapping=true, useidcoordinates=true)
small_bench_pfam["getseq2pdb"] = @benchmarkable getseq2pdb($aln_mapping)
small_bench_pfam["msacolumn2pdbresidue_sifts"] = @benchmarkable msacolumn2pdbresidue($aln_mapping, "BUB1_YEAST/291-355", "4BL0", "B", "PF08171","../../data/4bl0.xml")
small_bench_pfam["msacolumn2pdbresidue_sifts_gzipped"] = @benchmarkable msacolumn2pdbresidue($aln_mapping, "BUB1_YEAST/291-355", "4BL0", "B", "PF08171","../../data/4bl0.xml.gz")
small_bench_pfam["read_PDBML_gzipped"] = @benchmarkable read("../../data/4BL0.xml.gz", PDBML)
small_bench_pfam["residue_list_to_dict"] = @benchmarkable residuesdict($pdb_residues,"1","B","ATOM",All)
small_bench_pfam["msaresidues"] = @benchmarkable msaresidues($aln_mapping, $resdict, $col2res)
small_bench_pfam["hasresidues"] = @benchmarkable hasresidues($aln_mapping, $col2res)
small_bench_pfam["contact_map"] = @benchmarkable msacontacts($aln_mapping, $resdict, $col2res)
small_bench_pfam["buslje09"] = @benchmarkable buslje09($aln_mapping)
small_bench_pfam["AUC"] = @benchmarkable AUC($ZMIp, $cmap)

# --------------------------------------------------------------------------- #

function SetUp!(;module_msa::Bool=true,module_information::Bool=true,module_pdb::Bool=true,module_pfam::Bool=true)
    if module_msa
        tune!(small_bench_msa)
        JLD.save("small_bench_msa.jld", "small_bench_msa", params(small_bench_msa))
    end
    if module_information
        tune!(small_bench_information)
        JLD.save("small_bench_information.jld", "small_bench_information", params(small_bench_information))
    end
    if module_pdb
        tune!(small_bench_mitos_pdb)
        JLD.save("small_bench_mitos_pdb.jld", "small_bench_mitos_pdb", params(small_bench_mitos_pdb))
    end
    if module_pfam
        tune!(small_bench_pfam)
        JLD.save("small_bench_pfam.jld", "small_bench_pfam", params(small_bench_pfam))
    end
end

function Run!(;module_msa::Bool=true,module_information::Bool=true,module_pdb::Bool=true,module_pfam::Bool=true)
    loadparams!(small_bench_msa, JLD.load("small_bench_msa.jld", "small_bench_msa"), :evals, :samples);
    loadparams!(small_bench_information, JLD.load("small_bench_information.jld", "small_bench_information"), :evals, :samples);
    loadparams!(small_bench_mitos_pdb, JLD.load("small_bench_mitos_pdb.jld", "small_bench_mitos_pdb"), :evals, :samples);
    loadparams!(small_bench_pfam, JLD.load("small_bench_pfam.jld", "small_bench_pfam"), :evals, :samples);
    bench = Dict{String,BenchmarkGroup}()
    if module_msa
        bench["small_bench_msa"] = run(small_bench_msa)
    end
    if module_information
        bench["small_bench_information"] = run(small_bench_information)
    end
    if module_pdb
        bench["small_bench_mitos_pdb"] = run(small_bench_mitos_pdb)
    end
    if module_pfam
        bench["small_bench_pfam"] = run(small_bench_pfam)
    end
    bench
end

#@elapsed include("SmallBenchmark.jl") # 17 s
#@elapsed SetUp!() # 400 s
#@elapsed small_result = Run!() # 338 s
#JLD.save("result_SmallBenchmark.jld", "small_result", small_result);
