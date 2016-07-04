using MIToS: MSA, Pfam, PDB, Information

# Download Pfam Stockholm
# msa = downloadpfam('PF08171')

macro output_timings(t,name)
    quote
        @printf "MIToS,%s,%f\n" $name minimum($t)
    end
end

macro timeit(ex,name)
    quote
        t = Float64[]
        for i in 1:6
            e = 1000*(@elapsed $(esc(ex)));
            if i > 1
                # warm up on first iteration
                push!(t, e)
            end
        end
        @output_timings t $name
    end
end

const aln_mapping = read("../data/PF08171.sth", Stockholm, generatemapping=true, useidcoordinates=true);
const col2res = msacolumn2pdbresidue(aln_mapping, "BUB1_YEAST/291-355", "4BL0", "B", "PF08171","../data/4bl0.xml.gz");
const pdb_residues = read("../data/4BL0.xml", PDBML);
const resdict = @residuesdict pdb_residues model "1" chain "B" group "ATOM" residue "*";
const cmap = msacontacts(aln_mapping, resdict, col2res);
const ZMIp, MIp = buslje09(aln_mapping, samples=0, clustering=false, lambda=0.0, usegap=true);

@timeit read("../data/PF08171.sth", Stockholm, Matrix{Residue}) "Read Pfam Stockholm MSA"
@timeit read("../data/PF08171.sth", Stockholm) "Read MSA and annotations"
@timeit read("../data/PF08171.sth", Stockholm, generatemapping=true, useidcoordinates=true) "Read MSA and annotations, generate coordinates"
@timeit percentidentity(aln_mapping) "Percent Identity Matrix"
@timeit msacolumn2pdbresidue(aln_mapping, "BUB1_YEAST/291-355", "4BL0", "B", "PF08171","../data/4bl0.xml.gz") "SIFTS residue level mapping"
@timeit read("../data/4BL0.xml", PDBML) "Read PDBML"
@timeit msacontacts(aln_mapping, resdict, col2res) "Protein Contact Map"
@timeit APC!(estimateincolumns(aln_mapping, ResidueCount{Int,2,true}, MutualInformation{Float64}())) "Mutual Information APC"
@timeit AUC($MIp, $cmap) "AUC (ROC) for contact prediction, MIp"
