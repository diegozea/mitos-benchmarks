# mitos-benchmarks
Benchmarking MIToS against various languages and packages

## Pipeline benchmark

|                                                         | MIToS     | Prody           | Bio3D   | BioJulia | BioPython |
|---------------------------------------------------------|-----------|-----------------|---------|----------|-----------|
| Language                                                | Julia     | Python/C        | R       | Julia    | Python    |
| License                                                 | MIT       | MIT             | GPLv2   | MIT      | Biopython |
| **Download Pfam MSA**                                   | Stockholm | Stockholm/FASTA | FASTA   | ✗        | ✗         |
| Read Pfam Stockholm [ms]                                | 1.35      | 0.60            | ✗       | ✗        | 7.32      |
| Read MSA and annotations [ms]                           | 2.07      | ✗               | ✗       | ✗        | ✗         |
| **Read MSA and annotations, generate coordinates [ms]** | 11.16     | ✗               | ✗       | ✗        | ✗         |
| Percent Identity Matrix [ms]                            | 6.87      | 5.71            | 457.00  | ✗        | ✗         |
| **SIFTS residue level mapping [ms]**                    | 40.34     | ✗               | ✗       | ✗        | ✗         |
| **Read PDBML [s]**                                      | 0.5271    | ✗               | ✗       | ✗        | ✗         |
| **Protein Contact Map [ms]**                            | 0.99      | ✗               | 5114.00 | ✗        | ✗         |
| **Mutual Information APC (MIp) [ms]**                   | 23.67     | 13.13           | ✗       | ✗        | ✗         |
| **AUC (ROC) for contact prediction, MIp [ms]**          | 2.65      | ✗               | ✗       | ✗        | ✗         |

### Installations

- [**MIToS**](http://diegozea.github.io/MIToS.jl/): `Pkg.add("MIToS")`
- [**ProDy & Evol**](http://prody.csb.pitt.edu/): `pip install -U ProDy`
- [**Bio3D**](http://thegrantlab.org/bio3d/): `install.packages("bio3d")`
- [**BioJulia**](http://biojulia.github.io/Bio.jl/latest/): `Pkg.add("Bio")`
- [**BioPython**](http://biopython.org/): `pip install numpy` and `pip install biopython`

## MIToS benchmarks

#### MSA module

| hobohmI | |
|---|---|
| 10 | 55.36 μs |
| 20 | 73.48 μs |
| 30 | 124.34 μs |
| 40 | 203.35 μs |
| 50 | 464.75 μs |
| 60 | 797.88 μs |
| 70 | 1.04 ms |
| 80 | 1.07 ms |
| 90 | 802.98 μs |

| output | |
|---|---|
| Stockholm_gzipped | 1.15 ms |
| Stockholm_ungzipped | 539.15 μs |
| FASTA_gzipped | 589.09 μs |
| FASTA_ungzipped | 174.26 μs |

| identity | |
|---|---|
| matrix_Float16 | 5.15 ms |
| matrix_Float32 | 4.81 ms |
| matrix_Float64 | 4.80 ms |
| matrix_BigFloat | 8.43 ms |
| mean | 13.27 ms |

| input | |
|---|---|
| Stockholm_gzipped | 3.75 ms |
| Stockholm_gzipped_mapping | 9.43 ms |
| Stockholm_ungzipped | 1.42 ms |
| Stockholm_ungzipped_mapping | 7.07 ms |
| FASTA_gzipped | 473.97 μs |
| FASTA_gzipped_mapping | 1.88 ms |
| FASTA_ungzipped | 378.86 μs |
| FASTA_ungzipped_mapping | 1.78 ms |

#### PDB module

| output | |
|---|---|
| pdb_PDBFile_gzipped | 79.18 ms |
| pdb_PDBFile_ungzipped | 33.06 ms |
| xml_PDBML_gzipped | 78.37 ms |
| xml_PDBML_ungzipped | 59.85 ms |

| input | |
|---|---|
| pdb_PDBFile_gzipped | 86.12 ms |
| pdb_PDBFile_ungzipped | 48.35 ms |
| xml_PDBML_gzipped | 395.45 ms |
| xml_PDBML_ungzipped | 385.28 ms |

#### Information module

| highlevel | |
|---|---|
| ZBLMIp_ | 59.32 s |
| Buslje09_ | 4.35 s |

| lowlevel | |
|---|---|
| count_col | 1.17 μs |
| count_col_col | 5.05 μs |
| count_col_col_col | 571.06 μs |
| probabilities_col | 3.01 μs |
| probabilities_col_col | 12.20 μs |
| probabilities_col_col_col | 604.04 μs |
| count_col_clusters | 2.27 μs |
| count_col_col_clusters | 7.47 μs |
| count_col_col_col_clusters | 576.96 μs |
| blosum_pseudofrequencies | 614.25 μs |
| probabilities_blosum | 636.79 μs   


| estimateincolumns | |  
|---|---|  
| Entropy_Count_Gapped | 16.66 ms |  
| Entropy_Probability_Gapped | 93.68 μs |  
| MI_Count | 22.15 ms |  
| MI_Probability | 30.68 ms |  
| MI_Count_Gapped | 23.13 ms  |  
| MI_Probability_Gapped | 31.63 ms |  
| Entropy_Count | 17.18 ms |  
| Entropy_Probability | 96.86 μs |  

#### Pfam/pipeline  
  
|   |   |  
|---|---|  
| read_pfam_gzipped | 9.41 ms |  
| getseq2pdb | 21.99 μs |  
| hasresidues | 7.91 μs |  
| msacolumn2pdbresidue_sifts | 28.16 ms |  
| msacolumn2pdbresidue_sifts_gzipped | 46.50 ms |  
| read_PDBML_gzipped | 657.89 ms |  
| residue_list_to_dict | 458.38 μs |  
| msaresidues | 23.85 μs |  
| contact_map | 591.21 μs |  
| buslje09 | 4.15 s |  
| AUC | 2.20 ms |  
