# Benchmarking Julia's MIToS... 
**...against various languages and packages.**

The following benchmarks are conducted using the *master* branch of MIToS (using *Julia 0.4.5*).
We use as an example the *Pfam PF08171* (208 sequences, 68 columns without inserts) and the *PDB 4BL0* (1133 residues, 6408 atoms).  

A more detailed benchmark of the **MIToS' PDB module** can be found in the [**pdb-benchmarks**](https://github.com/jgreener64/pdb-benchmarks) repository.  

## Pipeline benchmark

Here we show the number of *seconds* or *milliseconds* that takes the **common steps (in bold)** in a pipeline from a Pfam MSA to the calculation of the MIp contact prediction performance. We only include in the benchmark the capabilities that are directly provide by the packages like a single function or method. MIToS was designed to perform this kind of operations, it takes a different approach to other packages. That makes difficult the comparison. MIToS is closer to *ProDy/Evol* in terms of capabilities. *Prody* is a *Python* package but their parsing, mutual information and other functions are written in C. MIToS is completely written in Julia, which has a performance approximately between [1 and 2 times C](http://julialang.org/benchmarks/).  

These times are the minimum time that takes 5 executions of the same function in the following computer:
```
  System: Linux (x86_64-linux-gnu)
  CPU: Intel(R) Core(TM) i5 CPU         750  @ 2.67GHz
```

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

These benchmarks were run on a computer with Ubuntu’s and one i7 CPU.  
The following times are useful to choose the fastest method signatures.  
This benchmark will be used to improve MIToS performance in the near future.  

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
