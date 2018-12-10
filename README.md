# Benchmarking Julia's MIToS... 
**...against various languages and packages.**

The following benchmarks are conducted using MIToS v2.3.1 and Julia 1.0.2.
We use as an example the *Pfam PF08171* (208 sequences, 68 columns without inserts) and the *PDB 4BL0* (1133 residues, 6408 atoms).  

A more detailed benchmark of the **MIToS' PDB module** can be found in the [**pdb-benchmarks**](https://github.com/jgreener64/pdb-benchmarks) repository.  

## Pipeline benchmark

Here we show the number of *seconds* or *milliseconds* that takes the **common steps (in bold)** in a pipeline from a Pfam MSA to the calculation of the MIp contact prediction performance. We only include in the benchmark the capabilities that are directly provide by the packages like a single function or method. MIToS was designed to perform this kind of operations, it takes a different approach to other packages. That makes difficult the comparison. MIToS is closer to *ProDy/Evol* in terms of capabilities. *Prody* is a *Python* package but their parsing, mutual information and other functions are written in C. MIToS is completely written in Julia, which has a performance approximately between [1 and 2 times C](http://julialang.org/benchmarks/).  

These times are the minimum time that takes 5 executions of the same function in the following computer:
```
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Core(TM) i7-8550U CPU @ 1.80GHz
  Memory: 15.387805938720703 GB
```

|                                                         | MIToS     | Prody           | Bio3D   | BioJulia | BioPython |
|---------------------------------------------------------|-----------|-----------------|---------|----------|-----------|
| Language                                                | Julia     | Python/C        | R       | Julia    | Python    |
| License                                                 | MIT       | MIT             | GPLv2   | MIT      | Biopython |
| **Download Pfam MSA**                                   | Stockholm | Stockholm/FASTA | FASTA   | ✗        | ✗         |
| Read Pfam Stockholm [ms]                                | 0.46      | 0.70            | ✗       | ✗        | 3.99      |
| Read MSA and annotations [ms]                           | 0.83      | ✗               | ✗       | ✗        | ✗         |
| **Read MSA and annotations, generate coordinates [ms]** | 9.45      | ✗               | ✗       | ✗        | ✗         |
| Percent Identity Matrix [ms]                            | 2.06      | 4.38            | 267.00  | ✗        | ✗         |
| **SIFTS residue level mapping [s]**                     | 0.03      | ✗               | ✗       | ✗        | ✗         |
| **Read PDBML [s]**                                      | 0.29      | ✗               | ✗       | NA       | ✗         |
| **Protein Contact Map [ms]**                            | 0.63      | ✗               | 5026.00 | NA       | ✗         |
| **Mutual Information APC (MIp) [ms]**                   | 4.77      | 8.50            | ✗       | ✗        | ✗         |
| **AUC (ROC) for contact prediction, MIp [ms]**          | 0.79      | ✗               | ✗       | ✗        | ✗         |

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
| 40 | 201.389 μs |  
| 80 | 876.749 μs |  
| 20 | 76.812 μs |  
| 10 | 56.748 μs |  
| 70 | 865.224 μs |  
| 90 | 694.931 μs |  
| 50 | 412.437 μs |  
| 30 | 130.826 μs |  
| 60 | 670.980 μs |  
  
| output | |  
|---|---|  
| Stockholm_ungzipped | 2.102 ms |  
| FASTA_gzipped | 2.677 ms |  
| FASTA_ungzipped | 1.852 ms |  
| Stockholm_gzipped | 3.443 ms |  
  
| identity | |  
|---|---|  
| matrix_Float64 | 1.899 ms |  
| matrix_BigFloat | 3.034 ms |  
| matrix_Float16 | 2.223 ms |  
| mean | 434.556 ms |  
| matrix_Float32 | 1.918 ms |  
  
| input | |  
|---|---|  
| Stockholm_ungzipped | 835.052 μs |  
| FASTA_gzipped | 914.909 μs |  
| Stockholm_gzipped_mapping | 12.116 ms |  
| FASTA_ungzipped | 181.769 μs |  
| Stockholm_ungzipped_mapping | 11.562 ms |  
| Stockholm_gzipped | 1.315 ms |  
  
#### PDB module  
  
| output | |  
|---|---|  
| xml_PDBML_gzipped | 108.046 ms |  
| pdb_PDBFile_ungzipped | 18.712 ms |  
| xml_PDBML_ungzipped | 63.930 ms |  
| pdb_PDBFile_gzipped | 61.506 ms |  
  
| input | |  
|---|---|  
| xml_PDBML_gzipped | 264.028 ms |  
| pdb_PDBFile_ungzipped | 12.043 ms |  
| xml_PDBML_ungzipped | 245.476 ms |  
| pdb_PDBFile_gzipped | 19.502 ms |  
  
#### Information module  
  
| highlevel | |  
|---|---|  
| ZBLMIp_ | 13.228 s |  
| Buslje09_ | 1.379 s |  
  
| mapcolfreq! | |  
|---|---|  
| Entropy_Count_Gapped | 51.292 μs |  
| MI_Count | 4.056 ms |  
| Entropy_Probability | 49.799 μs |  
| MI_Probability | 5.639 ms |  
| Entropy_Count | 48.138 μs |  
| Entropy_Probability_Gapped | 50.423 μs |  
| MI_Count_Gapped | 4.533 ms |  
| MI_Probability_Gapped | 6.312 ms |  
  
| lowlevel | |  
|---|---|  
| probabilities_col_col | 2.997 μs |  
| count_col | 604.861 ns |  
| probabilities_col_col_col | 64.766 μs |  
| count_col_col | 2.219 μs |  
| count_col_col_clusters | 2.304 μs |  
| count_col_col_col | 39.995 μs |  
| probabilities_blosum | 103.877 μs |  
| count_col_clusters | 593.821 ns |  
| probabilities_col | 632.840 ns |  
| count_col_col_col_clusters | 39.601 μs |  

#### Pfam/pipeline  
  
|   |   |  
|---|---|  
| read_pfam_gzipped | 11.937 ms |  
| contact_map | 413.881 μs |  
| getseq2pdb | 7.306 μs |  
| msacolumn2pdbresidue_sifts | 35.733 ms |  
| residue_list_to_dict | 158.207 μs |  
| hasresidues | 17.626 μs |  
| AUC | 788.969 μs |  
| msaresidues | 21.869 μs |  
| read_PDBML_gzipped | 261.228 ms |  
| buslje09 | 1.380 s |  
| msacolumn2pdbresidue_sifts_gzipped | 37.308 ms |  

