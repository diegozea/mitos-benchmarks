# Benchmarking Julia's MIToS...
**...against various languages and packages.**

These benchmarks now target current toolchains:

- Julia ≥ 1.9 with **MIToS 3.x**
- Python ≥ 3.11 with **ProDy 2.x** and **Biopython 1.8x**
- R ≥ 4.x with the latest **Bio3D** from CRAN

The scripts are kept small and print `[BENCH]` lines so results remain comparable across 
releases. We use as examples the *Pfam PF08171* (208 sequences, 68 columns without 
inserts), *PF00089/PF16957* alignments, and the *PDB 4BL0* (1133 residues, 6408 atoms).  

## Pipeline benchmark

Here we show the number of *seconds* or *milliseconds* that takes the 
**common steps (in bold)** in a pipeline from a Pfam MSA to the calculation of the 
MIp contact prediction performance. We only include in the benchmark the capabilities that 
are directly provide by the packages like a single function or method. MIToS was designed 
to perform this kind of operations, it takes a different approach to other packages. That 
makes difficult the comparison. MIToS is closer to *ProDy/Evol* in terms of capabilities. 
*Prody* is a *Python* package but their parsing, mutual information and other functions 
are written in C. MIToS is completely written in Julia, which has a performance 
approximately between [1 and 2 times C](http://julialang.org/benchmarks/).  

These times are the minimum time that takes 5 executions of the same function in the 
following computer:
```
  OS: Linux (Debian 6.1.94-1-amd64, kernel 6.1.0-22-amd64)
  CPU: 2 x Intel(R) Xeon(R) Silver 4316 @ 2.30GHz (80 threads total)
  Memory: 1.0 TiB
```

|                                                         | MIToS  | ProDy  | Bio3D        | BioJulia | Biopython     |
|---------------------------------------------------------|--------|--------|--------------|----------|---------------|
| Language                                                | Julia  | Python/C | R          | Julia    | Python        |
| License                                                 | MIT    | MIT    | GPLv2        | MIT      | Biopython     |
| **Download Pfam MSA**                                   | Stockholm | Stockholm/FASTA | FASTA | N/A      | N/A           |
| Read Pfam Stockholm [ms]                                | 0.35 ms | 4.08 ms | N/A         | N/A      | 2.28 ms       |
| Read MSA and annotations [ms]                           | 0.58 ms | N/A    | N/A         | N/A      | 2.28 ms       |
| **Read MSA and annotations, generate coordinates [ms]** | 3.13 ms | N/A    | N/A         | N/A      | N/A           |
| Percent Identity Matrix [ms]                            | 1.65 ms | 4.63 ms | 181.00 ms†† | N/A      | 16.71 ms      |
| **SIFTS residue level mapping [s]**                     | 0.02 s  | N/A    | N/A         | N/A      | N/A           |
| **Read PDBML [s]**                                      | 0.25 s  | N/A    | N/A         | N/A      | N/A           |
| **Protein Contact Map [ms]**                            | 0.37 ms | N/A    | 2156.00 ms  | 38.64 ms | N/A           |
| **Mutual Information APC (MIp) [ms]**                   | 4.93 ms | 5.09 ms | N/A        | N/A      | 413.65 ms     |
| **AUC (ROC) for contact prediction, MIp [ms]**          | 0.09 ms | N/A    | N/A         | N/A      | N/A           |

Run everything from the repository root with the conda env activated:

```
conda activate mitos-benchmarks
julia --project=. run_benchmark.jl
```

### Running the MIToS pipeline benchmark

The MIToS pipeline (`MIToS/Pipeline.jl`) benchmarks the end-to-end contact-prediction steps:

1) Activate the conda env (for R/Python tools) and ensure Julia has MIToS and ROCAnalysis installed.
2) From the repository root:

```
conda activate mitos-benchmarks
julia --project=. MIToS/Pipeline.jl
```

This produces CSV-style lines:

```
MIToS,<step>,<time_ms>
```

where `<step>` covers:
- Read Pfam Stockholm MSA
- Read MSA and annotations
- Read MSA and annotations, generate coordinates
- Percent Identity Matrix
- SIFTS residue level mapping
- Read PDBML
- Protein Contact Map
- Mutual Information APC
- AUC (ROC) for contact prediction, MIp

Use these values to compare against the summary table above.

To collect all pipeline timings in one shot:

```
conda activate mitos-benchmarks
julia --project=. run_pipeline.jl
```

This runs the Pipeline scripts for MIToS, ProDy, Bio3D, and Biopython and prints the same CSV-style rows.

> Notes: MIToS runs `buslje09` with `samples=0` (no shuffling), `clustering=false`, 
  `lambda=0`, and `maxgap=1.0` to mirror ProDy’s “all columns, no clustering” behavior. 
  ProDy uses `buildMutinfoMatrix` plus `applyMutinfoCorr(corr="prod")`. Results are 
  comparable as “MIp-style” but not bit-for-bit identical.  
> †† Bio3D PID uses PF08171 here; PF00089 remains too large for `seqidentity` in this 
  environment (PF16957 PID ≈ 0.37 s in the `[BENCH]` output).  
> N/A = capability not available or not benchmarked in this suite.

### Installations

- Julia project for these benchmarks (MIToS, ROCAnalysis, BioStructures, Graphs, MetaGraphs): `cd mitos-benchmarks && julia --project=. -e 'using Pkg; Pkg.instantiate()'`
- Conda (recommended unified env): `conda env create -f environment.yml`
- [**MIToS**](http://diegozea.github.io/MIToS.jl/): `using Pkg; Pkg.add("MIToS")`
- [**ProDy**](http://prody.csb.pitt.edu/): `python -m pip install -U prody`
- [**Biopython**](http://biopython.org/): `python -m pip install -U biopython numpy`
- [**Bio3D**](http://thegrantlab.org/bio3d/): `install.packages("bio3d")`
- [**BioJulia**](http://biojulia.github.io/Bio.jl/latest/): `Pkg.add("Bio")` (used only in legacy scripts)

### Tests

Small correctness checks using a toy alignment are available:

- Python (ProDy + Biopython): `python -m unittest tests/test_python_metrics.py`
- MIToS: `julia MIToS/TestPID.jl`
- Bio3D: `Rscript Bio3D/TestPID.R`
- Biopython PID benchmarking can be capped via `BIOPYTHON_PID_MAXSEQ`; by default the full alignment is used for comparability.

### Optional structural benchmarks (BioStructures + Bio3D)

- BioStructures (Julia): `julia --project=. BioJulia/BioStructuresBench.jl` computes Cα 
  contact maps, distance maps, and an optional contact graph for 4BL0 chain B. The bundled 
  project includes `BioStructures`, `Graphs`, and `MetaGraphs`; the script skips 
  cleanly if these packages are absent.
- Bio3D structure (R): `Rscript Bio3D/Structure.R` reads 4BL0, builds a Cα contact 
  map (8 Å), distance matrix, NMA-based DCCM, and a contact-filtered correlation network. 
  Requires `bio3d` (already in the conda env) and `igraph` (added to `environment.yml`; 
  the network step is skipped if it is missing).
- These optional sections are wired into `run_benchmark.jl` (run with `julia --project=.`) 
  and will emit `[SKIP]` messages rather than failing when dependencies are unavailable.

## MIToS benchmarks

These benchmarks were run on the machine described above (Debian, dual Xeon Silver 4316, 80 threads).  
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
