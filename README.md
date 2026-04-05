# DPRO - Primer Design and Off-Target Analysis Tool

A comprehensive C++ tool for bluetongue virus (BTV) primer design, optimization, and off-target screening. DPRO provides three specialized modes for bioinformatics workflows: PDR optimization, dimer checking, and off-target analysis with advanced filtering capabilities.

## Features

- **PDR Optimization (`pdr`)**: Risk-based primer design region optimization with customizable amplicon constraints
- **Dimer Analysis (`dimer`)**: Detection and minimization of primer-primer interactions  
- **Off-target Screening (`off`)**: High-performance off-target search with thermodynamic filtering
- **Advanced Filtering**: dG-based primer quality control with detailed removal reporting
- **Multi-threaded Processing**: Optimized for large reference genomes (livestock host screening)
- **Comprehensive Reporting**: Detailed tables showing filtered primers and summary statistics

## Installation

### Prerequisites
- C++17 compatible compiler (GCC 7+, Clang 5+)
- Standard C++ libraries (`<thread>`, `<atomic>`, `<unordered_map>`)
- POSIX-compliant system (Linux/macOS)

### Build
```bash
git clone <repo_url>
cd dpro
make
# or
g++ -std=c++17 -O3 -pthread -o bin/dpro src/*.cpp
```

The executable will be generated as `bin/dpro`.

## Usage

```bash
./bin/dpro <mode> -i <input> -o <output> [options]
```

### Modes
- `pdr` - PDR optimization mode
- `dimer` - Dimer analysis mode  
- `off` - Off-target search mode
- `all` - Run complete pipeline (PDR → dimer → off-target)

### Core Arguments
```bash
-i <file>    Input file (required)
-o <file>    Output file/directory (required)
-r <file>    Reference genome file (optional for off mode)
-h, --help   Show help message
```

## Mode-Specific Parameters

### PDR Optimization (`pdr`)
```bash
-Ln <int>     Amplicon max length (default: 420)
-Lx <int>     Amplicon min length (default: 252)  
-Lp <int>     PDR length (default: 40)
-Ux <double>  Max risk threshold (default: 10000.0)
-Un <double>  Min risk threshold (default: 0.1)
-S <int>      Random seed (default: 42)
```

### Dimer Analysis (`dimer`)
```bash
-I <int>   Optimization iterations (default: 1000)
-S <int>   Random seed (default: 42)
```

### Off-target Search (`off`)
```bash
# Search Parameters
-k <int>       Kmer length (default: 15)
-H <int>       Hamming threshold (default: 3)  
-H2 <int>      Hamming threshold odd (default: 3)
-t <int>       Number of threads (default: 8)
-C <int>       Chunk size (default: 100000)
-B <int>       I/O block size (default: 8192)

# Thermodynamic Parameters  
-G <double>     dG threshold for filtering (default: -20000.0)
--mv <double>   Monovalent concentration mM (default: 50.0)
--dv <double>   Divalent concentration mM (default: 1.5) 
--dntp <double> dNTP concentration mM (default: 0.6)
--dna <double>  DNA concentration nM (default: 50.0)
--temp <double> Temperature °C (default: 60.0)
```

## Examples

### BTV Primer Screening Workflow
```bash
# Complete pipeline for BTV segment 1
./bin/dpro all \
  -i ./scripts/aligned_sequences/segment_1.aligned.fasta \
  -o ./bin/tmp \
  -r ./scripts/ref/livestock_combined.fna

# Off-target screening only
./bin/dpro off \
  -i primers.fasta \
  -o results/ \
  -r livestock_genomes.fna \
  -G -25000.0 \
  -t 16

# PDR optimization with custom parameters
./bin/dpro pdr \
  -i sequences.fasta \
  -o optimized_primers.txt \
  -Ln 300 \
  -Lx 150 \
  -Ux 5000.0
```

## Input Formats

### PDR/Dimer Modes
- **Input**: FASTA format with target sequences
- **Output**: Optimized primer sets with coordinates and scores

### Off-target Mode  
- **Input**: FASTA format primer sequences
- **Reference**: FASTA format reference genomes
- **Output**: Filtered primers with off-target analysis

## Output Structure

### Off-target Search Output
The tool generates comprehensive filtering reports:

**TABLE A1: Filtered Primer Pairs**
```
Output  Index  Direction  dG        Sequence                    Reason
------  -----  ---------  --------  --------------------------  --------------
3       1      F          -17132.7  AATTGGCCTTAATTGGCCTT...     Forward primer
3       1      R          -17132.7  CCTTAACCGGAATTGGCCAA...     Reverse primer
```

**TABLE A2: Filtered Standalone Primers**  
```
Output  Index  Type  dG        Sequence                    Reason
------  -----  ----  --------  --------------------------  ---------------
2       39     R     -17129.2  ATGCGATCGATCGATCGATC...     Standalone right
7       39     L     -20314.9  AAATTTGGGCCCAAATTTGG...     Standalone left
```

**TABLE B: Output Summary**
```
Output  Original (P/L/R)  Removed (P/L/R)  Final (P/L/R)  Status
------  ----------------  ---------------  -------------  ----------
0       10/63/5           0/0/0            10/63/5        No removal
3       10/28/45          2/3/8            8/25/37        Filtered
```

## Performance Features

- **Multi-threaded processing** with configurable thread counts
- **Memory-efficient chunked processing** for large genomes  
- **Progress reporting** with processing rates (characters/second)
- **Optimized for livestock genome screening** (cattle, sheep, goat)

## Reference Genome Support

The tool is optimized for screening against livestock host genomes:
- **Cattle**: ARS-UCD1.2 (~2.7GB)
- **Sheep**: ARS-UI_Ramb_v2.0 (~2.9GB)  
- **Goat**: ARS1 (~2.9GB)

If no reference file is provided for off-target mode:
```
No ref file provided, skip off target search
```

## License

This project is licensed under the MIT License.

## Citation

If you use DPRO in your research, please cite:
Scaling Variant-Aware Multiplex Primer Design
Yunheng Han, Christina Boucher
bioRxiv 2026.02.03.703607; doi: https://doi.org/10.64898/2026.02.03.703607

## Support

For questions or issues, please open an issue on the GitHub repository or contact the development team.
