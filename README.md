# DPRO - Primer Design and Off-Target Analysis Tool

A comprehensive C++ tool for primer design optimization and off-target screening. DPRO provides an integrated pipeline for primer analysis workflows including design region optimization, dimer checking, and off-target analysis with advanced filtering capabilities.

## Features

- **Integrated Pipeline**: Complete primer analysis from design through off-target screening
- **Advanced Filtering**: Thermodynamic-based primer quality control with detailed removal reporting
- **Multi-threaded Processing**: Optimized for large reference genomes and high-throughput analysis
- **Comprehensive Reporting**: Detailed tables showing filtered primers and summary statistics
- **Flexible Parameters**: Customizable thermodynamic and search parameters

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
./bin/dpro all -i <input> -o <output> [options]
```

### Core Arguments
```bash
-i <file>    Input sequence file (required)
-o <file>    Output directory (required)
-r <file>    Reference genome file (optional)
-h, --help   Show help message
```

## Parameters

### Search Parameters
```bash
-k <int>       Kmer length for indexing (default: 15)
-H <int>       Hamming distance threshold (default: 3)  
-H2 <int>      Alternative Hamming threshold (default: 3)
-t <int>       Number of threads (default: 8)
-C <int>       Processing chunk size (default: 100000)
-B <int>       I/O block size (default: 8192)
```

### Optimization Parameters
```bash
-Ln <int>     Maximum amplicon length (default: 420)
-Lx <int>     Minimum amplicon length (default: 252)  
-Lp <int>     Primer design region length (default: 40)
-Ux <double>  Maximum risk threshold (default: 10000.0)
-Un <double>  Minimum risk threshold (default: 0.1)
-I <int>      Optimization iterations (default: 1000)
-S <int>      Random seed (default: 42)
```

### Thermodynamic Parameters  
```bash
-G <double>     dG threshold for filtering (default: -20000.0)
--mv <double>   Monovalent concentration mM (default: 50.0)
--dv <double>   Divalent concentration mM (default: 1.5) 
--dntp <double> dNTP concentration mM (default: 0.6)
--dna <double>  DNA concentration nM (default: 50.0)
--temp <double> Temperature °C (default: 60.0)
```

## Examples

### Basic Pipeline Analysis
```bash
# Complete pipeline with reference genome screening
./bin/dpro all \
  -i sequences.fasta \
  -o results/ \
  -r reference_genome.fna

# Pipeline without off-target screening
./bin/dpro all \
  -i sequences.fasta \
  -o results/

# Custom parameters for stringent filtering
./bin/dpro all \
  -i sequences.fasta \
  -o results/ \
  -r genome.fna \
  -G -25000.0 \
  -t 16 \
  -H 2
```

### Advanced Configuration
```bash
# High-throughput analysis
./bin/dpro all \
  -i large_dataset.fasta \
  -o analysis_results/ \
  -r combined_genomes.fna \
  -t 32 \
  -C 200000 \
  -B 16384

# Custom thermodynamic conditions
./bin/dpro all \
  -i primers.fasta \
  -o filtered_results/ \
  -r host_genome.fna \
  --temp 65.0 \
  --mv 100.0 \
  --dna 25.0
```

## Input Format

- **Input**: FASTA format with target sequences
- **Reference**: FASTA format reference genomes (optional)

## Output Structure

The pipeline generates comprehensive analysis reports:

### TABLE A1: Filtered Primer Pairs
```
Output  Index  Direction  dG        Sequence                    Reason
------  -----  ---------  --------  --------------------------  --------------
3       1      F          -17132.7  AATTGGCCTTAATTGGCCTT...     Forward primer
3       1      R          -17132.7  CCTTAACCGGAATTGGCCAA...     Reverse primer
```

### TABLE A2: Filtered Standalone Primers  
```
Output  Index  Type  dG        Sequence                    Reason
------  -----  ----  --------  --------------------------  ---------------
2       39     R     -17129.2  ATGCGATCGATCGATCGATC...     Standalone right
7       39     L     -20314.9  AAATTTGGGCCCAAATTTGG...     Standalone left
```

### TABLE B: Output Summary
```
Output  Original (P/L/R)  Removed (P/L/R)  Final (P/L/R)  Status
------  ----------------  ---------------  -------------  ----------
0       10/63/5           0/0/0            10/63/5        No removal
3       10/28/45          2/3/8            8/25/37        Filtered
6       10/19/83          2/3/12           8/16/71        Filtered
```

### Processing Progress
```
Searching reference.fna (threads: 8)
------------------------------------------------------------
   Contigs      Chunks    Candidates    Time(s)  Rate(ch/s)
------------------------------------------------------------
      1250        1000            78        30      2.5M
      2890        2400           156        60      2.8M
------------------------------------------------------------
Found 234 primer candidates
```

## Performance Features

- **Multi-threaded processing** with configurable thread counts
- **Memory-efficient chunked processing** for large genomes  
- **Progress reporting** with processing rates (characters/second)
- **Optimized indexing** for fast sequence searching

## Pipeline Stages

1. **Design Region Optimization** - Risk-based primer selection
2. **Dimer Analysis** - Detection and minimization of primer interactions  
3. **Off-target Screening** - Reference genome analysis with thermodynamic filtering
4. **Quality Control** - Advanced filtering with detailed reporting

If no reference file is provided:
```
No ref file provided, skip off target search
```

## File Output

- Optimized primer sequences
- Quality control reports  
- Off-target analysis results
- Processing statistics and logs

## License

This project is licensed under the MIT License.

## Support

For questions or issues, please open an issue on the GitHub repository.
