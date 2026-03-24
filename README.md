# Primer Design and Off-Target Analysis Tool

This repository contains a C++ program for designing primers, optimizing PDRs (Primer Design Regions), checking for dimer formation, and performing off-target search for primers. The tool provides flexible modes and parameters for bioinformatics workflows.

## Features

- **PDR Optimization (`pdr`)**: Optimize primer design regions based on risk scoring.
- **Dimer Optimization (`dimer`)**: Minimize undesired primer-primer interactions.
- **Off-target Search (`off`)**: Identify potential off-target binding sites in reference sequences.
- **Customizable Parameters**: Supports tuning of amplicon length, risk thresholds, thermodynamic parameters, and multithreading.

## Installation

### Prerequisites
- C++17 compatible compiler
- `primer3` (required for thermodynamic calculations)
- Standard C++ libraries
- CMake or Makefile for building (optional)

### Build
```bash
git clone <repo_url>
cd primer-design-tools
mkdir build && cd build
cmake ..
make
```

The executable will be generated, e.g., `primer_tool`.

## Usage

```bash
Usage:
  primer_tool <mode> -i <input> -o <output> [options]

Modes (required):
  pdr       Run PDR optimization
  dimer     Run dimer optimization
  off       Run off-target search

I/O (required):
  -i <file> Input file
  -o <file> Output file
  -r <file> Reference file (required for off mode)
```

### PDR Optimization Parameters (`mode = pdr`)
```text
-S <int>     Random seed (default 42)
-Ln <int>    Amplicon max length (default 420)
-Lx <int>    Amplicon min length (default 252)
-Lp <int>    PDR length (default 40)
-Ux <double> Max risk (default 10000)
-Un <double> Min risk (default 0.1)
```

### Dimer Optimization Parameters (`mode = dimer`)
```text
-I <int> Iterations for solver (default 1000)
```

### Off-target Search Parameters (`mode = off`)
```text
-k <int>        K-mer length (default 8)
-H <int>        Hamming threshold (default 2)
-G <double>     Delta G threshold (default 1e8)
-t <int>        Number of threads (default 8)
-C <int>        Chunk size (default 100000000)
-B <int>        Block size (default 134217728)
--mv <double>   Monovalent conc (default 50.0)
--dv <double>   Divalent conc (default 1.5)
--dntp <double> dNTP conc (default 0.6)
--dna <double>  DNA conc (default 200.0)
--temp <double> Temperature in °C (default 37.0)
```

### Examples

```bash
# PDR optimization
./primer_tool pdr -i in.txt -o out.txt

# Dimer optimization
./primer_tool dimer -i in.txt -o out.txt -I 500

# Off-target search
./primer_tool off -i primers.fa -r ref.fa -o hits.txt -t 16
```

## Input File Formats

- **PDR / Dimer Modes**: Tab-delimited text file with primer sequences and relevant data.
- **Off-target Mode**: FASTA format for primers (`-i`) and reference sequences (`-r`).

## Output

- **PDR Mode**: Optimized primer design regions.
- **Dimer Mode**: Optimized primer indices.
- **Off-target Mode**: List of off-target hits with coordinates and scores.

## License

This project is licensed under the MIT License. See `LICENSE` for details.

## References

- Primer3: [http://primer3.org](http://primer3.org)
- Risk-based primer design and dimer detection methods