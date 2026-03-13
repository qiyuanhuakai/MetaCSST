# MetaCSST - Metagenomic Complex Sequence Scanning Tool

**Type**: Bioinformatics CLI Tool | **Languages**: C++, Perl | **Domain**: DGR Prediction

## Overview
Tool to predict Diversity-Generating Retroelements (DGRs) in sequenced genomes and metagenomic datasets using Generalized Hidden Markov Models (GHMM).

## Structure
```
./
├── MetaCSSTmain          # Main executable: DGR prediction
├── MetaCSSTsub           # Sub-executable: TR/VR/RT identification
├── *.config              # GHMM configuration files
├── src/                  # Source code (C++ headers + Perl scripts)
├── align/                # Alignment matrices for GHMM training
├── addition/             # Training/test datasets
└── example/              # Example input/output
```

## Entry Points

| Executable | Source | Purpose |
|------------|--------|---------|
| `MetaCSSTmain` | `src/main.cpp` | Full DGR prediction pipeline |
| `MetaCSSTsub` | `src/sub.cpp` | Identify substructures (TR/VR/RT only) |

## Workflow

1. **Build HMM models** from config files (arg.config references TR/VR/RT.config)
2. **Scan sequences** using Viterbi algorithm in `ghmm.h`
3. **Call VRs** via `callVR.pl` (pairs TRs with variant repeats)
4. **Remove duplicates** via `removeRepeat.pl`

## Build Commands

```bash
# Compile main executable
g++ -lpthread src/main.cpp -o MetaCSSTmain

# Compile sub executable
g++ -lpthread src/sub.cpp -o MetaCSSTsub
```

## Key Algorithms

- **Pattern scoring**: Position Weight Matrices (PWM) in `pattern.score[][]`
- **State transitions**: Markov chains with start/end probabilities
- **Viterbi search**: Dynamic programming for optimal path finding
- **Sequence scanning**: Bidirectional (positive/negative strand)

## Conventions

### C++ Code Style
- `S`, `N`, `P`, `M`, `D` - Capitalized constants in `fun.h`
- Classes: `HMM`, `HMM_class`, `SCAN`
- Structs: lowercase with underscores (`pattern`, `sub_hmm`)
- Memory: Manual `calloc/free` (no smart pointers)

### Config Format
```ini
[motif]
motif=align/CLASS.align
cov=0.9       # Coverage cutoff
len=7         # Min pattern length
score=0.3     # State score cutoff
ratio=0.9     # Training score ratio
gap=20        # Max gap between states
```

### Data Flow
1. Input: FASTA format (preprocess with `chomp.pl`)
2. Output: GTF-like tabular format
3. Intermediate: Split files in `tmp/` directories

## Threading Model
- pthread-based multi-threading
- Input file splitting via `split()` in `fun.h`
- Thread-safe: Each thread writes to separate output files

## Critical Code Paths

| Component | File | Lines | Complexity |
|-----------|------|-------|------------|
| HMM scan | `ghmm.h:scanSeqSingle()` | ~150 | High - Viterbi algorithm |
| Class builder | `ghmm.h:buildHMM()` | ~250 | High - PWM construction |
| DGR scan | `ghmm.h:SCAN::scanSeq()` | ~150 | High - Multi-state search |
| VR caller | `callVR.pl` | ~427 | Medium - Perl string ops |

## Anti-Patterns (Don't Do)

- **Don't** modify `ghmm.h` constants without understanding memory implications
- **Don't** skip `chomp.pl` preprocessing on non-standard FASTA
- **Don't** ignore return values from `calloc()` - no null checks in current code
- **Don't** change config `ratio` values without re-training

## Debugging Tips

- Check `align.txt` and `score.txt` in output dirs for PWM debugging
- Use single thread (`-thread 1`) for deterministic output
- VR search parameters hardcoded in `callVR.pl` lines 17-21

## Dependencies
- Linux OS (memory: 2G+ for multi-threading)
- Perl 5.8.5+
- gcc 4.1.2+ with pthread support

## Example Usage
```bash
./MetaCSSTmain -build arg.config -in example/hv29.fa -out results -thread 4
perl src/callVR.pl results/raw.gtf example/hv29.fa out_tmp
perl src/removeRepeat.pl out_tmp final.gtf
```
