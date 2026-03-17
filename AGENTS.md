# MetaCSST - Metagenomic Complex Sequence Scanning Tool

**Type**: Bioinformatics CLI Tool | **Languages**: C++, Python | **Domain**: DGR Prediction

## Overview
Tool to predict Diversity-Generating Retroelements (DGRs) in sequenced genomes and metagenomic datasets using Generalized Hidden Markov Models (GHMM).

## Structure
```
./
├── MetaCSSTmain            # Main executable: DGR prediction
├── MetaCSSTsub             # Sub-executable: TR/VR/RT identification
├── config.json/toml/yaml   # Single-file GHMM configurations
├── src/                    # C++20 source + Python post-processing
├── align/                  # Alignment matrices for GHMM training
├── addition/               # Training/test datasets
└── example/                # Example input/output
```

## Entry Points
| Executable | Source | Purpose |
|------------|--------|---------|
| `MetaCSSTmain` | `src/main_modern.cpp` | Full DGR prediction pipeline |
| `MetaCSSTsub` | `src/sub_modern.cpp` | Identify substructures (TR/VR/RT only) |

## Workflow
1. Build HMM models from a single config file (`config.json` / `config.toml` / `config.yaml`)
2. Scan sequences using Viterbi algorithm in `src/ghmm_modern.hpp`
3. Call VRs and remove duplicates via `src/call_vr.py`

## Build Commands
```bash
make modern
make verify
```

## Dependencies
- Linux OS (memory: 2G+ for multi-threading)
- Python 3.9+ with `numpy`, `numba`, `pyfastx`
- g++ with C++20 support
