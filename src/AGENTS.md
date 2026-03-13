# MetaCSST Source Code

**Domain**: GHMM implementation for DGR prediction

## Files

| File | Lines | Purpose |
|------|-------|---------|
| `main.cpp` | 251 | Main executable - DGR prediction pipeline |
| `sub.cpp` | 196 | Sub-executable - TR/VR/RT identification only |
| `ghmm.h` | 1036 | **Core**: HMM classes, Viterbi algorithm, scanning |
| `fun.h` | 242 | Utility functions (string ops, sorting, file split) |
| `callVR.pl` | 427 | Perl: Search VRs from TR predictions |
| `removeRepeat.pl` | 40 | Perl: Deduplicate TR-VR pairs |
| `chomp.pl` | ~20 | FASTA preprocessing |
| `callORF.pl` | ~50 | ORF calling utility |
| `coden.txt` | - | Codon table for ORF calling |

## Class Hierarchy

```
HMM              # Single GHMM model
├── pattern      # PWM state data
├── box          # HMM state definition
├── sub_hmm      # Sub-structure (TR/VR/RT)
└── OUT          # Scan results

HMM_class        # Cluster of HMM models
└── HMM[] hmm    # Multiple motif models

SCAN             # Main DGR scanner
└── HMM_class[3] # TR, VR, RT clusters
```

## Key Functions

| Function | Location | Purpose |
|----------|----------|---------|
| `buildHMM()` | `ghmm.h:392` | Build PWM from alignment |
| `scanSeqSingle()` | `ghmm.h:197` | Viterbi scan (positive strand) |
| `scanSeqFull()` | `ghmm.h:361` | Bidirectional scan |
| `split()` | `fun.h:84` | File splitting for threading |
| `complementary()` | `fun.h:108` | DNA reverse complement |

## Constants (from `fun.h`)

```cpp
#define S 100000   // Max search space
#define N 20000000 // Max sequence length
#define P 1000     // Max pattern length
#define M 1000     // Max pattern count
#define D 1000     // Max path length
```

## Memory Model
- All memory manually allocated via `calloc()`
- No smart pointers or RAII
- Caller responsible for `free()`
- Pattern: Allocate in init, use, caller frees result

## Threading
- pthread-based (`-lpthread`)
- File split → parallel scan → result merge
- Each thread writes to separate temp file

## Config Parsing
Config files parsed in `HMM_class::init()`:
```cpp
// Format: name=value
// Sections: [motif] starts new model
// Multiple [motif] sections = multiple HMM models
```
