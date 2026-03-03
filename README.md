# FauxFastQ
A tool for generating simulated sequencing errors of an input DNA sequence.

## Description

This program generates simulated sequencing data by introducing controlled edits (substitutions, insertions, deletions) at every position in a reference DNA sequence. It's designed for testing bioinformatics pipelines and generating synthetic datasets with known ground truth.

## Edit Modes

- **Deletion** (`-d`): Removes the specified number of bases at each position
- **Insertion** (`-i`): Inserts random bases at each position
- **Replacement** (`-r`): Substitutes bases with random nucleotides
- **Multiple modes**: Combine with `--append` to include all specified edit types in one file

| Argument | Type | Description |
|----------|------|-------------|
| `reference_fa` | path | Input reference FASTA file |
| `output_forward` | path | Output file for forward strand reads |
| `output_reverse` | path | Output file for reverse strand reads |
| `-n, --edit_length` | int | Length of edits in base pairs (default: 1) |
| `--unidirectional` | flag | Simulate unidirectional sequencing only |
| `--reference_forward` | path | Reference FASTQ file (forward) for quality score sampling |
| `--reference_reverse` | path | Reference FASTQ file (reverse) for quality score sampling |
| `-s, --seed` | int | Random seed for reproducible edits |
| `-d, --deletion` | flag | Enable deletion mode |
| `-i, --insertion` | flag | Enable insertion mode |
| `-r, --replacement` | flag | Enable replacement mode |
| `--append` | flag | Combine multiple edit modes in one file |
| `--sample_size` | int | Number of qualities to sample from reference files (default: 500) |
| `--use_quality` | bool | Output FASTQ (True) or FASTA (False) format (default: True) |
| `-l, --read_length` | int | Manually set read length (overrides reference read length) |

## Examples

See 'examples' directory for example usage.

## Requirements

- Python 3.6+
- pyfastx
