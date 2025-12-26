# Remap: Assembly-Mediated Long-Read Sequence Re-mapping

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://www.python.org/)
[![Bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)

**Remap** is a projection-based re-mapping framework designed to improve the alignment consistency of long-read sequencing data (PacBio HiFi and ONT). 

Standard read-to-reference alignment often suffers from ambiguity in repetitive regions (e.g., segmental duplications) and breakpoint jitter in structural variants (SVs). Remap addresses this by using **assembly contigs as a bridge**: reads are first mapped to their local assembly, and then projected onto the reference genome. This approach significantly improves the precision and recall of downstream SV and SNP/Indel calling.

## Key Features

* **Assembly-Mediated Projection:** Resolves multi-mapping ambiguity by utilizing the longer context of assembled contigs.
* **Exact Anchor Consistency:** Uses strict one-to-one mapping anchors to prevent breakpoint drifting.
* **High Precision in Complex Regions:** Specifically optimized for low-mappability regions, segmental duplications, and extreme GC regions.
* **Versatile Support:** Works with PacBio HiFi, ONT, and R10+ (HQLR) data.
* **Compatible Output:** Generates standard BAM files compatible with existing callers like Sniffles2, SVIM, and DeepVariant.

## Method Overview

Instead of aligning reads directly to the reference ($R \to G$), Remap decomposes the problem into two easier steps:

1.  **Read-to-Assembly ($R \to A$):** Reads are aligned to de novo assembled contigs.
2.  **Assembly-to-Reference ($A \to G$):** Contigs are aligned to the reference genome.
3.  **Coordinate Projection:** Valid anchors are projected from $R$ to $G$ via $A$, and gaps are filled using local alignment.

## Installation

### Prerequisites
* Python >= 3.10
* [minimap2](https://github.com/lh3/minimap2)
* [samtools](http://www.htslib.org/)
* [hifiasm](https://github.com/chhylp123/hifiasm) (if raw reads are provided for assembly)

### Option 1: Using Conda (Recommended)

We provide a  environment file.

```bash
git clone https://github.com/micahvista/Remap.git
cd Remap
conda env create -f remap_environment.yml
conda activate remap_env
pip install .
```
# Usage
Remap requires an input read set, a reference genome, and a working directory. It can perform de novo assembly internally or use a pre-computed assembly.

# Quick Start
Run full pipeline (Assembly + Remap) for PacBio HiFi data:
```bash
remap -i "reads/*.fastq.gz" \
      -r reference.fasta \
      -w ./workdir \
      -o ./output/sample_name \
      -d hifi \
      -t 16
```

Run full pipeline (Assembly + Remap) for ONT R10.4 data:
```bash
remap -i "reads/*.fastq.gz" \
      -r reference.fasta \
      -w ./workdir \
      -o ./output/sample_name \
      -d hqlr \
      -t 16
```

| Argument | Required | Description |
| :--- | :---: | :--- |
| `-i`, `--inputreadpath` | **Yes** | Path to input reads (supports wildcards, e.g., `*.fq.gz`). |
| `-r`, `--refpath` | **Yes** | Path to the reference genome (FASTA). |
| `-w`, `--workdir` | **Yes** | Working directory for intermediate files. |
| `-o`, `--outputprefix` | **Yes** | Prefix for the final output BAM file. |
| `-d`, `--datatype` | **Yes** | Data type: `ont` (R9.x), `hqlr` (R10.x), or `hifi` (PacBio HiFi). |
| `-t`, `--threads` | No | Number of threads (default: 8). |
| `-a`, `--asmdata` | No | Path to existing assembly (FASTA or GFA). If omitted, `hifiasm` is run. |
| `--localasm` | No | Enable local assembly mode for reads containing SVs. |
| `--asm2ref` | No | Path to existing Assembly-to-Reference BAM. |
| `--read2asm` | No | Path to existing Read-to-Assembly BAM. |
| `--rutg` | No | Use raw unitig graph for assembly (default: False). |

# Output
The pipeline generates a standard BAM file containing the re-mapped reads.
## Merged Output: 
  {outputprefix}all.bam (Contains both re-mapped and rescued reads).
## Tags: 
Includes standard SAM tags plus custom tags for alignment consistency.
# Method Overview
### Read-to-Assembly ($R \to A$): 
Reads are aligned to high-quality contigs to resolve local ambiguities.
### Assembly-to-Reference ($A \to G$): 
Contigs are aligned to the reference, spanning repetitive regions with uniqueness.
### Projection: 
Coordinates are projected mathematically from $R$ to $G$ via $A$ using exact anchors.
### Refinement: 
Gaps between anchors are filled using local alignment (Smith-Waterman) to ensure base-level precision.
# Performance Benchmarks
Based on the GIAB HG002 T2T benchmark:
Structural Variants: Remap increases F1 scores by +0.04 ~ +0.05 over direct mapping.
Segmental Duplications: Significant recall improvement ($\Delta R \approx +0.15$ for ONT) in these difficult regions.
Small Indels: Achieves an F1 score of 0.82 (vs 0.77 for direct mapping) for HiFi data22.
# Citation
If you use Remap in your research, please cite:

This project is licensed under the MIT License.
