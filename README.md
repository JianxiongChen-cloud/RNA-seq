# RNA-seq
RNA-Seq Processing Pipeline
# RNA-seq Analysis of NPM1 Knockdown in HCT116 Cells

## Project Overview

This repository contains the RNA-seq analysis pipeline for investigating transcriptomic changes induced by **NPM1 knockdown** in the human colorectal cancer cell line **HCT116**.

The workflow includes:

- Raw read quality control
- Genome alignment
- Gene quantification
- Differential expression analysis
- Functional enrichment analysis
- Publication-quality visualization

---

## Experimental Design

### Cell Line
- HCT116

### Groups
- **NC**: Negative control shRNA
- **KD**: NPM1 knockdown

### Biological Replicates

| Group | Samples |
|------|---------|
| NC | NC1, NC2, NC3 |
| KD | KD1, KD2, KD3 |

During quality assessment, **KD2** was identified as a potential outlier based on mapping metrics and PCA clustering. Therefore, the main downstream analysis was performed using:

- KD1
- KD3
- NC1
- NC2
- NC3

---

## Computational Workflow

```text
FASTQ
→ FastQC / MultiQC
→ STAR alignment
→ featureCounts
→ DESeq2
→ GO enrichment
→ Hallmark GSEA
→ Visualization
