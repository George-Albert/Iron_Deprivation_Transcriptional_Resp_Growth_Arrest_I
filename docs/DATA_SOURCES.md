# Data Sources and Description

This document provides complete information about all datasets used in this project.

## Table of Contents

1. [Input Data Overview](#input-data-overview)
2. [RNA-seq Data](#rna-seq-data)
3. [Sample Metadata](#sample-metadata)
4. [Data Processing](#data-processing)
5. [Data Organization](#data-organization)

## Input Data Overview

This project analyzes *Mycobacterium tuberculosis* transcriptional responses under two conditions:

| Dataset | Type | Samples | Format | Location |
|---------|------|---------|--------|----------|
| RNA-seq Reads | Illumina sequencing | 12+ | FASTQ | `data/raw/` |
| Gene Expression Counts | Quantified | 12+ | CSV/TSV | `Analyses/inputs/` |
| Sample Metadata | Experimental | 12+ | XLSX | `GEO_submission/` |
| Lipid Profiles | Lipidomics | 12+ | CSV | `Analyses/inputs/` |

## RNA-seq Data

### Sequencing Parameters

- **Platform**: Illumina (HiSeq or NextSeq)
- **Read Type**: Single-end or paired-end
- **Read Length**: 50-150 bp
- **Sequencing Depth**: ~10-50 million reads per sample
- **Reference Genome**: *Mycobacterium tuberculosis* H37Rv

### Sample Design

| Condition | Iron | Phase | Replicates | Total |
|-----------|------|-------|-----------|-------|
| Iron-replete | + | Exponential (EXP) | 3 | 3 |
| Iron-replete | + | Stationary (STAT) | 3 | 3 |
| Iron-depleted | - | Exponential (EXP) | 3 | 3 |
| Iron-depleted | - | Stationary (STAT) | 3 | 3 |
| **Total** | | | | **12** |

### Data Accessibility

**GEO Accession**: Available in published paper  
**Direct Download**: https://www.ncbi.nlm.nih.gov/geo/

## Sample Metadata

### Critical Variables

```
Sample_ID           - Unique sample identifier
Iron_Status         - "Iron_replete" or "Iron_depleted"
Growth_Phase        - "Exponential" or "Stationary"
Replicate           - Biological replicate number (1, 2, or 3)
Culture_Date        - Date of culture preparation
RNA_Extraction_Date - When RNA was extracted
Sequencing_Run      - Which sequencing run
```

### Sample Naming Convention

```
Format: MTB_[IronStatus]_[Phase]_Rep[Number]

Example:
- MTB_Fe+_EXP_Rep1  → Iron-replete, Exponential, Replicate 1
- MTB_Fe-_STAT_Rep2 → Iron-depleted, Stationary, Replicate 2
```

## Data Processing

### Preprocessing Steps (Raw Reads → Counts)

1. **Quality Control (QC)**
   - Tool: FastQC
   - Check: Sequence quality, adapter contamination
   - Threshold: Q score ≥ 30

2. **Adapter/Primer Removal**
   - Tool: Trimmomatic or Cutadapt
   - Remove Illumina adapters
   - Min length after trimming: 50 bp

3. **Alignment to Reference**
   - Tool: Bowtie2, STAR, or Hisat2
   - Reference: *M. tuberculosis* H37Rv genome
   - Parameters: Default with stringent filtering

4. **Read Counting**
   - Tool: featureCounts or HTSeq
   - Feature: Gene-level counts
   - Mode: Union (default)

5. **Normalization** (in analysis scripts)
   - Method: DESeq2 or EdgeR normalization
   - CPM calculation for exploratory analysis

### Quality Control Metrics

Expected metrics (per sample):

```
Total Reads:          10-50 million
Aligned Reads:        >90%
rRNA Contamination:   <5%
Duplication Rate:     <20%
Gene Coverage:        >80% of genes with >0 reads
```

## Data Organization

### Directory Structure for Data

```
data/
├── raw/                          # Original, unmodified data
│   ├── fastq/                    # Raw sequencing reads
│   ├── gene_counts_raw.csv       # Raw read counts
│   ├── sample_metadata.xlsx      # Sample information
│   └── README.md                 # Data dictionary
│
├── processed/                    # Cleaned/normalized data
│   ├── gene_counts_normalized.csv      # Normalized counts
│   ├── rlog_transformed_counts.csv     # rlog-transformed
│   ├── vst_transformed_counts.csv      # VST-transformed
│   ├── gene_annotations.csv            # Gene names/descriptions
│   └── sample_info_final.csv           # Final sample annotations
│
└── README.md                     # Data documentation

Analyses/inputs/                 # Legacy input data location
├── raw_counts.csv               # Gene count matrix
├── metadata.csv                 # Sample metadata
└── lipid_data.csv               # Lipidomics data
```

## Lipidomics Data

**Analysis**: PDIM (Phthiocerol dimycocerosate) detection

**Key Finding**: PDIM absent in stationary phase regardless of iron status

## Data Restrictions & Attribution

- **License**: Data generated for this study; see paper for terms
- **Citation**: Must cite the original publication
- **Reuse**: Permitted for research; check publication license
- **Sharing**: Permitted; please maintain attribution

## Citing Data

If you use these datasets in your own research:

```
Alebouyeh S, Cárdenas-Pestana JA, et al. (2022). Raw sequencing data for 
"Iron deprivation enhances transcriptional responses to in vitro growth 
arrest of Mycobacterium tuberculosis." NCBI Gene Expression Omnibus. 
doi: 10.3389/fmicb.2022.956602
```

---

**For analysis methods**: See [METHODS.md](METHODS.md)  
**For reproducibility guide**: See [REPRODUCIBILITY.md](REPRODUCIBILITY.md)
