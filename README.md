# GRN_inference_scRNAseq
GRN inference by pyscenic for scRNAseq data 

## Overview

This pipeline performs end-to-end analysis of single-cell RNA sequencing data, including preprocessing, gene regulatory network inference using pySCENIC, and visualization of results.

## Pipeline Components

### 1. Preprocessing (`src/preprocessing.py`)

The preprocessing module handles initial data cleanup and preparation:

#### Input Data Requirements
- Raw counts matrix in h5ad format
- Cell type annotations in `adata.obs['BroadCellType']`
- Gene names in `adata.var_names`

#### Processing Steps
1. Cell type filtering (optional)
2. Total count normalization (target sum: 1e4)
3. Log transformation
4. Highly variable gene selection
5. Optional transcription factor filtering

#### Output Files
- `preprocessed_data.h5ad`: Processed AnnData object
- `preprocessing_report.html`: Quality control metrics

### 2. pySCENIC Analysis (`scripts/run_pyscenic.sh`)

Performs regulatory network inference in three steps:

#### GRN Inference
- Input: Expression matrix (loom format)
- Output: `adjacencies.tsv`
- Metrics: Correlation and importance scores for TF-target relationships

#### Regulon Prediction
- Input: `adjacencies.tsv`
- Output: `reg.csv`
- Contains: Predicted regulons with target genes and motif information

#### AUCell Scoring
- Input: Expression matrix and regulons
- Output: `adata_processed_output.loom`
- Metrics: Regulon activity scores per cell

### 3. Visualization (`src/visualization.py`)

Generates multiple visualization types:

#### UMAP Embeddings
- Based on regulon activity scores
- Colored by metadata (e.g., mutation status)
- Output: `umap_plot.pdf`

#### Regulon Activity Heatmap
- Groups cells by specified metadata
- Shows top variable regulons
- Output: `regulon_heatmap.pdf`

## Output Files

### Directory Structure
```
output_dir/
├── preprocessing/
│   ├── preprocessed_data.h5ad
│   └── preprocessing_report.html
├── pyscenic/
│   ├── adjacencies.tsv
│   ├── reg.csv
│   └── adata_processed_output.loom
├── visualization/
│   ├── umap_plot.pdf
│   └── regulon_heatmap.pdf
└── logs/
    └── pipeline_TIMESTAMP.log
```

### File Descriptions

#### preprocessed_data.h5ad
- Normalized and filtered expression data
- Highly variable genes
- Cell annotations

#### adjacencies.tsv
Columns:
- TF: Transcription factor name
- Target: Target gene name
- Importance: Regulatory importance score
- Correlation: Expression correlation
- Status: Activation/repression prediction

#### reg.csv
Columns:
- Regulon name
- Target genes
- Motif information
- Enrichment scores

#### adata_processed_output.loom
Contents:
- Raw expression matrix
- Normalized data
- Regulon AUC scores
- Cell metadata

## Usage Examples

### Basic Run
```bash
python main.py --config config.yaml
```

### Cell Type-Specific Analysis
```bash
# Modify config.yaml:
preprocessing:
  cell_type: "DAN"
  n_top_genes: 2000
```

### Custom Visualization
```bash
# Modify config.yaml:
visualization:
  group_by: "Mutation"
  top_n_regulons: 50
```
