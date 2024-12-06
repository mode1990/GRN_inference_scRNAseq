#!/bin/bash

# Script to run pySCENIC analysis steps
# Usage: ./run_pyscenic.sh <input_loom> <output_dir>

set -e

INPUT_LOOM=$1
OUTPUT_DIR=$2

if [ -z "$INPUT_LOOM" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: $0 <input_loom> <output_dir>"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

echo "Starting pySCENIC analysis..."

# Step 1: GRN inference
echo "Running GRN inference..."
docker run -v $(pwd):/data aertslab/pyscenic:0.12.1 pyscenic grn \
    /data/"$INPUT_LOOM" \
    /data/allTFs_hg38.txt \
    -o /data/"$OUTPUT_DIR"/adjacencies.tsv \
    --num_workers 4

# Step 2: Regulon prediction
echo "Running regulon prediction..."
docker run -v $(pwd):/data aertslab/pyscenic:0.12.1 pyscenic ctx \
    /data/"$OUTPUT_DIR"/adjacencies.tsv \
    /data/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
    --annotations_fname /data/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname /data/"$INPUT_LOOM" \
    --output /data/"$OUTPUT_DIR"/reg.csv \
    --mask_dropouts \
    --num_workers 4 \
    > /data/"$OUTPUT_DIR"/pyscenic_ctx_stdout.txt

# Step 3: AUCell scoring
echo "Running AUCell scoring..."
docker run -v $(pwd):/data aertslab/pyscenic:0.12.1 pyscenic aucell \
    /data/"$INPUT_LOOM" \
    /data/"$OUTPUT_DIR"/reg.csv \
    --output /data/"$OUTPUT_DIR"/adata_processed_output.loom \
    --num_workers 4 \
    > /data/"$OUTPUT_DIR"/pyscenic_aucell_stdout.txt

echo "pySCENIC analysis complete!"
