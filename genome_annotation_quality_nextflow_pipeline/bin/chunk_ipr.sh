#!/bin/bash
# Usage: bash split_fasta.sh <input.fasta> <chunk_size> <output_dir>

INPUT=${1:-./after_461.fasta}
CHUNK_SIZE=${2:-20}
OUTDIR=${3:-/fasta_chunks}

mkdir -p $OUTDIR

# Split into chunks of N sequences (each FASTA record = 2 lines if unwrapped,
# but we use awk to handle multi-line FASTA correctly)
awk -v chunk=$CHUNK_SIZE -v outdir=$OUTDIR '
  /^>/ {
    seq_count++
    if ((seq_count - 1) % chunk == 0) {
      file_idx = int((seq_count - 1) / chunk)
      filename = sprintf("%s/chunk_%04d.fasta", outdir, file_idx)
    }
    print > filename
    next
  }
  { print > filename }
' $INPUT

echo "Split complete. Chunks written to $OUTDIR"
ls $OUTDIR | wc -l
echo "chunk files created"


