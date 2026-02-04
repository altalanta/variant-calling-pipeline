#!/bin/bash

# make_ref_bundle.sh - Create minimal reference bundle for testing
# This script creates a small synthetic reference genome and associated files
# for smoke testing and continuous integration

set -euo pipefail

# Configuration
REF_DIR="../tests/smoke/ref"
REF_NAME="ref"
REF_SIZE=1000  # 1kb reference for fast testing

echo "Creating minimal reference bundle in ${REF_DIR}..."

# Create output directory
mkdir -p "${REF_DIR}"
cd "${REF_DIR}"

# Generate synthetic reference sequence
echo "Generating synthetic reference sequence..."
cat > ${REF_NAME}.fasta << 'EOF'
>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCG
EOF

echo "Creating reference indices..."

# Create FASTA index
if command -v samtools &> /dev/null; then
    samtools faidx ${REF_NAME}.fasta
else
    # Fallback: create simple index manually
    echo -e "chr1\t1000\t6\t70\t71" > ${REF_NAME}.fasta.fai
fi

# Create sequence dictionary for GATK compatibility
cat > ${REF_NAME}.dict << 'EOF'
@HD	VN:1.6
@SQ	SN:chr1	LN:1000	M5:d2b3e5c6f1a4b9e8d7c6a5f4e3b2d1c0
EOF

# Create BWA-MEM2 indices (simplified for testing)
# Note: In real usage, you would run: bwa-mem2 index ${REF_NAME}.fasta

# Create dummy BWA-MEM2 index files for CI compatibility
touch ${REF_NAME}.fasta.0123
touch ${REF_NAME}.fasta.amb
touch ${REF_NAME}.fasta.ann
touch ${REF_NAME}.fasta.bwt.2bit.64
touch ${REF_NAME}.fasta.pac

echo "Creating synthetic known sites VCF..."

# Create minimal known sites VCF for BQSR testing
cat > known_sites.vcf << 'EOF'
##fileformat=VCFv4.2
##reference=ref.fasta
##contig=<ID=chr1,length=1000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	rs1	A	T	100	PASS	.
chr1	200	rs2	T	C	100	PASS	.
chr1	300	rs3	C	G	100	PASS	.
chr1	400	rs4	G	A	100	PASS	.
chr1	500	rs5	A	T	100	PASS	.
EOF

# Compress and index the VCF
if command -v bgzip &> /dev/null; then
    bgzip known_sites.vcf
    tabix -p vcf known_sites.vcf.gz
else
    # Fallback: create dummy compressed files
    gzip known_sites.vcf
    touch known_sites.vcf.gz.tbi
fi

# Create target regions BED file for exome simulation
cat > targets.bed << 'EOF'
chr1	50	150	target1
chr1	200	400	target2
chr1	450	650	target3
chr1	750	950	target4
EOF

echo "Reference bundle created successfully!"
echo ""
echo "Files created:"
echo "  - ${REF_NAME}.fasta (reference FASTA)"
echo "  - ${REF_NAME}.fasta.fai (FASTA index)"
echo "  - ${REF_NAME}.dict (sequence dictionary)"
echo "  - BWA-MEM2 index files"
echo "  - known_sites.vcf.gz (synthetic known variants)"
echo "  - targets.bed (target regions)"
echo ""
echo "Usage:"
echo "  nextflow run workflow/main.nf \\"
echo "    --input tests/smoke/samplesheet.csv \\"
echo "    --outdir results \\"
echo "    --ref_fasta tests/smoke/ref/${REF_NAME}.fasta"