# Benchmarking and Validation

## Overview

This document describes how to benchmark the variant calling pipeline against known truth sets and how to evaluate its performance using standard metrics.

## GIAB Benchmarking with hap.py

### Setup

The Genome in a Bottle (GIAB) consortium provides high-confidence variant calls for several reference samples. We recommend using NA12878 (HG001) for initial validation.

#### Download GIAB Truth Sets

```bash
# Create benchmarking directory
mkdir -p benchmarks/giab
cd benchmarks/giab

# Download NA12878 (HG001) truth set for GRCh38
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi

# Download high-confidence regions
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed

# Download reference FASTQ files (if needed)
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
```

### Run Pipeline on GIAB Data

```bash
# Create samplesheet for NA12878
cat > na12878_samplesheet.csv << EOF
sample_id,read1,read2
HG001,benchmarks/giab/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz,benchmarks/giab/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
EOF

# Run pipeline
nextflow run workflow/main.nf \
  -profile docker \
  --input na12878_samplesheet.csv \
  --outdir results_giab \
  --ref_fasta reference/hg38.fa \
  --known_sites_vcf reference/dbsnp_151.hg38.vcf.gz \
  --run_bqsr
```

### Install hap.py

```bash
# Using conda
conda install -c bioconda hap.py

# Or using Docker
docker pull pkrusche/hap.py:latest
```

### Run hap.py Comparison

```bash
# Basic comparison
hap.py \
  benchmarks/giab/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  results_giab/HG001/variants/HG001.snv_indel.vcf.gz \
  -f benchmarks/giab/HG001_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  -r reference/hg38.fa \
  -o benchmark_results/hg001_comparison \
  --engine=vcfeval \
  --pass-only

# Stratified analysis
hap.py \
  benchmarks/giab/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  results_giab/HG001/variants/HG001.snv_indel.vcf.gz \
  -f benchmarks/giab/HG001_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  -r reference/hg38.fa \
  -o benchmark_results/hg001_stratified \
  --stratification stratification/hg38_stratification.tsv \
  --engine=vcfeval \
  --pass-only
```

### Expected Performance Metrics

#### SNP Performance (GIAB NA12878)
- **Sensitivity (Recall)**: >99.5%
- **Precision (PPV)**: >99.8%
- **F1-score**: >99.6%
- **False Positive Rate**: <0.1 per Mb

#### Indel Performance (GIAB NA12878)
- **Sensitivity (Recall)**: >98.5%
- **Precision (PPV)**: >99.0%
- **F1-score**: >98.7%
- **False Positive Rate**: <0.5 per Mb

## Custom Benchmarking

### Performance Test Suite

Create a test suite for different scenarios:

```bash
# Create test configurations
mkdir -p benchmarks/custom

# WGS test (30X)
cat > benchmarks/custom/wgs_30x.config << EOF
params {
  input = 'benchmarks/data/wgs_30x_samplesheet.csv'
  expected_sensitivity_snp = 0.995
  expected_precision_snp = 0.998
  expected_sensitivity_indel = 0.985
  expected_precision_indel = 0.990
}
EOF

# Exome test (100X)
cat > benchmarks/custom/exome_100x.config << EOF
params {
  input = 'benchmarks/data/exome_100x_samplesheet.csv'
  expected_sensitivity_snp = 0.998
  expected_precision_snp = 0.999
  expected_sensitivity_indel = 0.990
  expected_precision_indel = 0.995
}
EOF
```

### Automated Benchmarking Workflow

```bash
#!/bin/bash
# benchmark_pipeline.sh

set -euo pipefail

TRUTH_VCF=$1
QUERY_VCF=$2
CONFIDENT_BED=$3
REFERENCE=$4
OUTPUT_PREFIX=$5

echo "Running hap.py comparison..."

# Run comparison
hap.py \
  "$TRUTH_VCF" \
  "$QUERY_VCF" \
  -f "$CONFIDENT_BED" \
  -r "$REFERENCE" \
  -o "$OUTPUT_PREFIX" \
  --engine=vcfeval \
  --pass-only \
  --threads 8

# Parse results
python << EOF
import pandas as pd
import json

# Read hap.py summary
summary = pd.read_csv('${OUTPUT_PREFIX}.summary.csv')

# Extract key metrics
snp_metrics = summary[summary['Type'] == 'SNP'].iloc[0]
indel_metrics = summary[summary['Type'] == 'INDEL'].iloc[0]

results = {
    'snp_sensitivity': snp_metrics['METRIC.Recall'],
    'snp_precision': snp_metrics['METRIC.Precision'],
    'snp_f1': snp_metrics['METRIC.F1_Score'],
    'indel_sensitivity': indel_metrics['METRIC.Recall'],
    'indel_precision': indel_metrics['METRIC.Precision'],
    'indel_f1': indel_metrics['METRIC.F1_Score']
}

# Save results
with open('${OUTPUT_PREFIX}.metrics.json', 'w') as f:
    json.dump(results, f, indent=2)

print("Benchmarking Results:")
print(f"SNP Sensitivity: {results['snp_sensitivity']:.4f}")
print(f"SNP Precision: {results['snp_precision']:.4f}")
print(f"INDEL Sensitivity: {results['indel_sensitivity']:.4f}")
print(f"INDEL Precision: {results['indel_precision']:.4f}")
EOF
```

## Performance Monitoring

### Runtime Benchmarking

Track pipeline performance across different configurations:

```bash
# Create runtime benchmarking script
cat > benchmarks/runtime_benchmark.sh << 'EOF'
#!/bin/bash

# Runtime benchmarking for different data sizes
declare -A datasets=(
  ["small"]="1M reads"
  ["medium"]="10M reads" 
  ["large"]="100M reads"
  ["xlarge"]="1B reads"
)

for size in "${!datasets[@]}"; do
  echo "Benchmarking ${size} dataset (${datasets[$size]})"
  
  start_time=$(date +%s)
  
  nextflow run workflow/main.nf \
    -profile docker \
    --input "benchmarks/data/${size}_samplesheet.csv" \
    --outdir "results_bench_${size}" \
    --ref_fasta reference/hg38.fa
  
  end_time=$(date +%s)
  runtime=$((end_time - start_time))
  
  echo "${size},${datasets[$size]},${runtime}" >> runtime_results.csv
done
EOF

chmod +x benchmarks/runtime_benchmark.sh
```

### Resource Usage Monitoring

```bash
# Monitor resource usage during execution
cat > benchmarks/monitor_resources.sh << 'EOF'
#!/bin/bash

# Start monitoring
mkdir -p monitoring
nohup top -b -d 5 > monitoring/cpu_usage.log &
TOP_PID=$!

nohup iostat -x 5 > monitoring/io_usage.log &
IOSTAT_PID=$!

# Run pipeline
nextflow run workflow/main.nf \
  -profile docker \
  --input samplesheet.csv \
  --outdir results \
  --ref_fasta reference/hg38.fa

# Stop monitoring
kill $TOP_PID $IOSTAT_PID

# Parse monitoring data
python scripts/parse_monitoring.py monitoring/ > resource_report.txt
EOF
```

## Quality Metrics

### Variant Quality Distributions

```bash
# Extract quality metrics from VCF
bcftools query -f '%QUAL\t%DP\t%GQ\n' results/sample/variants/sample.snv_indel.vcf.gz > variant_metrics.txt

# Plot quality distributions
python << EOF
import matplotlib.pyplot as plt
import pandas as pd

# Load metrics
df = pd.read_csv('variant_metrics.txt', sep='\t', names=['QUAL', 'DP', 'GQ'])

# Create plots
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Quality score distribution
axes[0].hist(df['QUAL'], bins=50, alpha=0.7)
axes[0].set_xlabel('Variant Quality (QUAL)')
axes[0].set_ylabel('Count')
axes[0].set_title('Variant Quality Distribution')

# Depth distribution
axes[1].hist(df['DP'], bins=50, alpha=0.7)
axes[1].set_xlabel('Read Depth (DP)')
axes[1].set_ylabel('Count')
axes[1].set_title('Read Depth Distribution')

# Genotype quality distribution
axes[2].hist(df['GQ'], bins=50, alpha=0.7)
axes[2].set_xlabel('Genotype Quality (GQ)')
axes[2].set_ylabel('Count')
axes[2].set_title('Genotype Quality Distribution')

plt.tight_layout()
plt.savefig('variant_quality_metrics.png', dpi=300)
print("Quality metrics plot saved to variant_quality_metrics.png")
EOF
```

### Ti/Tv Ratio Analysis

```bash
# Calculate transition/transversion ratio
bcftools stats results/sample/variants/sample.snv_indel.vcf.gz | grep "TSTV"

# Expected Ti/Tv ratios:
# - Whole genome: ~2.0-2.1
# - Exome: ~2.8-3.0
```

## Comparison with Other Callers

### Multi-caller Comparison

```bash
# Run multiple variant callers for comparison
# (This would require implementing additional caller modules)

# Compare results using hap.py
for caller in gatk freebayes deepvariant; do
  hap.py \
    truth.vcf.gz \
    results_${caller}/sample.vcf.gz \
    -f confident.bed \
    -r reference.fa \
    -o comparison_${caller} \
    --engine=vcfeval
done

# Generate comparison report
python scripts/compare_callers.py comparison_*.summary.csv > caller_comparison.html
```

## Continuous Integration Benchmarking

### CI Performance Tests

Add to `.github/workflows/benchmark.yml`:

```yaml
name: Performance Benchmark

on:
  push:
    branches: [main]
  schedule:
    - cron: '0 2 * * 0'  # Weekly

jobs:
  benchmark:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      
      - name: Run GIAB benchmark
        run: |
          # Download test data
          ./benchmarks/download_giab_test.sh
          
          # Run pipeline
          nextflow run workflow/main.nf -profile docker --input giab_test.csv
          
          # Run hap.py comparison
          ./benchmarks/run_happly_comparison.sh
          
          # Check performance thresholds
          python benchmarks/check_performance.py happly_results.json
      
      - name: Upload benchmark results
        uses: actions/upload-artifact@v3
        with:
          name: benchmark-results
          path: benchmark_results/
```

## Performance Optimization

### Identifying Bottlenecks

```bash
# Use Nextflow execution reports to identify slow processes
nextflow run workflow/main.nf -with-report execution_report.html

# Analyze trace file for resource usage
python scripts/analyze_trace.py trace.txt > performance_analysis.txt
```

### Resource Tuning

```bash
# Test different resource allocations
for cpus in 4 8 16 32; do
  for memory in 16 32 64 128; do
    echo "Testing ${cpus} CPUs, ${memory}GB RAM"
    
    nextflow run workflow/main.nf \
      --max_cpus $cpus \
      --max_memory "${memory}.GB" \
      -profile docker \
      --input test_samplesheet.csv
      
    # Record runtime
    echo "${cpus},${memory},$(get_runtime)" >> resource_optimization.csv
  done
done
```

## Expected Benchmarking Results

### Performance Targets

#### Accuracy (vs GIAB truth sets)
- **SNP F1-score**: >99.6%
- **Indel F1-score**: >98.5%
- **Ti/Tv ratio**: 2.0-2.1 (WGS), 2.8-3.0 (Exome)

#### Runtime (30X WGS, 16 cores)
- **Total runtime**: 4-8 hours
- **Alignment**: 2-4 hours
- **Variant calling**: 1-2 hours
- **Post-processing**: <30 minutes

#### Resource Usage
- **Peak memory**: 16-32 GB
- **Storage**: 5x input size (temporary)
- **Network I/O**: Minimal (local processing)

### Regression Testing

Run benchmarks on every release to ensure:
1. No decrease in accuracy metrics
2. No significant increase in runtime
3. Consistent resource usage
4. Reproducible results across platforms