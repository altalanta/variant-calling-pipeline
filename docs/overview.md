# Pipeline Overview

## Introduction

The Variant Calling Pipeline is a comprehensive, containerized bioinformatics workflow that implements GATK4 best practices for germline variant calling. It processes paired-end sequencing data from FASTQ files through to filtered, annotated VCF files with comprehensive quality control reporting.

## Workflow Steps

The pipeline implements the following main steps:

### 1. Quality Control and Preprocessing
- **FastP**: Adapter trimming, quality filtering, and read statistics
- **FastQC**: Additional quality metrics (optional)
- Generates per-sample QC reports

### 2. Read Alignment
- **BWA-MEM2**: Alignment to reference genome with read group assignment
- **SAMtools**: Coordinate sorting and indexing
- Supports both WGS and targeted/exome sequencing

### 3. Duplicate Marking
- **SAMtools markdup**: Identifies and marks PCR/optical duplicates
- Outputs CRAM format for space efficiency
- Generates duplicate metrics for QC

### 4. Base Quality Score Recalibration (Optional)
- **GATK BaseRecalibrator**: Models systematic errors in base quality scores
- **GATK ApplyBQSR**: Applies recalibration to improve accuracy
- Requires known variant sites (dbSNP, Mills indels)

### 5. Variant Calling
- **GATK HaplotypeCaller**: Local assembly-based variant calling
- Generates both GVCF and VCF outputs
- Supports interval-based calling for targeted sequencing

### 6. Joint Genotyping (Optional)
- **GATK GenomicsDBImport**: Consolidates GVCFs across samples
- **GATK GenotypeGVCFs**: Joint genotyping for population-scale calling
- Improves sensitivity for rare variants

### 7. Variant Filtering and Annotation
- **GATK VariantFiltration**: Hard filtering based on quality metrics
- **BCFtools norm**: Variant normalization and multiallelic splitting
- Outputs filtered VCF with quality annotations

### 8. Quality Control Reporting
- **MultiQC**: Aggregates QC metrics across all pipeline steps
- Generates comprehensive HTML report

## Pipeline Architecture

```
FASTQ Files
    ↓
FastP (QC + trimming)
    ↓
BWA-MEM2 (alignment)
    ↓
SAMtools (sort + index)
    ↓
SAMtools markdup
    ↓
BQSR (optional)
    ↓
GATK HaplotypeCaller
    ↓
Joint Genotyping (optional)
    ↓
Variant Filtering
    ↓
MultiQC Report
```

## Key Features

### Scalability
- **Single-sample mode**: Process individual samples independently
- **Joint calling mode**: Population-scale analysis with improved sensitivity
- **Resource management**: Configurable CPU/memory requirements

### Flexibility
- **WGS support**: Whole genome sequencing analysis
- **Exome/targeted**: Support for capture-based sequencing
- **Optional steps**: Enable/disable BQSR, joint calling, quality control steps

### Reproducibility
- **Containerized**: All tools packaged in Docker container
- **Version control**: Pinned tool versions for reproducible results
- **Parameter tracking**: Complete parameter logging and version reporting

### Quality Control
- **Comprehensive QC**: FastP, markdup metrics, GATK metrics
- **MultiQC integration**: Unified reporting across all steps
- **Variant QC**: Quality-based filtering with configurable thresholds

## Performance Considerations

### Resource Requirements
- **Minimum**: 8 GB RAM, 4 CPU cores
- **Recommended**: 32 GB RAM, 8+ CPU cores
- **Storage**: ~5x input data size for intermediate files

### Optimization Features
- **CRAM output**: Compressed alignment format reduces storage
- **Parallel processing**: Multi-threaded tools where supported
- **Nextflow caching**: Resumes from failed/interrupted runs

### Runtime Estimates
- **30X WGS (2x150bp)**: 4-8 hours on 16-core system
- **Exome**: 1-2 hours on 8-core system
- **Smoke test**: ~10 minutes on 4-core system

## Output Organization

```
results/
├── sample_id/
│   ├── alignment/           # CRAM files
│   ├── variants/           # VCF files
│   ├── gvcf/              # GVCF files
│   └── qc/                # Sample-specific QC
├── multiqc/               # Aggregate QC report
├── joint_genotyping/      # Joint calling results
└── pipeline_info/         # Execution metadata
```

## Technology Stack

### Core Technologies
- **Nextflow**: Workflow management (DSL2)
- **Docker**: Containerization
- **Conda**: Package management within containers

### Bioinformatics Tools
- **BWA-MEM2 2.2.1**: Next-generation BWA for faster alignment
- **GATK 4.4.0.0**: Variant calling and manipulation
- **SAMtools 1.17**: BAM/CRAM manipulation
- **BCFtools 1.17**: VCF processing
- **FastP 0.23.4**: Read preprocessing
- **MultiQC 1.15**: Quality control reporting

## Best Practices Implementation

The pipeline implements established best practices:

1. **GATK Best Practices**: Follows Broad Institute recommendations
2. **Read group assignment**: Proper sample/library/lane tracking
3. **Duplicate marking**: Handles PCR/optical duplicates
4. **Base recalibration**: Corrects systematic base quality errors
5. **Hard filtering**: Quality-based variant filtering
6. **Variant normalization**: Standardizes variant representation

## Customization Options

### Configurable Parameters
- Quality filtering thresholds
- Resource allocation (CPU/memory)
- Reference genome and known sites
- Target regions for exome/panel data
- Output file formats and compression

### Profile Support
- **docker**: Standard Docker execution
- **singularity**: HPC-compatible containers
- **conda**: Native conda environment
- **local**: Local tool installation

## Integration and Extensibility

### Upstream Integration
- Compatible with standard FASTQ inputs
- Supports various sequencing platforms
- Handles different read lengths and insert sizes

### Downstream Analysis
- Standard VCF output compatible with:
  - Variant annotation (VEP, SnpEff)
  - Population genetics analysis
  - Clinical interpretation tools
  - Visualization tools (IGV, etc.)

### Extension Points
- Modular Nextflow design allows easy addition of:
  - Additional QC metrics
  - Alternative variant callers
  - Structural variant calling
  - Copy number analysis