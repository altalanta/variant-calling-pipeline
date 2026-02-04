# Inputs and Outputs

## Input Requirements

### 1. Samplesheet (Required)

**Format**: CSV file with header
**Parameter**: `--input`
**Example**: `--input samplesheet.csv`

#### Required Columns:
- `sample_id`: Unique identifier for the sample
- `read1`: Path to forward read FASTQ file (R1)
- `read2`: Path to reverse read FASTQ file (R2)

#### Optional Columns:
- `platform`: Sequencing platform (default: ILLUMINA)
- `library`: Library name (default: sample_id)
- `lane`: Sequencing lane (default: 1)

#### Example Samplesheet:
```csv
sample_id,read1,read2,platform,library,lane
NA12878,/data/NA12878_R1.fastq.gz,/data/NA12878_R2.fastq.gz,ILLUMINA,WGS_lib1,001
NA24385,/data/NA24385_R1.fastq.gz,/data/NA24385_R2.fastq.gz,ILLUMINA,WGS_lib2,002
```

### 2. Reference Genome (Required)

**Parameter**: `--ref_fasta`
**Format**: FASTA file with BWA-MEM2 indices
**Example**: `--ref_fasta /ref/hg38.fa`

#### Required Index Files:
- `reference.fa.fai` (samtools index)
- `reference.dict` (GATK dictionary)
- `reference.fa.0123` (BWA-MEM2 index)
- `reference.fa.amb` (BWA-MEM2 index)
- `reference.fa.ann` (BWA-MEM2 index)
- `reference.fa.bwt.2bit.64` (BWA-MEM2 index)
- `reference.fa.pac` (BWA-MEM2 index)

### 3. Known Sites VCF (Optional)

**Parameter**: `--known_sites_vcf`
**Format**: Compressed VCF with tabix index
**Example**: `--known_sites_vcf /ref/dbsnp_151.hg38.vcf.gz`
**Purpose**: Required for Base Quality Score Recalibration (BQSR)

#### Required Index:
- `known_sites.vcf.gz.tbi` (tabix index)

### 4. Target Regions (Optional)

**Parameter**: `--targets_bed`
**Format**: BED file (0-based coordinates)
**Example**: `--targets_bed /data/exome_targets.bed`
**Purpose**: Restricts analysis to specific genomic regions

#### BED Format Example:
```
chr1    1000    2000    exon_1
chr1    5000    6000    exon_2
chr2    10000   15000   exon_3
```

## Output Structure

### Directory Layout

```
results/
├── {sample_id}/
│   ├── alignment/
│   │   ├── {sample_id}.cram
│   │   ├── {sample_id}.cram.crai
│   │   └── {sample_id}.markdup.metrics.txt
│   ├── variants/
│   │   ├── {sample_id}.snv_indel.vcf.gz
│   │   ├── {sample_id}.snv_indel.vcf.gz.tbi
│   │   ├── {sample_id}.vcf.gz (unfiltered)
│   │   └── {sample_id}.vcf.gz.tbi
│   ├── gvcf/
│   │   ├── {sample_id}.g.vcf.gz
│   │   └── {sample_id}.g.vcf.gz.tbi
│   └── qc/
│       ├── {sample_id}.fastp.html
│       └── {sample_id}.fastp.json
├── multiqc/
│   ├── multiqc_report.html
│   └── multiqc_data/
├── joint_genotyping/ (if --run_joint enabled)
│   ├── joint_genotyped.vcf.gz
│   └── joint_genotyped.vcf.gz.tbi
└── pipeline_info/
    ├── execution_report.html
    ├── execution_timeline.html
    ├── execution_trace.txt
    ├── pipeline_dag.html
    └── software_versions.yml
```

## Output File Descriptions

### Alignment Files

#### `{sample_id}.cram`
- **Format**: CRAM (compressed SAM)
- **Content**: Aligned, sorted, duplicate-marked reads
- **Size**: ~30-50% smaller than BAM format
- **Usage**: Input for downstream analysis, visualization

#### `{sample_id}.cram.crai`
- **Format**: CRAM index
- **Content**: Index for random access to CRAM file
- **Usage**: Required for tools that access CRAM files

#### `{sample_id}.markdup.metrics.txt`
- **Format**: Tab-delimited text
- **Content**: Duplicate marking statistics
- **Metrics**: 
  - Total reads
  - Duplicate reads/pairs
  - Duplication rate
  - Optical duplicates

### Variant Files

#### `{sample_id}.snv_indel.vcf.gz`
- **Format**: Compressed VCF (v4.2)
- **Content**: Filtered SNVs and indels
- **Filtering**: Hard-filtered variants meeting quality thresholds
- **Usage**: Main variant output for analysis

#### `{sample_id}.snv_indel.vcf.gz.tbi`
- **Format**: Tabix index
- **Content**: Index for random access to VCF file
- **Usage**: Required for tools that query VCF regions

#### `{sample_id}.g.vcf.gz` (GVCF)
- **Format**: Genomic VCF (GVCF)
- **Content**: Variant and reference block information
- **Usage**: Input for joint genotyping across multiple samples

### Quality Control Files

#### `{sample_id}.fastp.html`
- **Format**: HTML report
- **Content**: Read quality metrics, filtering statistics
- **Metrics**:
  - Read count before/after filtering
  - Quality score distributions
  - Adapter content
  - GC content
  - Length distributions

#### `{sample_id}.fastp.json`
- **Format**: JSON
- **Content**: Machine-readable QC metrics
- **Usage**: Input for MultiQC aggregation

#### `multiqc_report.html`
- **Format**: HTML report
- **Content**: Aggregated QC metrics across all samples
- **Sections**:
  - General statistics
  - FastP results
  - Alignment statistics
  - Duplicate metrics
  - Variant calling statistics

## File Size Estimates

### WGS (30X coverage, human genome)
- **Input FASTQ**: 60-80 GB per sample
- **CRAM file**: 15-20 GB per sample
- **VCF file**: 100-500 MB per sample
- **GVCF file**: 1-3 GB per sample

### Exome (100X coverage)
- **Input FASTQ**: 6-12 GB per sample
- **CRAM file**: 1-3 GB per sample
- **VCF file**: 10-50 MB per sample
- **GVCF file**: 100-300 MB per sample

### Test Data (Smoke test)
- **Input FASTQ**: <1 MB per sample
- **CRAM file**: <1 MB per sample
- **VCF file**: <1 MB per sample

## Data Formats and Standards

### FASTQ Format
- **Encoding**: Phred+33 quality scores
- **Compression**: gzip compression recommended
- **Naming**: Consistent R1/R2 suffixes for paired reads

### VCF Format
- **Version**: VCF v4.2
- **Compression**: bgzip compression
- **Indexing**: Tabix indexing required
- **Fields**: Standard VCF fields plus GATK annotations

### CRAM Format
- **Reference**: Reference-compressed alignment format
- **Indexing**: CRAI index format
- **Compression**: ~50% reduction vs BAM files

## Parameter Files

### Generated Parameters File
The pipeline automatically generates a parameters file documenting the exact settings used:

**Location**: `pipeline_info/params.json`

**Content**:
```json
{
  "input": "/path/to/samplesheet.csv",
  "ref_fasta": "/path/to/reference.fa",
  "known_sites_vcf": "/path/to/dbsnp.vcf.gz",
  "outdir": "results",
  "run_bqsr": true,
  "run_joint": false,
  "run_fastp": true
}
```

### Software Versions File
All tool versions are recorded for reproducibility:

**Location**: `pipeline_info/software_versions.yml`

**Format**: YAML with tool versions used in each process

## Data Access Patterns

### Sequential Access
- MultiQC processing
- Variant filtering
- Report generation

### Random Access
- IGV visualization (requires indices)
- Region-specific analysis
- Variant lookup

### Streaming Access
- Alignment (BWA-MEM2 → SAMtools)
- Variant calling pipeline

## Storage Recommendations

### Working Directory
- **Size**: 3-5x input data size
- **Type**: Fast local storage (SSD preferred)
- **Cleanup**: Automatically cleaned by Nextflow

### Results Directory
- **Size**: 1-2x input data size
- **Type**: Network storage acceptable
- **Backup**: Recommend backup of final results

### Reference Data
- **Size**: 3-4 GB for human genome + indices
- **Type**: Fast storage preferred (shared across runs)
- **Caching**: Can be shared across multiple pipeline runs