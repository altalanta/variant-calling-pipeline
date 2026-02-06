# Variant Calling Pipeline

[![CI](https://github.com/your-org/variant-calling-pipeline/workflows/CI/badge.svg)](https://github.com/your-org/variant-calling-pipeline/actions)

A complete, containerized bioinformatics pipeline for germline variant calling from paired-end FASTQ files to annotated VCF files.

## Overview

This pipeline implements best-practice germline variant calling using GATK4, with comprehensive quality control and reproducible containerized execution. It processes paired-end sequencing data through alignment, duplicate marking, base quality score recalibration, variant calling, and filtering.

## Quick Start

```bash
# Clone the repository
git clone https://github.com/your-org/variant-calling-pipeline.git
cd variant-calling-pipeline

# Build the reference bundle (for testing)
cd references && bash make_ref_bundle.sh && cd ..

# Run the pipeline with test data
nextflow run workflow/main.nf \
    -profile docker \
    --input tests/smoke/samplesheet.csv \
    --outdir results \
    --ref_fasta tests/smoke/ref/ref.fasta
```

## Requirements

- **Nextflow** (≥22.10.0)
- **Docker** (or Singularity)
- **8+ GB RAM** recommended
- **Reference genome** and indices (see [references documentation](docs/references.md))

## Entropy-Augmented Variant Discovery

Beyond the standard GATK best-practice workflow, the pipeline now emits an entropy-aware variant layer designed to surface variants missed by conventional haplotype assembly. Each sample automatically runs:

1. **Complex region detection** — Sliding windows (configurable via `--entropy_window_size`/`--entropy_window_step`) quantify mismatch/indel densities, soft-clip burden, and base/CIGAR entropies to flag difficult loci.
2. **Local read re-assembly** — Reads spanning flagged windows are reassembled with a lightweight De Bruijn graph (`--entropy_kmer`, `--entropy_max_contigs`) to recover complex indels, MNVs, and micro-variants.
3. **Entropy-weighted scoring** — Read-level allele balance, base qualities, MAPQ, entropy z-scores, and assembly support feed a configurable linear score (`--entropy_weight_*`) to classify calls into baseline-supported, entropy-supported novel, or low-confidence categories.
4. **Reporting** — Three VCFs (baseline copy, augmented, novel-only) plus a JSON summary (Ti/Tv, indel spectra, region densities, overlap statistics) are written under `variants/entropy/` for downstream benchmarking.

## Command Examples

### Basic run with test data
```bash
nextflow run workflow/main.nf \
    -profile docker \
    --input tests/smoke/samplesheet.csv \
    --outdir results \
    --ref_fasta tests/smoke/ref/ref.fasta
```

### Full WGS run with BQSR and joint genotyping
```bash
nextflow run workflow/main.nf \
    -profile docker \
    --input samplesheet.csv \
    --outdir results \
    --ref_fasta reference/hg38.fasta \
    --known_sites_vcf reference/dbsnp_151.hg38.vcf.gz \
    --run_bqsr \
    --run_joint
```

### Exome sequencing with target regions
```bash
nextflow run workflow/main.nf \
    -profile docker \
    --input samplesheet.csv \
    --outdir results \
    --ref_fasta reference/hg38.fasta \
    --targets_bed targets/exome_targets.bed \
    --known_sites_vcf reference/dbsnp_151.hg38.vcf.gz \
    --run_bqsr
```

## Input Requirements

### Samplesheet Format

The input samplesheet must be a CSV file with the following columns:

| Column     | Required | Description                    |
|------------|----------|--------------------------------|
| sample_id  | Yes      | Unique sample identifier       |
| read1      | Yes      | Path to R1 FASTQ file         |
| read2      | Yes      | Path to R2 FASTQ file         |
| platform   | No       | Sequencing platform (default: ILLUMINA) |
| library    | No       | Library name (default: sample_id) |
| lane       | No       | Sequencing lane (default: 1)  |

Example:
```csv
sample_id,read1,read2,platform,library,lane
NA12878,data/NA12878_R1.fastq.gz,data/NA12878_R2.fastq.gz,ILLUMINA,WGS_lib1,1
NA24385,data/NA24385_R1.fastq.gz,data/NA24385_R2.fastq.gz,ILLUMINA,WGS_lib2,2
```

### Reference Files

- **Reference FASTA**: BWA-indexed reference genome
- **Known sites VCF** (optional): dbSNP or Mills indels for BQSR
- **Target BED** (optional): Regions to restrict analysis

## Outputs

The pipeline produces the following outputs for each sample:

```
results/
├── {sample_id}/
│   ├── alignment/
│   │   ├── {sample_id}.cram           # Aligned reads
│   │   └── {sample_id}.cram.crai      # CRAM index
│   ├── variants/
│   │   ├── {sample_id}.snv_indel.vcf.gz    # Filtered variants
│   │   ├── {sample_id}.snv_indel.vcf.gz.tbi # VCF index
│   │   └── entropy/
│   │       ├── {sample_id}.entropy.baseline.vcf.gz
│   │       ├── {sample_id}.entropy.augmented_variants.vcf.gz
│   │       ├── {sample_id}.entropy.novel_entropy_variants.vcf.gz
│   │       ├── {sample_id}.entropy.low_confidence_variants.vcf.gz
│   │       └── {sample_id}.entropy.entropy_summary.json
│   └── qc/
│       └── {sample_id}_fastqc.html
├── multiqc/
│   └── multiqc_report.html            # Aggregate QC report
└── run_metadata/
    ├── params.json                    # Run parameters
    └── versions.txt                   # Tool versions
```

## Tool Versions

- **BWA-MEM2**: 2.2.1
- **SAMtools**: 1.17
- **GATK**: 4.4.0.0
- **bcftools**: 1.17
- **fastp**: 0.23.4
- **MultiQC**: 1.15

Full version manifest available in [`containers/tool_versions.lock`](containers/tool_versions.lock).

## Documentation

- [Pipeline Overview](docs/overview.md)
- [Inputs and Outputs](docs/inputs_outputs.md)
- [Reference Requirements](docs/references.md)
- [Benchmarking](docs/benchmarking.md)
- [Troubleshooting](docs/troubleshooting.md)

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{variant_calling_pipeline,
  title = {Variant Calling Pipeline},
  author = {Your Name},
  year = {2024},
  url = {https://github.com/your-org/variant-calling-pipeline}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please read our contributing guidelines and submit pull requests to the `dev` branch.

## Support

For questions and support:
- Open an [issue](https://github.com/your-org/variant-calling-pipeline/issues)
- Check the [troubleshooting guide](docs/troubleshooting.md)
