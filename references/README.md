# Reference Files

This directory contains scripts and documentation for setting up reference genomes and associated files required by the variant calling pipeline.

## Required Reference Files

For a complete analysis, you will need:

1. **Reference FASTA file** (indexed with BWA-MEM2)
2. **Known sites VCF** (for Base Quality Score Recalibration - optional)
3. **Target BED file** (for exome/targeted sequencing - optional)

## Quick Setup for Testing

To generate a small reference bundle for testing and CI:

```bash
bash make_ref_bundle.sh
```

This will create a minimal reference in `../tests/smoke/ref/` suitable for the smoke test.

## Human Genome References

For production analyses with human data, use one of these standard references:

### GRCh38/hg38 (Recommended)
```bash
# Download reference FASTA
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Create BWA-MEM2 index
bwa-mem2 index hg38.fa

# Download known sites for BQSR
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi

# Download Mills and 1000G gold standard indels
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
```

### GRCh37/hg19 (Legacy)
```bash
# Download reference FASTA
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz

# Create BWA-MEM2 index
bwa-mem2 index hg19.fa

# Download known sites
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg19/v0/dbsnp_138.b37.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg19/v0/dbsnp_138.b37.vcf.gz.tbi
```

## File Structure

Your reference directory should look like:

```
references/
├── hg38.fa                              # Reference FASTA
├── hg38.fa.fai                          # FASTA index
├── hg38.fa.0123                         # BWA-MEM2 indices
├── hg38.fa.amb                          # BWA-MEM2 indices  
├── hg38.fa.ann                          # BWA-MEM2 indices
├── hg38.fa.bwt.2bit.64                  # BWA-MEM2 indices
├── hg38.fa.pac                          # BWA-MEM2 indices
├── dbsnp_151.hg38.vcf.gz                # Known sites VCF
├── dbsnp_151.hg38.vcf.gz.tbi            # VCF index
├── mills_1000g.hg38.vcf.gz              # Mills indels
└── mills_1000g.hg38.vcf.gz.tbi          # VCF index
```

## Reference Index Generation

The pipeline expects BWA-MEM2 indices. To generate them:

```bash
# Index reference FASTA
bwa-mem2 index reference.fasta

# Index for GATK (if not present)
samtools faidx reference.fasta

# Create sequence dictionary for GATK
gatk CreateSequenceDictionary -R reference.fasta -O reference.dict
```

## Target Regions (Exome/Panel)

For exome or targeted sequencing, provide a BED file with target regions:

```bash
# Example BED format (0-based coordinates)
chr1    1000    2000    target_1
chr1    5000    6000    target_2
chr2    1000    3000    target_3
```

## Storage and Performance

- **Reference files can be large** (3GB for hg38)
- **Use fast storage** (SSD) for reference files to improve pipeline performance
- **Index files must be in the same directory** as the reference FASTA

## CI/Test Reference

The smoke test uses a minimal synthetic reference created by `make_ref_bundle.sh`. This reference:

- Contains a single small contig (1000 bp)
- Includes synthetic variants for testing
- Has all required indices pre-built
- Is suitable for continuous integration testing