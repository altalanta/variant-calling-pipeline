# Troubleshooting Guide

## Common Issues and Solutions

### Installation and Setup Issues

#### 1. Nextflow Installation Problems

**Problem**: Nextflow not found or version incompatible
```bash
nextflow: command not found
# OR
N E X T F L O W  ~  version 20.10.0
Your local Nextflow version is too old. Please update to version 22.10.0 or later.
```

**Solution**:
```bash
# Install latest Nextflow
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# Or update existing installation
nextflow self-update

# Verify version
nextflow -version
```

#### 2. Docker Permission Issues

**Problem**: Docker daemon not accessible
```
Got permission denied while trying to connect to the Docker daemon socket
```

**Solution**:
```bash
# Add user to docker group
sudo usermod -aG docker $USER

# Log out and log back in, then verify
docker run hello-world

# Alternative: use sudo for docker commands
sudo nextflow run workflow/main.nf -profile docker
```

#### 3. Docker Image Build Failures

**Problem**: Container build fails during CI or manual build
```
ERROR: failed to solve: process "/bin/sh -c mamba env create -f /tmp/environment.yml" did not complete successfully
```

**Solution**:
```bash
# Check Docker daemon is running
systemctl status docker

# Clean Docker cache and rebuild
docker system prune -af
docker build --no-cache -t variantcalling/pipeline:latest containers/

# Check available disk space
df -h

# Alternative: pull pre-built image
docker pull variantcalling/pipeline:latest
```

### Reference Genome Issues

#### 4. Missing Reference Indices

**Problem**: BWA-MEM2 index files not found
```
[ERROR] Cannot read reference genome index files
```

**Solution**:
```bash
# Check if indices exist
ls reference/hg38.fa.*

# Generate missing indices
cd reference/
bwa-mem2 index hg38.fa
samtools faidx hg38.fa
gatk CreateSequenceDictionary -R hg38.fa -O hg38.dict
```

#### 5. Reference File Paths

**Problem**: Reference file not found
```
[ERROR] Reference FASTA file not found: /path/to/reference.fa
```

**Solution**:
```bash
# Use absolute paths
--ref_fasta /full/path/to/reference.fa

# OR use relative path from run directory
--ref_fasta reference/hg38.fa

# Verify file exists
ls -la $(readlink -f reference/hg38.fa)
```

#### 6. FASTA Sequence Dictionary Issues

**Problem**: GATK requires sequence dictionary
```
A USER ERROR has occurred: FASTA file provided lacks a sequence dictionary
```

**Solution**:
```bash
# Create sequence dictionary
gatk CreateSequenceDictionary \
  -R reference.fa \
  -O reference.dict

# Verify dictionary exists
ls reference.dict
```

### Input Data Issues

#### 7. Samplesheet Format Errors

**Problem**: Samplesheet parsing fails
```
[ERROR] Samplesheet validation failed: Missing required column 'sample_id'
```

**Solution**:
```bash
# Check samplesheet format
head samplesheet.csv

# Required format:
sample_id,read1,read2
NA12878,/path/to/sample_R1.fastq.gz,/path/to/sample_R2.fastq.gz

# Check for common issues:
# - No header line
# - Wrong column names
# - Missing columns
# - Extra whitespace
```

#### 8. FASTQ File Access Issues

**Problem**: Cannot read FASTQ files
```
[ERROR] FASTQ file not found or not readable: /path/to/reads_R1.fastq.gz
```

**Solution**:
```bash
# Check file exists and is readable
ls -la /path/to/reads_R1.fastq.gz

# Check file permissions
chmod 644 /path/to/reads_R1.fastq.gz

# Test file integrity
zcat /path/to/reads_R1.fastq.gz | head -4

# Use absolute paths in samplesheet
/full/path/to/reads_R1.fastq.gz
```

#### 9. Corrupted FASTQ Files

**Problem**: FASTQ parsing errors
```
[ERROR] Malformed FASTQ file: quality string length differs from sequence length
```

**Solution**:
```bash
# Validate FASTQ format
seqtk sample -s 42 reads_R1.fastq.gz 1000 | head -20

# Check for common issues:
# - Truncated files
# - Mixed line endings
# - Non-standard quality encoding

# Re-download or regenerate files if corrupted
```

### Memory and Resource Issues

#### 10. Out of Memory Errors

**Problem**: Java heap space or system memory exhausted
```
java.lang.OutOfMemoryError: Java heap space
# OR
Killed (signal 9)
```

**Solution**:
```bash
# Increase Java memory allocation
export _JAVA_OPTIONS="-Xmx8g"

# Reduce parallel processes
nextflow run workflow/main.nf -profile docker --max_cpus 4 --max_memory 16.GB

# Monitor memory usage
htop

# Check available memory
free -h
```

#### 11. Disk Space Issues

**Problem**: No space left on device
```
ERROR: No space left on device
```

**Solution**:
```bash
# Check disk usage
df -h

# Clean up Nextflow work directory
rm -rf work/

# Clean up Docker containers and images
docker system prune -af

# Move work directory to larger disk
nextflow run workflow/main.nf -w /path/to/large/disk/work

# Clean up intermediate files more aggressively
nextflow run workflow/main.nf -profile docker --cleanup_workdir
```

### Pipeline Execution Issues

#### 12. Process Failures

**Problem**: Individual process fails
```
Process `ALIGN_BWA_MEM2` failed with status 1
```

**Solution**:
```bash
# Check process logs
less .nextflow.log

# Find work directory for failed process
grep "ALIGN_BWA_MEM2" .nextflow.log | grep "work/"

# Examine stderr/stdout in work directory
cat work/12/34567890abcdef/.command.err
cat work/12/34567890abcdef/.command.out

# Check command that was run
cat work/12/34567890abcdef/.command.sh

# Resume from checkpoint
nextflow run workflow/main.nf -profile docker -resume
```

#### 13. Container Issues

**Problem**: Container execution fails
```
[ERROR] Container execution failed: docker: Cannot connect to Docker daemon
```

**Solution**:
```bash
# Check Docker status
systemctl status docker
sudo systemctl start docker

# Test container manually
docker run --rm variantcalling/pipeline:latest bwa-mem2 version

# Use Singularity instead
nextflow run workflow/main.nf -profile singularity

# Use conda instead
nextflow run workflow/main.nf -profile conda
```

#### 14. Network Issues

**Problem**: Cannot download references or containers
```
[ERROR] Failed to download container image
```

**Solution**:
```bash
# Check internet connectivity
curl -I https://hub.docker.com

# Set proxy if needed
export HTTP_PROXY=http://proxy:port
export HTTPS_PROXY=http://proxy:port

# Use local mirror or registry
docker pull local-registry/variantcalling/pipeline:latest

# Pre-download containers
docker pull variantcalling/pipeline:latest
```

### Output and Results Issues

#### 15. Missing Output Files

**Problem**: Expected output files not generated
```
[ERROR] Expected output file not found: results/sample/variants/sample.snv_indel.vcf.gz
```

**Solution**:
```bash
# Check pipeline execution log
grep -i "error\|failed\|exception" .nextflow.log

# Verify all processes completed
grep "Completed process" .nextflow.log

# Check if process was skipped
grep "Skipping process" .nextflow.log

# Run with verbose output
nextflow run workflow/main.nf -profile docker --verbose

# Check publishDir settings in modules.config
```

#### 16. Empty or Invalid VCF Files

**Problem**: VCF file is empty or corrupted
```
[ERROR] No variants found in VCF file
# OR
[ERROR] Invalid VCF header
```

**Solution**:
```bash
# Check VCF file content
zcat results/sample/variants/sample.snv_indel.vcf.gz | head -20

# Validate VCF format
bcftools view results/sample/variants/sample.snv_indel.vcf.gz | head

# Check if variants were filtered out
zcat results/sample/variants/sample.vcf.gz | grep -v "^#" | wc -l

# Examine filtering statistics
grep "VariantFiltration" .nextflow.log

# Adjust filtering parameters
--filter_qual_threshold 10  # Lower quality threshold
```

### Performance Issues

#### 17. Slow Pipeline Execution

**Problem**: Pipeline runs much slower than expected

**Solution**:
```bash
# Check resource allocation
nextflow run workflow/main.nf -with-report report.html
# Open report.html to analyze resource usage

# Increase CPU allocation for specific processes
# Edit workflow/conf/base.config:
process.withName:ALIGN_BWA_MEM2.cpus = 16

# Use faster storage for work directory
nextflow run workflow/main.nf -w /tmp/nextflow-work

# Enable process parallelization
nextflow run workflow/main.nf --max_cpus 32
```

#### 18. High Memory Usage

**Problem**: Processes using excessive memory

**Solution**:
```bash
# Monitor memory usage
htop
# OR
ps aux --sort=-%mem | head

# Reduce memory allocation for Java tools
export _JAVA_OPTIONS="-Xmx4g"

# Split large files
# For large samples, consider splitting FASTQ files

# Use memory-efficient formats
# Pipeline already uses CRAM instead of BAM
```

### Debugging Techniques

#### General Debugging Approach

1. **Check the log file**:
```bash
less .nextflow.log
tail -f .nextflow.log  # Follow log in real-time
```

2. **Enable debug mode**:
```bash
nextflow run workflow/main.nf -profile docker --debug
```

3. **Use resume functionality**:
```bash
nextflow run workflow/main.nf -profile docker -resume
```

4. **Run with stub-run for testing**:
```bash
nextflow run workflow/main.nf -profile docker -stub-run
```

5. **Test individual processes**:
```bash
# Create test data
echo "sample_id,read1,read2" > test.csv
echo "test,test_R1.fq.gz,test_R2.fq.gz" >> test.csv

# Run specific process only
nextflow run workflow/main.nf -profile docker --input test.csv -entry QC_FASTQ
```

### Getting Help

#### 1. Check Documentation
- [Pipeline Overview](overview.md)
- [Inputs and Outputs](inputs_outputs.md)
- [Reference Setup](references.md)

#### 2. Search Known Issues
```bash
# Check GitHub issues
# https://github.com/your-org/variant-calling-pipeline/issues

# Search Nextflow documentation
# https://www.nextflow.io/docs/latest/

# Check tool-specific documentation
# GATK: https://gatk.broadinstitute.org/
# BWA-MEM2: https://github.com/bwa-mem2/bwa-mem2
```

#### 3. Report New Issues

When reporting issues, please include:

1. **Pipeline version**: `git describe --tags`
2. **Nextflow version**: `nextflow -version`
3. **System information**: `uname -a`
4. **Command used**: Full command line
5. **Error message**: Complete error message
6. **Log file**: Relevant sections of `.nextflow.log`
7. **Configuration**: Any custom configuration files used

**Example issue report**:
```
## Description
Pipeline fails during BWA-MEM2 alignment step

## Environment
- Pipeline version: v1.0.0
- Nextflow version: 22.10.6
- OS: Ubuntu 20.04 LTS
- Docker version: 20.10.12

## Command
nextflow run workflow/main.nf -profile docker --input samplesheet.csv --ref_fasta hg38.fa

## Error
Process `ALIGN_BWA_MEM2` failed with status 1

## Log excerpt
[... relevant log lines ...]

## Additional context
Reference genome is GRCh38/hg38 downloaded from UCSC
```

### Advanced Troubleshooting

#### 1. Interactive Debugging

```bash
# Enter failed process work directory
cd work/12/34567890abcdef

# Run command interactively in container
docker run -it --rm -v $(pwd):$(pwd) -w $(pwd) variantcalling/pipeline:latest bash

# Execute the failed command step by step
source .command.sh
```

#### 2. Resource Monitoring

```bash
# Monitor system resources during execution
nohup top -b -d 5 > resource_monitor.log &
nohup iostat -x 5 > io_monitor.log &

# Run pipeline
nextflow run workflow/main.nf -profile docker

# Analyze resource usage
grep "load average" resource_monitor.log
```

#### 3. Network Debugging

```bash
# Test container registry connectivity
nslookup registry-1.docker.io
curl -I https://registry-1.docker.io/v2/

# Test download speeds
wget --spider --server-response http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
```

This troubleshooting guide covers the most common issues encountered when running the variant calling pipeline. For issues not covered here, please check the GitHub repository issues page or open a new issue with detailed information.