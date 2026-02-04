FROM condaforge/mambaforge:4.14.0-0

# Dockerfile for Variant Calling Pipeline
LABEL maintainer="Pipeline Developer <dev@pipeline.com>"
LABEL description="Container for germline variant calling pipeline"
LABEL version="1.0.0"

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        curl \
        wget \
        unzip \
        git \
        procps \
        ca-certificates && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Create conda environment and install tools with pinned versions
COPY containers/environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml && \
    mamba clean --all -f -y

# Activate conda environment
ENV PATH="/opt/conda/envs/variantcalling/bin:$PATH"
ENV CONDA_DEFAULT_ENV="variantcalling"

# Create working directory
WORKDIR /work

# Copy scripts
COPY scripts/ /usr/local/bin/

# Make scripts executable
RUN chmod +x /usr/local/bin/*.py

# Add labels for tool versions
LABEL bwa-mem2.version="2.2.1"
LABEL samtools.version="1.17"
LABEL gatk4.version="4.4.0.0"
LABEL bcftools.version="1.17"
LABEL fastp.version="0.23.4"
LABEL multiqc.version="1.15"
LABEL python.version="3.9"

# Default command
CMD ["bash"]