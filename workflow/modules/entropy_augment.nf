process ENTROPY_AUGMENT {
    tag "$meta.sample_id"
    label 'process_high'

    conda "conda-forge::python=3.10 bioconda::pysam=0.22.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'variantcalling/pipeline:latest' }"

    input:
    tuple val(meta), path(bam), path(bai), path(vcf), path(tbi)
    path(ref_fasta)

    output:
    tuple val(meta), path("*.baseline.vcf.gz"), path("*.baseline.vcf.gz.tbi"), emit: baseline
    tuple val(meta), path("*.augmented_variants.vcf.gz"), path("*.augmented_variants.vcf.gz.tbi"), emit: augmented
    tuple val(meta), path("*.novel_entropy_variants.vcf.gz"), path("*.novel_entropy_variants.vcf.gz.tbi"), emit: novel
    tuple val(meta), path("*.low_confidence_variants.vcf.gz"), path("*.low_confidence_variants.vcf.gz.tbi"), emit: low_conf
    tuple val(meta), path("*.entropy_summary.json"), emit: summary
    path "versions.yml"                                              , emit: versions

    when:
    params.run_entropy_augment

    script:
    def prefix = task.ext.prefix ?: "${meta.sample_id}.entropy"
    def args = [
        "--window-size ${params.entropy_window_size}",
        "--window-step ${params.entropy_window_step}",
        "--min-window-reads ${params.entropy_min_window_reads}",
        "--min-mapq ${params.entropy_min_mapq}",
        "--kmer ${params.entropy_kmer}",
        "--max-contigs ${params.entropy_max_contigs}",
        "--entropy-z ${params.entropy_z}",
        "--error-z ${params.entropy_error_z}",
        "--novel-score-threshold ${params.entropy_novel_score}",
        "--low-conf-threshold ${params.entropy_low_conf}",
        "--weight-allele-balance ${params.entropy_weight_allele_balance}",
        "--weight-baseq ${params.entropy_weight_baseq}",
        "--weight-mapq ${params.entropy_weight_mapq}",
        "--weight-entropy ${params.entropy_weight_entropy}",
        "--weight-assembly ${params.entropy_weight_assembly}"
    ].join(' ')

    """
    python3 ${projectDir}/scripts/entropy_variant_pipeline.py \\
        --bam $bam \\
        --reference $ref_fasta \\
        --baseline-vcf $vcf \\
        --baseline-vcf-index $tbi \\
        --output-prefix ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
        pysam: \$(python3 -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """
}
