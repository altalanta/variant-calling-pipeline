process MARKDUP {
    tag "$meta.sample_id"
    label 'process_medium'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'variantcalling/pipeline:latest' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.cram"), path("*.cram.crai"), emit: bam
    tuple val(meta), path("*.metrics.txt")               , emit: metrics
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"
    
    """
    samtools collate \\
        -o ${prefix}.collated.bam \\
        $bam

    samtools fixmate \\
        -m \\
        ${prefix}.collated.bam \\
        ${prefix}.fixmate.bam

    samtools sort \\
        -o ${prefix}.sorted.bam \\
        ${prefix}.fixmate.bam

    samtools markdup \\
        $args \\
        -f ${prefix}.markdup.metrics.txt \\
        --write-index \\
        --output-fmt CRAM \\
        ${prefix}.sorted.bam \\
        ${prefix}.cram

    # Clean up intermediate files
    rm ${prefix}.collated.bam ${prefix}.fixmate.bam ${prefix}.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}