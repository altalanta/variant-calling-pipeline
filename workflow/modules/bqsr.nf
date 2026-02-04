process BQSR {
    tag "$meta.sample_id"
    label 'process_medium'

    conda "bioconda::gatk4=4.4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0' :
        'variantcalling/pipeline:latest' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(ref_fasta)
    path(known_sites)

    output:
    tuple val(meta), path("*.cram"), path("*.cram.crai"), emit: bam
    tuple val(meta), path("*.recal_data.table")          , emit: table
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"
    def known_sites_cmd = known_sites ? "--known-sites $known_sites" : ""
    
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        BaseRecalibrator \\
        -R $ref_fasta \\
        -I $bam \\
        $known_sites_cmd \\
        -O ${prefix}.recal_data.table \\
        $args

    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        ApplyBQSR \\
        -R $ref_fasta \\
        -I $bam \\
        --bqsr-recal-file ${prefix}.recal_data.table \\
        -O ${prefix}.cram \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}