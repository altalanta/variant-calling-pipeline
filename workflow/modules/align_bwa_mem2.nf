process ALIGN_BWA_MEM2 {
    tag "$meta.sample_id"
    label 'process_high'

    conda "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' :
        'variantcalling/pipeline:latest' }"

    input:
    tuple val(meta), path(reads)
    path(ref_fasta)

    output:
    tuple val(meta), path("*.sam"), emit: sam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"
    def read_group = "@RG\\tID:${meta.sample_id}\\tSM:${meta.sample_id}\\tPL:${meta.platform ?: 'ILLUMINA'}\\tLB:${meta.library ?: meta.sample_id}\\tPU:${meta.lane ?: '1'}"
    
    """
    INDEX=`find -L . -name "*.amb" | sed 's/.amb\$//'`

    bwa-mem2 mem \\
        $args \\
        -t $task.cpus \\
        -R \"$read_group\" \\
        \$INDEX \\
        ${reads[0]} \\
        ${reads[1]} \\
        > ${prefix}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
    END_VERSIONS
    """
}