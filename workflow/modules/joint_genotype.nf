process JOINT_GENOTYPE {
    tag "joint_genotyping"
    label 'process_high'

    conda "bioconda::gatk4=4.4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0' :
        'variantcalling/pipeline:latest' }"

    input:
    path(gvcf_files)
    path(ref_fasta)

    output:
    path("joint_genotyped.vcf.gz"), path("joint_genotyped.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def gvcf_list = gvcf_files.collect { "-V $it" }.join(' ')
    
    """
    # Create genomics database
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        GenomicsDBImport \\
        $gvcf_list \\
        --genomicsdb-workspace-path genomicsdb \\
        --intervals chr1:1-1000000 \\
        --tmp-dir ./tmp

    # Joint genotyping
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        GenotypeGVCFs \\
        -R $ref_fasta \\
        -V gendb://genomicsdb \\
        -O joint_genotyped.vcf.gz \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}