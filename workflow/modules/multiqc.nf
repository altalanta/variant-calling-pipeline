process MULTIQC {
    label 'process_single'

    conda "bioconda::multiqc=1.15"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.15--pyhdfd78af_0' :
        'variantcalling/pipeline:latest' }"

    input:
    path(multiqc_files)
    path(multiqc_config)
    path(multiqc_custom_config)
    path(software_versions)

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def config = multiqc_config ? "--config $multiqc_config" : ''
    def custom_config = multiqc_custom_config ? "--config $multiqc_custom_config" : ''
    
    """
    multiqc \\
        -f \\
        $args \\
        $config \\
        $custom_config \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}