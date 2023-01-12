process CheckQC {
    tag {"CheckQC ${identifier}"}
    label 'CheckQC'
    container = 'container_url'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        val(identifier)
        path(input_files)
     
    output:
        path("${identifier}_summary.csv", emit: qc_output)

    script:
        """
        python ${baseDir}/CustomModules/CheckQC/check_qc.py ${params.qc_settings_path} '.' ${identifier} ${input_files} 
        """
}
