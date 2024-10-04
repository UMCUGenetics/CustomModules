process CheckQC {
    tag {"CheckQC ${identifier}"}
    label 'CheckQC'
    container = 'docker.io/umcugenbioinf/checkqc:1.0.0'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        val(identifier)
        path(input_files)

    output:
        path("${identifier}_summary.csv", emit: qc_output)

    script:
        """
        python ${projectDir}/CustomModules/CheckQC/check_qc.py ${params.qc_settings_path} '.' ${identifier} ${input_files}
        """
}
