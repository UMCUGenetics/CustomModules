process CheckQC {
    tag {"CheckQC ${identifier}"}
    label 'CheckQC'
    container = 'docker.io/umcugenbioinf/checkqc:1.0.0'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        val(identifier)
        path(input_files, stageAs: "?/*")
     
    output:
        path("${identifier}_checkqc_summary.csv", emit: qc_output)

    script:
        """
        python ${baseDir}/CustomModules/CheckQC/check_qc.py ${params.qc_settings_path} '.' ${identifier}_checkqc ${input_files} 
        """
}
