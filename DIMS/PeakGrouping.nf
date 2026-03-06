process PeakGrouping {
    tag "DIMS PeakGrouping ${hmdbpart_file}"
    label 'PeakGrouping'
    container = 'docker://umcugenbioinf/dims:1.4'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(hmdbpart_file)
       path(averagedpeaks_file)
       each path(pattern_file)

    output:
       path '*_identified.RData', emit: grouped_identified

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/PeakGrouping.R \
                --hmdbpart_file $hmdbpart_file \
                --preprocessing_scripts_dir $params.preprocessing_scripts_dir \
                --ppm $params.ppm
        """
}
