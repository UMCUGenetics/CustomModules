process PeakGrouping {
    tag "DIMS PeakGrouping ${hmdbpart_file}"
    label 'PeakGrouping'
    container = 'ghcr.io/umcugenetics/dims:v1.4.0'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(hmdbpart_file)
       path(averagedpeaks_file)
       each path(pattern_file)

    output:
       path '*_identified.RData', emit: grouped_identified

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/PeakGrouping.R $hmdbpart_file $params.preprocessing_scripts_dir $params.ppm
        """
}
