process PeakFinding {
    tag "DIMS PeakFinding ${rdata_file}"
    label 'PeakFinding'
    container = 'ghcr.io/umcugenetics/dims:v1.4.0'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(rdata_file)
       each path(sample_techreps)

    output:
       path '*tive.RData', emit: peaklist_rdata, optional: true
       path '*tive.txt', optional: true

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/PeakFinding.R $rdata_file $params.resolution $params.preprocessing_scripts_dir
        """
}
