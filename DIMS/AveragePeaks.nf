process AveragePeaks {
    tag "DIMS AveragePeaks"
    label 'AveragePeaks'
    container = 'ghcr.io/umcugenetics/dims:v1.4.0'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(rdata_files)
       tuple val(sample_id), val(tech_reps), val(scanmode)

    output:
       path 'AvgPeaks_*.RData'

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/AveragePeaks.R $sample_id $tech_reps $scanmode $params.preprocessing_scripts_dir
        """
}
