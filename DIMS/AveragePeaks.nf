process AveragePeaks {
    tag "DIMS AveragePeaks"
    label 'AveragePeaks'
    container = 'docker://umcugenbioinf/dims:1.3'
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
