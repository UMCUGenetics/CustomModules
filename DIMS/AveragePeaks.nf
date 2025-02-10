process AveragePeaks {
    tag "DIMS AveragePeaks"
    label 'AveragePeaks'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(rdata_files)
       path(replication_pattern)

    output:
       path 'AvgPeaks_*.RData'

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/AveragePeaks.R 
        """
}
