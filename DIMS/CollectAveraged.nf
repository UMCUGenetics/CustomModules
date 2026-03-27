process CollectAveraged {
    tag "DIMS CollectAveraged"
    label 'CollectAveraged'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(averaged_files)

    output:
       path('AvgPeaks*.RData'), emit: averaged_peaks

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/CollectAveraged.R
        """
}
