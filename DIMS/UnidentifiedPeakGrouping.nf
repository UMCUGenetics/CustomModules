process UnidentifiedPeakGrouping {
    tag "DIMS UnidentifiedPeakGrouping ${unidentified_spectrumpeaks_files}"
    label 'UnidentifiedPeakGrouping'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(unidentified_spectrumpeaks_files)
       each path(replication_pattern)

    output:
       path('*_Unidentified.txt')
       path('*_Unidentified.RData'), emit: grouped_unidentified

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/UnidentifiedPeakGrouping.R $unidentified_spectrumpeaks_files $params.resolution $params.ppm
        """
}
