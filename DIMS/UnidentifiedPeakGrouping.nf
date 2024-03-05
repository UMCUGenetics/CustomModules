process UnidentifiedPeakGrouping {
    tag "DIMS UnidentifiedPeakGrouping ${UnidentifiedSpectrumPeaks_file}"
    label 'UnidentifiedPeakGrouping'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(UnidentifiedSpectrumPeaks_file) // input files need to be linked, but called within R script
       path(replication_pattern)  // input files need to be linked, but called within R script

    output:
       path('*_Unidentified.txt')
       path('*_Unidentified.RData'), emit: grouped_unidentified

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/UnidentifiedPeakGrouping.R $params.resolution $params.ppm
        """
}
