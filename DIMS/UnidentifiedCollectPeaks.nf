process UnidentifiedCollectPeaks {
    tag "DIMS UnidentifiedCollectPeaks"
    label 'UnidentifiedCollectPeaks'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(SpectrumPeaks_file)
       path(PeakList_identified)

    output:
       path('SpectrumPeaks_*_Unidentified.RData')
       // path('SpectrumPeaks_negative_Unidentified.RData')
       // path('SpectrumPeaks_positive_Unidentified.RData')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/UnidentifiedCollectPeaks.R $SpectrumPeaks_file $params.ppm
        """
}
