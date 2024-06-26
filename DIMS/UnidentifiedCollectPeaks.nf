process UnidentifiedCollectPeaks {
    tag "DIMS UnidentifiedCollectPeaks"
    label 'UnidentifiedCollectPeaks'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(spectrumpeaks_file)
       path(peaklist_identified)

    output:
       path('*.RData')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/UnidentifiedCollectPeaks.R $params.ppm
        """
}
