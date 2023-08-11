process SpectrumPeakFinding {
    label 'SpectrumPeakFinding'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(RData_file)

    output:
       path 'negative.RData'
       path 'positive.RData'

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/SpectrumPeakFinding.R $RData_file
        """
}
