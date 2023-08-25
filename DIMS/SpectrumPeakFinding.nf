process SpectrumPeakFinding {
    label 'SpectrumPeakFinding'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(RData_file)
       path(replication_pattern)

    output:
       path 'SpectrumPeaks_*.RData'
       // path 'SpectrumPeaks_positive.RData'

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/SpectrumPeakFinding.R $replication_pattern
        """
}
