process SpectrumPeakFinding {
    tag "DIMS SpectrumPeakFinding"
    label 'SpectrumPeakFinding'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(rdata_files)
       path(replication_pattern)

    output:
       path 'SpectrumPeaks_*.RData'

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/SpectrumPeakFinding.R 
        """
}
