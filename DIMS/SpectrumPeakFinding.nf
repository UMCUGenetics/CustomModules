process SpectrumPeakFinding {
    tag "DIMS SpectrumPeakFinding"
    label 'SpectrumPeakFinding'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(RData_file) // input files need to be linked, but called within R script
       path(replication_pattern) // input files need to be linked, but called within R script

    output:
       path 'SpectrumPeaks_*.RData'

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/SpectrumPeakFinding.R 
        """
}
