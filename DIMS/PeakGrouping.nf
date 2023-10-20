process PeakGrouping {
    tag "DIMS PeakGrouping ${HMDBpart_file}"
    label 'PeakGrouping'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(HMDBpart_file)
       path(SpectrumPeak_file) // input files need to be linked, but called within R script
       path(pattern_file)

    output:
       path '*_peaks_used.RData', emit: peaks_used
       path '*_identified.RData', emit: grouped_identified

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/PeakGrouping.R $HMDBpart_file $params.ppm
        """
}
