process PeakGrouping {
    tag "DIMS PeakGrouping ${HMDBpart_file}"
    label 'PeakGrouping'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(HMDBpart_file)
       each path(SpectrumPeak_file) // input files need to be linked, but called within R script
       each path(pattern_file)      // Execute the process for each element in the input collection (HMDBpart_file)

    output:
       path '*_peaks_used.RData', emit: peaks_used
       path '*_identified.RData', emit: grouped_identified

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/PeakGrouping.R $HMDBpart_file $params.ppm
        """
}
