process PeakGrouping {
    tag "DIMS PeakGrouping ${hmdbpart_file}"
    label 'PeakGrouping'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(hmdbpart_file)
       each path(spectrumpeak_file)
       each path(pattern_file)

    output:
       path '*_peaks_used.RData', emit: peaks_used
       path '*_identified.RData', emit: grouped_identified

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/PeakGrouping.R $hmdbpart_file $params.ppm
        """
}
