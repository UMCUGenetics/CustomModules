process PeakGrouping {
    label 'PeakGrouping'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(HMDBpart_file)
       path(SpectrumPeak_file)
       path(pattern_file)

    output:
       // path '*_all.RData', emit: peaklist_all
       path '*_identified.RData', emit: grouped_identified

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/PeakGrouping.R $HMDBpart_file $SpectrumPeak_file $pattern_file $params.ppm
        """
}
