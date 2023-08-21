process PeakGroupingIdentified {
    label 'PeakGroupingIdentified'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(SpectrumPeak_file)
       path(HMDBpart_file)
       path(pattern_file)

    output:
       path '*_negative.RData'
       path '*_positive.RData'

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/PeakGroupingIdentified.R $SpectrumPeak_file $HMDBpart_file $pattern_file $params.resolution $params.ppm
        """
}
