process PeakGroupingIdentified {
    label 'PeakGroupingIdentified'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(HMDBpart_file)
       val(resolution)

    output:
       path 'negative.RData'
       path 'positive.RData'

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/PeakGroupingIdentified.R $HMDBpart_file $resolution $ppm
        """
}
