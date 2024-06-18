process FillMissing {
    tag "DIMS FillMissing ${peakgrouplist_file}"
    label 'FillMissing'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(peakgrouplist_file)
       each path(replication_pattern)

    output:
       path('*_filled.RData')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/FillMissing.R $peakgrouplist_file $params.scripts_dir $params.thresh $params.resolution $params.ppm
        """
}
