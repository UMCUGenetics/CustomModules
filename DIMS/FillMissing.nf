process FillMissing {
    label 'FillMissing'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(GroupedList_file)
       path(replication_pattern)

    output:
       path('*_filled.RData')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/FillMissing.R $GroupedList_file $params.scripts_dir $params.thresh $params.resolution $params.ppm
        """
}
