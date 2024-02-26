process FillMissing {
    tag "DIMS FillMissing ${GroupedList_file}"
    label 'FillMissing'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(GroupedList_file)
       each path(replication_pattern) // input files need to be linked, but called within R script

    output:
       path('*_filled.RData')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/FillMissing.R $GroupedList_file $params.scripts_dir $params.thresh $params.resolution $params.ppm
        """
}
