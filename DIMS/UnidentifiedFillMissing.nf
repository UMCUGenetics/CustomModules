process UnidentifiedFillMissing {
    tag "DIMS UnidentifiedFillMissing ${GroupedList_file}"
    label 'UnidentifiedFillMissing'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(GroupedList_file)
       path(replication_pattern) // input files need to be linked, but called within R script

    output:
       path('*_filled.RData')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/UnidentifiedFillMissing.R $GroupedList_file $params.scripts_dir $params.thresh $params.resolution $params.ppm
        """
}
