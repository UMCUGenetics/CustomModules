process UnidentifiedFillMissing {
    tag "DIMS UnidentifiedFillMissing ${groupedlist_file}"
    label 'UnidentifiedFillMissing'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(groupedlist_file)
       each path(replication_pattern)

    output:
       path('*_filled.RData')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/UnidentifiedFillMissing.R $groupedlist_file $params.scripts_dir $params.thresh $params.resolution $params.ppm
        """
}
