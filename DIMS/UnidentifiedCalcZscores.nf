process UnidentifiedCalcZscores {
    tag "DIMS UnidentifiedCalcZscores"
    label 'UnidentifiedCalcZscores'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(unidentified_filled_files)
       each path(replication_pattern)

    output:
       path('unidentified_outlist*.RData')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/UnidentifiedCalcZscores.R $params.scripts_dir $params.ppm $params.zscore
        """
}
