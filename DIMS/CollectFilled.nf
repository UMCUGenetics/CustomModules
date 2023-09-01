process CollectFilled {
    label 'CollectFilled'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(filled_files)
       path(replication_pattern)

    output:
       path('outlist*.RData')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/CollectFilled.R $params.scripts_dir $params.ppm $params.zscore
        """
}
