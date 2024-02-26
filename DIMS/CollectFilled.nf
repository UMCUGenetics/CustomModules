process CollectFilled {
    tag "DIMS CollectFilled"
    label 'CollectFilled'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(filled_files)
       each path(replication_pattern)  // input files need to be linked, but called within R script

    output:
       path('outlist*.txt')
       path('outlist*.RData'), emit: filled_pgrlist

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/CollectFilled.R $params.scripts_dir $params.ppm $params.zscore
        """
}
