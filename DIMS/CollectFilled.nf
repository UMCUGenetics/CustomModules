process CollectFilled {
    tag "DIMS CollectFilled"
    label 'CollectFilled'
    container = 'ghcr.io/umcugenetics/dims:v1.4.0'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(filled_files)
       each path(replication_pattern)

    output:
       path('peakgroup_list*.txt')
       path('peakgroup_list*.RData'), emit: filled_pgrlist

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/CollectFilled.R $params.preprocessing_scripts_dir $params.zscore
        """
}
