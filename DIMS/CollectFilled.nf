process CollectFilled {
    tag "DIMS CollectFilled"
    label 'CollectFilled'
    container = 'ghcr.io/umcugenetics/dims:v1.4.0'

    input:
       path(filled_files)
       each path(replication_pattern)
       val(preprocessing_scripts_dir)
       val(ppm)
       val(zscore)

    output:
       path('peakgroup_list*.txt')
       path('peakgroup_list*.RData'), emit: filled_pgrlist

    script:
        """
        Rscript ${projectDir}/CustomModules/DIMS/CollectFilled.R \
                $preprocessing_scripts_dir \
                $params.ppm \
                $params.zscore
        """
}
