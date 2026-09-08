process CollectFilled {
    tag "DIMS CollectFilled"
    label 'CollectFilled'
    container = 'docker://umcugenbioinf/dims:1.3'

    input:
       path(filled_files)
       each path(replication_pattern)
       val(preprocessing_scripts_dir)

    output:
       path('outlist*.txt'), emit: filled_txtfiles
       path('outlist*.RData'), emit: filled_pgrlist

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/CollectFilled.R \
                $preprocessing_scripts_dir \
                $params.ppm \
                $params.zscore
        """
}
