process CollectSumAdducts {
    tag "DIMS CollectSumAdducts"
    label 'CollectSumAdducts'
    container = 'docker://umcugenbioinf/dims:1.3'

    input:
       path(collect_files)
       val(preprocessing_scripts_dir)

    output:
       path('AdductSums_combined.RData'), emit: adductsums_combined
       tuple path('AdductSums_positive.RData'), path('AdductSums_negative.RData'), emit: adductsums_scanmodes

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/CollectSumAdducts.R \
                $preprocessing_scripts_dir
        """
}
