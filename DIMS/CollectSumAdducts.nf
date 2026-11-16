process CollectSumAdducts {
    tag "DIMS CollectSumAdducts"
    label 'CollectSumAdducts'
    container = 'ghcr.io/umcugenetics/dims:v1.4.0'

    input:
       path(collect_files)
       val(preprocessing_scripts_dir)

    output:
       path('AdductSums_combined.RData'), emit: adductsums_combined
       tuple path('AdductSums_positive.RData'), path('AdductSums_negative.RData'), emit: adductsums_scanmodes

    script:
        """
        Rscript ${projectDir}/CustomModules/DIMS/CollectSumAdducts.R \
                $preprocessing_scripts_dir
        """
}
