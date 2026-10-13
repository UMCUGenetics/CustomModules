process CollectSumAdducts {
    tag "DIMS CollectSumAdducts"
    label 'CollectSumAdducts'
    container = 'ghcr.io/umcugenetics/dims:v1.4.0'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(collect_files)

    output:
       path('AdductSums_combined.RData'), emit: adductsums_combined
       tuple path('AdductSums_positive.RData'), path('AdductSums_negative.RData'), emit: adductsums_scanmodes

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/CollectSumAdducts.R $params.preprocessing_scripts_dir
        """
}
