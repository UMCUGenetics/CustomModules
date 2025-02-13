process CollectSumAdducts {
    tag "DIMS CollectSumAdducts"
    label 'CollectSumAdducts'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(collect_files)

    output:
       path('AdductSums_combined.RData'), emit: adductsums_combined
       path('AdductSums_{positive, negative}.RData'), emit: adductsums_scanmodes

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/CollectSumAdducts.R $params.preprocessing_scripts_dir
        """
}
