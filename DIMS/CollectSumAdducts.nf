process CollectSumAdducts {
    tag "DIMS CollectSumAdducts"
    label 'CollectSumAdducts'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(collect_file) // input files need to be linked, but called within R script

    output:
       path('AdductSums_*.RData')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/CollectSumAdducts.R 
        """
}
