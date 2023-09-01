process SumAdducts {
    label 'SumAdducts'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(collect_file) // input files need to be linked, but called within R script
       path(replication_pattern)
       path(HMDBpart_main_file)

    output:
       path('*_SummedAdducts.RData')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/SumAdducts.R $HMDBpart_main_file $params.scripts_dir $params.zscore
        """
}
