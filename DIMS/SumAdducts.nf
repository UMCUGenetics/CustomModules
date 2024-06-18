process SumAdducts {
    tag "DIMS SumAdducts"
    label 'SumAdducts'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       each path(collect_files)
       each path(replication_pattern)
       path(HMDBpart_main_file)

    output:
       path('*_SummedAdducts.RData')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/SumAdducts.R $HMDBpart_main_file $params.scripts_dir $params.zscore
        """
}
