process SumAdducts {
    tag "DIMS SumAdducts ${hmdbpart_main_file}"
    label 'SumAdducts'
    container = 'docker://umcugenbioinf/dims:1.4'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       each path(collect_files)
       path(hmdbpart_main_file)

    output:
       path('*_SummedAdducts.RData')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/SumAdducts.R $hmdbpart_main_file $params.preprocessing_scripts_dir $params.zscore
        """
}
