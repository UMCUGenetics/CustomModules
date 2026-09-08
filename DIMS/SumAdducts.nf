process SumAdducts {
    tag "DIMS SumAdducts ${hmdbpart_main_file}"
    label 'SumAdducts'
    container = 'docker://umcugenbioinf/dims:1.3'

    input:
       each path(collect_files)
       path(hmdbpart_main_file)
       val(preprocessing_scripts_dir)

    output:
       path('*_SummedAdducts.RData')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/SumAdducts.R \
                $hmdbpart_main_file \
                $preprocessing_scripts_dir \
                $zscore
        """
}
