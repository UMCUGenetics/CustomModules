process SumAdducts {
    tag "DIMS SumAdducts ${hmdbpart_main_file}"
    label 'SumAdducts'
    container = 'ghcr.io/umcugenetics/dims:v1.4.0'

    input:
       each path(collect_files)
       path(hmdbpart_main_file)
       val(zscore)
       val(preprocessing_scripts_dir)

    output:
       path('*_SummedAdducts.RData')

    script:
        """
        Rscript ${projectDir}/CustomModules/DIMS/SumAdducts.R \
                $hmdbpart_main_file \
                $preprocessing_scripts_dir \
                $zscore
        """
}
