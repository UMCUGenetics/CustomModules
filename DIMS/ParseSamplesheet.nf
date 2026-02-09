process MakeInit {
    tag "DIMS MakeInit"
    label 'MakeInit'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(samplesheet) 

    output:
       path('init.RData')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/MakeInit.R $samplesheet $params.preprocessing_scripts_dir
        """
}
