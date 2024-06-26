process MakeInit {
    tag "DIMS MakeInit"
    label 'MakeInit'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(samplesheet) 
       val(nr_replicates)

    output:
       path('init.RData')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/MakeInit.R $samplesheet $nr_replicates 
        """
}
