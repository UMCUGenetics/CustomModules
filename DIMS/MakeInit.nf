process MakeInit {
    tag "DIMS MakeInit"
    label 'MakeInit'
    container = 'ghcr.io/umcugenetics/dims:v1.4.0'
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
