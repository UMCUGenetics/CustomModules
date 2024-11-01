process MakeInit {
    tag "DIMS MakeInit"
    label 'MakeInit'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(samplesheet) 
       val(nr_replicates)
       val(analysis_id)

    output:
       path('run_info.RData')
       path('init.RData')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/MakeInit.R $samplesheet $nr_replicates $analysis_id
        """
}
