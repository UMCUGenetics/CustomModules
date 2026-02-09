process ParseSamplesheet {
    tag "DIMS ParseSamplesheet"
    label 'ParseSamplesheet'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(samplesheet) 

    output:
       path('init.RData')
       path('technical_replicates.txt')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/ParseSamplesheet.R $samplesheet $params.preprocessing_scripts_dir
        """
}
