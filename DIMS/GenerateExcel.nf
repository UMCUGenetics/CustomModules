process GenerateExcel {
    tag {"DIMS GenerateExcel"}
    label 'GenerateExcel'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(collect_file) // input files need to be linked, but called within R script
       path(init_filepath) 

    output:
       path('AdductSums_*.RData')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/GenerateExcel.R $init_filepath $params.analysis_id $params.matrix $params.relevance_file $params.zscore 
        """
}
