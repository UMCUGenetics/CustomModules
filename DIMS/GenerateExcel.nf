process GenerateExcel {
    tag "DIMS GenerateExcel"
    label 'GenerateExcel'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(collect_file) // input files need to be linked, but called within R script
       path(init_filepath) 
       val(analysis_id) 
       path(relevance_file)

    output:
       path('AdductSums_*.RData')
       path('*.xlsx'), emit: excel_file 

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/GenerateExcel.R $init_filepath $analysis_id $params.matrix $relevance_file $params.zscore 
        """
}
