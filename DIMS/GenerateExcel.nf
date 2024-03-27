process GenerateExcel {
    tag "DIMS GenerateExcel"
    label 'GenerateExcel'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(collect_files)
       path(identified_files)
       path(init_file) 
       val(analysis_id) 
       path(relevance_file)

    output:
       path('AdductSums_*.txt')
       path('*IS_results.RData')
       path('*.xlsx'), emit: excel_file 
       path('plots'), emit: plot_files 

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/GenerateExcel.R $init_file $analysis_id $params.matrix $relevance_file $params.zscore 
        """
}
