process GenerateExcel {
    tag "DIMS GenerateExcel"
    label 'GenerateExcel'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(adductsums_combined)
       val(analysis_id) 
       path(relevance_file)

    output:
       path('outlist.RData'), emit: outlist_zscores
       path('{AdductSums_filtered_Zscores, AdductSums_filtered_robustZ, AdductSums_filtered_outliersremovedZ}.RData'), optional: true
       path('{AdductSums_filtered_Zscores, AdductSums_filtered_robustZ, AdductSums_filtered_outliersremovedZ}.txt'), optional: true
       path('*.xlsx'), emit: project_excel

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/GenerateExcel.R $analysis_id $relevance_file $params.zscore $params.export_scripts_dir
        """
}
