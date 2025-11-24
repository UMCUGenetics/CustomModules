process GenerateExcel {
    tag "DIMS GenerateExcel"
    label 'GenerateExcel'
    container = 'docker://umcugenbioinf/dims:1.4'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(adductsums_combined)
       val(analysis_id) 
       path(relevance_file)

    output:
       path("outlist.RData"), emit: outlist_zscores
       tuple path("AdductSums_filtered_Zscores.RData"), path("AdductSums_filtered_robustZ.RData"), path("AdductSums_filtered_outliersremovedZ.RData"), optional: true
       tuple path("AdductSums_filtered_Zscores.txt"), path("AdductSums_filtered_robustZ.txt"), path("AdductSums_filtered_outliersremovedZ.txt"), optional: true
       path("${analysis_id}.xlsx"), emit: project_excel
       path("Helix_${analysis_id}.xlsx"), optional: true

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/GenerateExcel.R \\
            --project $analysis_id \\
            --hmdb_rlvc_file $relevance_file \\
            --z_score $params.zscore \\
            --export_scripts_dir $params.export_scripts_dir \\
            --path_metabolite_groups $params.path_metabolite_groups
        """
}
