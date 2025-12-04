process GenerateQCOutput {
    tag "DIMS GenerateQCOutput"
    label 'GenerateQCOutput'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(outlist_zscores)
       path(adductsums_scanmodes)
       path(identified_files)
       path(init_file) 
       val(analysis_id)

    output:
       tuple path('*_IS_results.RData'), path('*_positive_control.RData'), optional: true
       tuple path('positive_controls_warning.txt'), path('missing_mz_warning.txt'), path('sample_names_nodata.txt'), optional: true
       tuple path('*_IS_SST.xlsx'), path('*_positive_control.xlsx'), optional: true
       path('plots/IS_*.png'), emit: plot_files
       path('Check_number_of_controls.txt'), optional: true

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/GenerateQCOutput.R $init_file \
                                                                 $analysis_id \
                                                                 $params.matrix \
                                                                 $params.sst_components_file \
                                                                 $params.export_scripts_dir
        """
}
