process GenerateQC {
    tag "DIMS GenerateQC"
    label 'GenerateQC'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(outlist_zscores)
       path(adductsums_scanmodes)
       path(identified_files)
       path(init_file) 
       val(analysis_id)

    output:
       path('{*_IS_results, *_positive_controls}.RData'), optional: true
       path('{positive_controls_warning, missing_mz_warning}.txt'), optional: true
       path('{*_IS_SST, *_positive_control}.xlsx'), optional: true
       path('plots/IS_*.png'), emit: plot_files

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/GenerateQC.R $init_file $analysis_id $params.matrix $params.zscore $params.sst_components_file $params.export_scripts_dir
        """
}