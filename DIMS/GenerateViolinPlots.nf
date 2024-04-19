process GenerateViolinPlots {
    tag "DIMS GenerateViolinPlots"
    label 'GenerateViolinPlots'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(excel_file)
       val(analysis_id)

    output:
       path('Diagnostics/*.pdf'), emit: diag_plot_files
       path('Other/*.pdf'), emit: other_plot_files
       path('dIEM_plots/*.pdf'), emit: diem_plot_files
       path('*.xlsx'), emit: excel_file
       path('*.csv'), emit: helix_file

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/GenerateViolinPlots.R $analysis_id $params.scripts_dir $params.zscore \
                                                                    $params.path_metabolite_groups \
                                                                    $params.file_ratios_metabolites \
                                                                    $params.file_expected_biomarkers_IEM \
                                                                    $params.file_explanation \
                                                                    $params.file_isomers
        """
}
