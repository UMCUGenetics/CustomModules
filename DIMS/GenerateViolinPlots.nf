process GenerateViolinPlots {
    tag "DIMS GenerateViolinPlots"
    label 'GenerateViolinPlots'
    container = 'docker://umcugenbioinf/dims:1.3'

    input:
       path(outlist_zscores)
       val(analysis_id)
       val(export_scripts_dir)
       val(path_metabolite_groups)
       val(file_ratios_metabolites)
       val(file_expected_biomarkers_IEM)
       val(file_explanation)

    output:
       path('Diagnostics/*.pdf'), emit: diag_plot_files, optional: true
       path('Other/*.pdf'), emit: other_plot_files, optional: true
       path('dIEM_plots/*.pdf'), emit: diem_plot_files, optional: true
       path('*.xlsx'), emit: excel_file, optional: true
       path('*.csv'), emit: helix_file, optional: true
       path('*.txt'), emit: txt_files, optional: true

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/GenerateViolinPlots.R \
                $analysis_id \
                $export_scripts_dir \
                $path_metabolite_groups \
                $file_ratios_metabolites \
                $file_expected_biomarkers_IEM \
                $file_explanation
        """
}
