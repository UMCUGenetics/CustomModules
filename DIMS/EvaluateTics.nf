process EvaluateTics {
    tag "DIMS EvaluateTics"
    label 'EvaluateTics'
    container = 'docker://umcugenbioinf/dims:1.3'

    input:
       path(rdata_file)
       path(tic_txt_files)
       path(init_file)
       val(analysis_id)
       path(highest_mz_file)
       path(trim_params_file)
       val(preprocessing_scripts_dir)

    output:
       path('*_repl_pattern.RData'),        emit: pattern_files
       path('replicates_per_sample.txt'),   emit: sample_techreps
       path('miss_infusions_negative.txt'), emit: miss_inf_neg
       path('miss_infusions_positive.txt'), emit: miss_inf_pos
       path('*_TICplots.pdf'),              emit: tic_plots_pdf

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/EvaluateTics.R \
                $init_file \
                $params.nr_replicates \
                $analysis_id \
                $params.matrix \
                $highest_mz_file \
                $trim_params_file \
                $preprocessing_scripts_dir
        """
}

