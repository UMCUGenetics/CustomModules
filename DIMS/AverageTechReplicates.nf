process AverageTechReplicates {
    tag "DIMS AverageTechReplicates"
    label 'AverageTechReplicates'
    container = 'docker://umcugenbioinf/dims:1.4'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(rdata_file)
       path(tic_txt_files)
       path(init_file)
       val(nr_replicates)
       val(analysis_id)
       val(matrix)
       path(highest_mz_file)
       path(breaks_file)

    output:
       path('*_repl_pattern.RData'), emit: pattern_files
       path('*_avg.RData'), emit: binned_files
       path('miss_infusions_negative.txt')
       path('miss_infusions_positive.txt')
       path('*_TICplots.pdf'), emit: tic_plots_pdf

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/AverageTechReplicates.R \\
            --init_filepath $init_file \\
            --nr_replicates $params.nr_replicates \\
            --run_name $analysis_id \\
            --matrix $matrix \\
            --highest_mz_file $highest_mz_file \\
            --breaks_filepath $breaks_file
        """
}


