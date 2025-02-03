process EvaluateTics {
    tag "DIMS EvaluateTics"
    label 'EvaluateTics'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(rdata_file)
       path(tic_txt_files)
       path(init_file)
       val(nr_replicates)
       val(analysis_id)
       val(matrix)
       path(highest_mz_file)
       path(trim_params_file)

    output:
       path('*_repl_pattern.RData'), emit: pattern_files
       path('miss_infusions_negative.txt')
       path('miss_infusions_positive.txt')
       path('*_TICplots.pdf')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/EvaluateTics.R $init_file \
                                                             $params.nr_replicates \
                                                             $analysis_id \
                                                             $matrix \
                                                             $highest_mz_file \
                                                             $trim_params_file
        """
}


