process AverageTechReplicates {
    tag "DIMS AverageTechReplicates"
    label 'AverageTechReplicates'
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

    output:
       path('*_repl_pattern.RData'), emit: pattern_files
       path('*_avg.RData'), emit: binned_files
       path('miss_infusions_negative.txt')
       path('miss_infusions_positive.txt')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/AverageTechReplicates.R $init_file \
                                                                      $params.nr_replicates \
                                                                      $analysis_id $matrix \
                                                                      $highest_mz_file
        """
}


