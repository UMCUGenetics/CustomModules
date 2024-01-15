process AverageTechReplicates {
    tag "DIMS AverageTechReplicates"
    label 'AverageTechReplicates'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(RData_file)     // input files need to be linked, but called within R script
       path(TIC_txt_files)  // input files need to be linked, but called within R script
       path(init_filepath)
       val(analysis_id)
       val(matrix)
       val(highest_mz)

    output:
       path('*_repl_pattern.RData'), emit: patterns
       path('*_avg.RData'), emit: binned
       path('miss_infusions_negative.txt')
       path('miss_infusions_positive.txt')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/AverageTechReplicates.R $init_filepath $params.nr_replicates $analysis_id $matrix
        """
}


