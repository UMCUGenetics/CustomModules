process AverageTechReplicates {
    tag "DIMS AverageTechReplicates"
    label 'AverageTechReplicates'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(RData_file)
       path(init_filepath)

    output:
       path('*_repl_pattern.RData'), emit: patterns
       path('*_avg.RData'), emit: binned
       path('miss_infusions_negative.txt')
       path('miss_infusions_positive.txt')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/AverageTechReplicates.R $init_filepath $params.nr_replicates
        """
}


