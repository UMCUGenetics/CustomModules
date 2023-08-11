process AverageTechReplicates {
    label 'AverageTechReplicates'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(RData_file)
       path(init_filepath)
       val(nr_replicates)

    output:
       path '*_*.RData', emit: binned
       path 'repl.pattern.negative.RData'
       path 'repl.pattern.positive.RData'
       path 'miss_infusions_neg.txt'
       path 'miss_infusions_pos.txt'

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/AverageTechReplicates.R $init_filepath $nr_replicates
        """
}


