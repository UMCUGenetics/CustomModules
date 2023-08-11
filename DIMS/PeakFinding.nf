process PeakFinding {
    label 'PeakFinding'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(RData_file)
       val(resolution)
       path(scripts)
       path(breaks_file)

    output:
       path '*_neg.RData'
       path '*_pos.RData'
       path 'repl.pattern.*.RData', emit: pattern
       path 'miss_infusions_neg.txt'
       path 'miss_infusions_pos.txt'

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/PeakFinding.R $RData_file $resolution $scripts $breaks_file
        """
}
