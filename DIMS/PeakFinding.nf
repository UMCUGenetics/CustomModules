process PeakFinding {
    label 'PeakFinding'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       // path(RData_file)
       // path(breaks_file)
       tuple(path(RData_file), path(breaks_file))

    output:
       path '*tive.RData'

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/PeakFinding.R $RData_file $breaks_file $params.resolution $params.scripts_dir
        """
}
