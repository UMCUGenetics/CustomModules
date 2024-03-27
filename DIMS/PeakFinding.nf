process PeakFinding {
    tag "DIMS PeakFinding ${rdata_file}"
    label 'PeakFinding'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       tuple(path(rdata_file), path(breaks_file))

    output:
       path '*tive.RData'

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/PeakFinding.R $rdata_file $breaks_file $params.resolution $params.scripts_dir
        """
}
