process UnidentifiedCalcZscores {
    tag "DIMS UnidentifiedCalcZscores"
    label 'UnidentifiedCalcZscores'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(unidentified_filled_files)
       path(replication_pattern)  // input files need to be linked, but called within R script

    output:
       path('unidentified_outlist*.txt')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/UnidentifiedCalcZscores.R $params.scripts_dir $params.ppm $params.zscore
        """
}
