process HMDBparts {
    tag "DIMS HMDBparts"
    label 'HMDBparts'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        path(hmdb_db_file)
        path(breaks_file)

    output:
        path('*.RData')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/HMDBparts.R $hmdb_db_file $breaks_file $params.hmdb_parts_files $params.standard_run $params.ppm
        """
}
