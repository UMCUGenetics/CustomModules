process HMDBparts {
    tag "DIMS HMDBparts"
    label 'HMDBparts'
    container = 'ghcr.io/umcugenetics/dims:v1.4.0'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        path(hmdb_db_file)
        path(breaks_file)

    output:
        path('*.RData')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/HMDBparts.R $hmdb_db_file $breaks_file $params.standard_run $params.hmdb_parts_files
        """
}
