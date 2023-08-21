process HMDBparts {
    // Custom process to cut HMDB db into parts
    tag {"DIMS HMDBparts"}
    label 'HMDBparts'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        // tuple(path(hmdb_db_file), path(breaks_file))
        path(hmdb_db_file)
        path(breaks_file)
        // val(standard_run)
        // val(ppm)

    output:
        path('*.RData')

    script:

        """
        Rscript ${baseDir}/CustomModules/DIMS/HMDBparts.R $hmdb_db_file $breaks_file $params.standard_run $params.ppm
        """
}
