process HMDBparts {
    tag "DIMS HMDBparts"
    label 'HMDBparts'
    container = 'docker://umcugenbioinf/dims:1.3'

    input:
        path(hmdb_db_file)
        path(breaks_file)
        val(hmdb_parts_files)

    output:
        path('*.RData')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/HMDBparts.R \
                $hmdb_db_file \
                $breaks_file \
                $standard_run \
                $hmdb_parts_files
        """
}
