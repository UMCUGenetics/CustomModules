process HMDBparts {
    tag "DIMS HMDBparts"
    label 'HMDBparts'
    container = 'docker://umcugenbioinf/dims:1.4'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        path(hmdb_db_file)
        path(breaks_file)

    output:
        path('*.RData')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/HMDBparts.R \
                --hmdb_db_file $hmdb_db_file \
                --breaks_file $breaks_file \
                --standard_run $params.standard_run \
                --hmdb_parts_path $params.hmdb_parts_files
        """
}
