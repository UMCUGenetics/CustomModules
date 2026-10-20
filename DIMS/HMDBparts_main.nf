process HMDBparts_main {
    tag "DIMS HMDBparts_main"
    label 'HMDBparts_main'
    container = 'ghcr.io/umcugenetics/dims:v1.4.0'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        path(hmdb_db_file)
        path(breaks_file)

    output:
        path('*.RData')

    script:

        """
        Rscript ${baseDir}/CustomModules/DIMS/HMDBparts_main.R $hmdb_db_file $breaks_file 
        """
}
