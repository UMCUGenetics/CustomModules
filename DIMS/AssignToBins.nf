process AssignToBins {
    label 'AssignToBins'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       tuple(file_id, path(mzML_filename))
       file(breaks_file)
       val(resolution)

    output:
       file 'AssignToBins.RData'

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/AssignToBins.R $mzML_filename $breaks_file $resolution
        """
}


