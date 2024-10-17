process AssignToBinsStitching {
    tag "DIMS AssignToBinsStitching ${file_id}"
    label 'AssignToBinsStitching'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       tuple(val(file_id), path(mzML_file), path(breaks_file))

    output:
       path("${file_id}.RData"), emit: rdata_file
       path("${file_id}_TIC.txt"), emit: tic_txt_file

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/AssignToBinsStitching.R $mzML_file $breaks_file $params.resolution $params.trim
        """
}


