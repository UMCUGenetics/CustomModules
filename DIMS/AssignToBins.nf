process AssignToBins {
    tag "DIMS AssignToBins ${file_id}"
    label 'AssignToBins'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       tuple val(file_id), path(mzML_filename) , path(breaks_file)

    output:
       path("${file_id}.RData"), emit: RData_files
       path("${file_id}_TIC.txt"), emit: TIC_txt_files

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/AssignToBins.R $mzML_filename $breaks_file $params.resolution
        """
}


