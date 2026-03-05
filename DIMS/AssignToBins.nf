process AssignToBins {
    tag "DIMS AssignToBins ${file_id}"
    label 'AssignToBins'
    container = 'docker://umcugenbioinf/dims:1.4'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       tuple(val(file_id), path(mzML_file), path(breaks_file), path(trim_params_file))

    output:
       path("${file_id}.RData"), emit: rdata_file
       path("${file_id}_TIC.txt"), emit: tic_txt_file

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/AssignToBins.R \
                --mzML_filepath $mzML_file \
                --breaks_filepath $breaks_file \
                --trim_params_filepath $trim_params_file
        """
}


