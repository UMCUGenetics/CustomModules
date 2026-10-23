process AssignToBins {
    tag "DIMS AssignToBins ${file_id}"
    label 'AssignToBins'
    container = 'ghcr.io/umcugenetics/dims:v1.4.0'

    input:
       tuple(val(file_id), path(mzML_file), path(breaks_file), path(trim_params_file))
       val(resolution)

    output:
       path("${file_id}.RData"), emit: rdata_file
       path("${file_id}_TIC.txt"), emit: tic_txt_file

    script:
        """
        Rscript ${projectDir}/CustomModules/DIMS/AssignToBins.R \
                $mzML_file \
                $breaks_file \
                $trim_params_file \
                $resolution
        """
}


