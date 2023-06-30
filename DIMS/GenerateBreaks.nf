process GenerateBreaks {
    label 'GenerateBreaks'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       tuple(file_id, path(mzML_file))


    output:
       file 'breaks.fwhm.RData' 

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/GenerateBreaks.R $mzML_file ./ $params.trim $params.resolution 
        """
}
