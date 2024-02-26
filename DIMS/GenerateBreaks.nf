process GenerateBreaks {
    tag "DIMS GenerateBreaks"
    label 'GenerateBreaks'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       tuple val(file_id), path(mzML_file)


    output:
       path('breaks.fwhm.RData'), emit: breaks
       path('highest_mz.RData'), emit: highest_mz

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/GenerateBreaks.R $mzML_file ./ $params.trim $params.resolution 
        """
}
