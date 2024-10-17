process GenerateBreaksStitching {
    tag "DIMS GenerateBreaksStitching"
    label 'GenerateBreaksStitching'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       tuple(val(file_id), path(mzML_file))


    output:
       path('breaks.fwhm.RData'), emit: breaks

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/GenerateBreaksStitching.R $mzML_file ./ $params.trim $params.resolution
        """
}
