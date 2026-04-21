process GenerateBreaks {
    tag "DIMS GenerateBreaks"
    label 'GenerateBreaks'
    container = 'docker://umcugenbioinf/dims:1.3'

    input:
       tuple(val(file_id), path(mzML_file))
       val(trim)
       val(resolution)
       val(preprocessing_scripts_dir)

    output:
       path('breaks.fwhm.RData'), emit: breaks
       path('trim_params.RData'), emit: trim_params
       path('highest_mz.RData'), emit: highest_mz

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/GenerateBreaks.R $mzML_file \
                                                               $trim \
                                                               $resolution \
                                                               $preprocessing_scripts_dir
        """
}
