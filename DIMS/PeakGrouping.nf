process PeakGrouping {
    tag "DIMS PeakGrouping ${hmdbpart_file}"
    label 'PeakGrouping'
    container = 'docker://umcugenbioinf/dims:1.3'

    input:
       path(hmdbpart_file)
       path(averagedpeaks_file)
       each path(pattern_file)
       val(preprocessing_scripts_dir)

    output:
       path '*_identified.RData', emit: grouped_identified

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/PeakGrouping.R \
                $hmdbpart_file \
                $preprocessing_scripts_dir \
                $ppm
        """
}
