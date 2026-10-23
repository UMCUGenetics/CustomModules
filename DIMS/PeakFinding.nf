process PeakFinding {
    tag "DIMS PeakFinding ${rdata_file}"
    label 'PeakFinding'
    container = 'ghcr.io/umcugenetics/dims:v1.4.0'

    input:
       path(rdata_file)
       each path(sample_techreps)
       val(resolution)
       val(preprocessing_scripts_dir)

    output:
       path '*tive.RData', emit: peaklist_rdata, optional: true
       path '*tive.txt', emit: peaklist_txt, optional: true

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/PeakFinding.R \
                $rdata_file \
                $resolution \
                $preprocessing_scripts_dir
        """
}
