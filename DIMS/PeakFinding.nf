process PeakFinding {
    tag "DIMS PeakFinding ${rdata_file}"
    label 'PeakFinding'
    container = 'docker://umcugenbioinf/dims:1.3'

    input:
       path(rdata_file)
       each path(sample_techreps)
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
