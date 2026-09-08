process FillMissing {
    tag "DIMS FillMissing ${peakgrouplist_file}"
    label 'FillMissing'
    container = 'docker://umcugenbioinf/dims:1.3'

    input:
       path(peakgrouplist_file)
       each path(replication_pattern)
       val(preprocessing_scripts_dir)

    output:
       path('*_filled.RData'), emit: filled_data

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/FillMissing.R \
                $peakgrouplist_file \
                $preprocessing_scripts_dir \
                $params.thresh
        """
}
