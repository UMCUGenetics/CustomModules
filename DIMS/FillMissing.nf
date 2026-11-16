process FillMissing {
    tag "DIMS FillMissing ${peakgrouplist_file}"
    label 'FillMissing'
    container = 'ghcr.io/umcugenetics/dims:v1.4.0'

    input:
       path(peakgrouplist_file)
       each path(replication_pattern)
       val(preprocessing_scripts_dir)
       val(threshold)

    output:
       path('*_filled.RData'), emit: filled_data

    script:
        """
        Rscript ${projectDir}/CustomModules/DIMS/FillMissing.R \
                $peakgrouplist_file \
                $preprocessing_scripts_dir \
                $threshold
        """
}
