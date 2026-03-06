process FillMissing {
    tag "DIMS FillMissing ${peakgrouplist_file}"
    label 'FillMissing'
    container = 'docker://umcugenbioinf/dims:1.4'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(peakgrouplist_file)
       each path(replication_pattern)

    output:
       path('*_filled.RData')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/FillMissing.R \
                --peakgrouplist_file $peakgrouplist_file \
                --preprocessing_scripts_dir $params.preprocessing_scripts_dir \
                --thresh_noise $params.thresh
        """
}
