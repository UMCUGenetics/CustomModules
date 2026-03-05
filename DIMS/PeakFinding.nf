process PeakFinding {
    tag "DIMS PeakFinding ${rdata_file}"
    label 'PeakFinding'
    container = 'docker://umcugenbioinf/dims:1.4'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(rdata_file)
       each path(sample_techreps)

    output:
       path '*tive.RData', emit: peaklist_rdata, optional: true
       path '*tive.txt', optional: true

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/PeakFinding.R \
                --rdata_file $rdata_file \
                --samplesheet $sample_techreps \
                --resolution $params.resolution \
                --preprocessing_scripts_dir $params.preprocessing_scripts_dir
        """
}
