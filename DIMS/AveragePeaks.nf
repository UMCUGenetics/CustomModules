process AveragePeaks {
    tag "DIMS AveragePeaks"
    label 'AveragePeaks'
    container = 'docker://umcugenbioinf/dims:1.3'

    input:
       path(rdata_files)
       tuple val(sample_id), val(tech_reps), val(scanmode)
       val(preprocessing_scripts_dir)

    output:
       path 'AvgPeaks_*.RData', emit: average_peaks

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/AveragePeaks.R \
                $sample_id \
                $tech_reps \
                $scanmode \
                $preprocessing_scripts_dir
        """
}
