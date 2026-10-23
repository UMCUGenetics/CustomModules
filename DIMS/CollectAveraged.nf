process CollectAveraged {
    tag "DIMS CollectAveraged"
    label 'CollectAveraged'
    container = 'ghcr.io/umcugenetics/dims:v1.4.0'

    input:
       path(averaged_files)

    output:
       path('AvgPeaks*.RData'), emit: averaged_peaks

    script:
        """
        Rscript ${projectDir}/CustomModules/DIMS/CollectAveraged.R
        """
}
