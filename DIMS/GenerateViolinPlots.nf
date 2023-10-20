process GenerateViolinPlots {
    tag {"DIMS GenerateViolinPlots"}
    label 'GenerateViolinPlots'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(excel_file) // input files need to be linked, but called within R script

    output:
       path('*.pdf')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/GenerateViolinPlots.R $params.scripts_dir $params.analysis_id $params.zscore 
        """
}
