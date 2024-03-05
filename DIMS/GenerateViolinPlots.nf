process GenerateViolinPlots {
    tag "DIMS GenerateViolinPlots"
    label 'GenerateViolinPlots'
    container = 'docker://umcugenbioinf/dims:1.3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       path(excel_file) // input files need to be linked, but called within R script
       val(analysis_id)

    output:
       // path('*.pdf') // pdf files are generated, but in different directory
       path('*.xlsx')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/GenerateViolinPlots.R $analysis_id $params.scripts_dir $params.zscore
        """
}
