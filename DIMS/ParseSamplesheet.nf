process ParseSamplesheet {
    tag "DIMS ParseSamplesheet"
    label 'ParseSamplesheet'
    container = 'docker://umcugenbioinf/dims:1.3'

    input:
       path(samplesheet) 

    output:
       path('init.RData'), emit: rdata_file
       path('replication_pattern.txt')

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/ParseSamplesheet.R $samplesheet $params.preprocessing_scripts_dir
        """
}
