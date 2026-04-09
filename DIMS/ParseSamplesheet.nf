process ParseSamplesheet {
    tag "DIMS ParseSamplesheet"
    label 'ParseSamplesheet'
    container = 'docker://umcugenbioinf/dims:1.3'

    input:
       path(samplesheet) 
       val(preprocessing_scripts_dir)

    output:
       path('init.RData'), emit: rdata_file
       path('replication_pattern.txt'), emit: repl_pattern_txtfile

    script:
        """
        Rscript ${baseDir}/CustomModules/DIMS/ParseSamplesheet.R $samplesheet $preprocessing_scripts_dir
        """
}
