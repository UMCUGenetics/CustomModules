process SingleIGV {
    // Custom process to run Single sample IGV analysis
    tag {"ExomeDepth SingleIGV ${sample_id}"}
    label 'ExomeDepth'
    label 'ExomeDepth_SingleIGV'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        tuple(val(sample_id), val(analysis_id), val(refset))

    output:
        path("*.xml", emit: Single_IGV_file)

    script:
        """
        source ${params.dx_resources_path}/${params.exomedepth_path}/venv/bin/activate
        python ${params.dx_resources_path}/${params.exomedepth_path}/igv_xml_session.py single_igv ./ ${sample_id} ${analysis_id} ${refset}
        """
}

process FamilyIGV {
    // Custom process to run Family IGV analysis
    tag {"ExomeDepth FamilyIGV ${sample_id}"}
    label 'ExomeDepth'
    label 'ExomeDepth_FamilyIGV'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        val(ped_file)
        val(analysis_id)
        path(bam_files)
        path(snv_vcf_files)
        path(cnv_vcf_files)
        path(igv_files)
        path(upd_files)
        path(baf_files)

    output:
        path("*.xml", emit: Family_IGV_file)

    script:
        def bam_files = bam_files.collect{"$it"}.join(" ")
        def snv_vcf_files = snv_vcf_files.collect{"$it"}.join(" ")
        def cnv_vcf_files = cnv_vcf_files.collect{"$it"}.join(" ")
        def igv_files = igv_files.collect{"$it"}.join(" ")
        def upd_files = upd_files.collect{"$it"}.join(" ")
        def baf_files = baf_files.collect{"$it"}.join(" ")
        """
        source ${params.dx_resources_path}/${params.exomedepth_path}/venv/bin/activate
        python ${params.dx_resources_path}/${params.exomedepth_path}/igv_xml_session.py family_igv ./ ${ped_file} ${analysis_id} \
            --bam_files $bam_files \
            --snv_vcf_files $snv_vcf_files \
            --cnv_vcf_files $cnv_vcf_files \
            --igv_files $igv_files \
            --upd_files $upd_files  \
            --baf_files $baf_files
        """
}
