process IGV {
    // Custom process to run BAF analysis
    tag {"BAF IGV ${output_name}"}
    label 'BAF'
    label 'BAF_IGV'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        tuple(val(output_name), path(vcf_files), path(vcf_idx_files))

    output:
        path("${output_name}_baf.igv", emit: BAF_IGV_files)

    script:
        """
        source ${params.baf_path}/venv/bin/activate
        python ${params.baf_path}/make_BAF_igv.py ${vcf_files} -o ${output_name}_baf.igv -c
        """
}
