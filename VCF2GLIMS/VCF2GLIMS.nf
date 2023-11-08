process VCF2GLIMS {
    tag {"VCF2GLIMS ${identifier}"}
    label 'VCF2GLIMS'
    container = 'ghcr.io/umcugenetics/vcf2glims:1.0.0'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        tuple(val(analysis_id), val(sample_id), path(vcf_file))

    output:
        path("${analysis_id}_${sample_id}.csv")

    script:
        """
        python ${baseDir}/CustomModules/VCF2GLIMS/vcf2glims.py ${analysis_id} ${vcf_file} > ${analysis_id}_${sample_id}.csv
        """
}
