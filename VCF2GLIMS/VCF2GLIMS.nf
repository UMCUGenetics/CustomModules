process VCF2GLIMS {
    tag {"VCF2GLIMS ${identifier}"}
    label 'VCF2GLIMS'
    container = 'docker.io/umcugenbioinf/vcf2glims:1.0.0'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        tuple(val(identifier), path(vcf_file))

    output:
        path("${identifier}.csv")

    script:
        """
        python ${baseDir}/CustomModules/VCF2GLIMS/vcf2glims.py ${vcf_file} > ${identifier}.csv
        """
}
