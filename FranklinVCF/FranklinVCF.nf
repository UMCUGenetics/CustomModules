process FranklinVCF {
    // Custom process to add FILTER status to INFO field as needed by Franklin software
    tag {"FranklinVCF ${input_vcf.name}"}
    label 'FranklinVCF'
    container 'ghcr.io/astral-sh/uv:python3.13-alpine'

    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        path(input_vcf)

    output:
        path("${input_vcf.baseName}.franklin.vcf")

    script:
        """
        uv run franklin_vcf.py $input_vcf ${input_vcf.baseName}.franklin.vcf
        """
}
