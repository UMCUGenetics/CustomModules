process FranklinVCF {
    // Custom process to add FILTER status to INFO field as needed by Franklin software
    tag {"FranklinVCF ${analysis_id}"}
    label 'FranklinVCF'
    container 'ghcr.io/astral-sh/uv:python3.13-bookworm'

    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        tuple val(analysis_id), path(input_vcf), path(input_vcf_idx)

    output:
        tuple val(analysis_id), path("${input_vcf.baseName}.franklin.vcf")

    script:
        """
        uv run --no-cache franklin_vcf.py $input_vcf ${input_vcf.baseName}.franklin.vcf
        """
}
