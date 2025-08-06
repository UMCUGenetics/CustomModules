process FranklinVCF {
    // Custom process to add FILTER status to INFO field as needed by Franklin software
    tag {"FranklinVCF ${sample_id}"}
    label 'FranklinVCF'
    //container = 'ghcr.io/umcugenetics/custommodules_gendercheck:1.0.0'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        path(vcf_file)

    output:
        tuple(path("*_corrected.vcf"), emit: FranklinVCF)

    script:
        output_file = "${vcf_file}_corrected.vcf"
        """
        python ${projectDir}/CustomModules/AddFilterToInfoField/franklinvcf.py \
            $vcf_file \
            $output_file
        """
}
