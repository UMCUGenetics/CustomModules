process AddFilterToInfoField {
    // Custom process to add FILTER status to INFO field
    tag {"AddFilterToInfoField ${sample_id}"}
    label 'AddFilterToInfoField'
    //container = 'ghcr.io/umcugenetics/custommodules_gendercheck:1.0.0'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        path(vcf_file)

    output:
        tuple(path("*_corrected.vcf"), emit: CorrectVCF)

    script:
        output_file = "${vcf_file}_corrected.vcf"
        """
        python ${projectDir}/CustomModules/AddFilterToInfoField/add_filter_to_infofield.py \
            $vcf_file \
            $output_file
        """
}
