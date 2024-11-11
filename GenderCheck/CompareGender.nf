process CompareGender {
    // Custom process to check gender of sample with known status
    tag {"CompareGender ${sample_id}"}
    label 'CompareGender'
    label 'CompareGender_Pysam'
    container = 'ghcr.io/umcugenetics/custommodules_gendercheck:1.0.1'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        tuple(val(sample_id), val(analysis_id), path(bam_file), path(bai_file), val(stated_gender))

    output:
        tuple(path("*gendercheck.txt"), emit: gendercheck_qc)

    script:
        """
        python ${projectDir}/CustomModules/GenderCheck/calculate_gender.py \
            ${sample_id} \
            ${analysis_id} \
            ${bam_file} \
            ./ \
            ${stated_gender} \
            $params.gendercheck_ratio_y \
            $params.gendercheck_mapping_qual \
            $params.gendercheck_locus_y
        """
}
