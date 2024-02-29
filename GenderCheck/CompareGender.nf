process CompareGender {
    // Custom process to check gender of sample with known status
    tag {"CompareGender ${sample_id}"}
    label 'CompareGender'
    label 'CompareGender_Pysam'
    container = 'quay.io/biocontainers/pysam:0.22.0--py38h15b938a_0'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        tuple(sample_id, analysis_id, path(bam_file), path(bai_file), val(true_gender))

    output:
        tuple(path("*gendercheck.txt"), emit: gendercheck_qc)

    script:
        """
        python ${baseDir}/CustomModules/GenderCheck/calculate_gender.py \
            ${sample_id} \
            ${analysis_id} \
            ${bam_file} \
            ./ \
            ${true_gender} \
            $params.gendercheck_ratio_y \
            $params.gendercheck_mapping_qual \
            $params.gendercheck_locus_y
        """
}
