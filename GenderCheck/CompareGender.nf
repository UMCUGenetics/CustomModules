process CompareGender {
    // Custom process to stats from ExomeDepth analysis
    tag {"CompareGender ${sample_id}"}
    label 'CompareGender'
    container = 'quay.io/biocontainers/pysam:0.22.0--py38h15b938a_0'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        val(sample_id)
        val(gender_clarity)
        path(tuple(bam_file, bai_file))

    output:
        path("${sample_id}_gender_check.txt")

    script:
        """
        python ${baseDir}/CustomModules/GenderCheck/calculate_gender.py ${bam_file} > ${sample_id}_gender_check.txt
        """
