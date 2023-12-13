process ClarityGender {
    // Custom process to get gender status of sample from Clarity LIMS
    tag {"ClarityGender ${sample_id}"}
    label 'ClarityGender'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        tuple(sample_id, analysis_id)

    output:
        path(stdout emit: clarity_gender)

    script:
        """
        python ${clarity_epp_path}/xxxxget_gender.pyxxxx sample_id analysis_id \
        """
