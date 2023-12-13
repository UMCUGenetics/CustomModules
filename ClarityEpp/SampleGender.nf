process SampleGender {
    // Custom process to run clarity_epp export sample_gender
    tag {"ClarityEpp SampleIndications ${sample_id}"}
    label 'ClarityEpp'
    label 'ClarityEpp_SampleGender'
    shell = ['/bin/bash', '-eo', 'pipefail']
    cache = false  //Disable cache to force a clarity export restarting the workflow.

    input:
        val(sample_id)

    output:
        tuple(val(sample_id), stdout)

    script:
        """
        source ${params.clarity_epp_path}/venv/bin/activate
        python ${params.clarity_epp_path}/clarity_epp.py export sample_indications \
        -a ${sample_id} | cut -f 2 | grep -v 'Indication' | tr -d '\n'
        """
}
