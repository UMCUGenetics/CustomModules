process CheckFingerprintVCF {
    // Custom process to check fingerprint vcf files
    tag {"CheckFingerprintVCF"}
    label 'CheckFingerprintVCF'
    container = 'docker://python:3.9.6'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        path(vcf_files)

    output:
        tuple(path('disapprovedVCFs'), path('approvedVCFs/*.vcf'), emit: vcf_files)
        path('logbook.txt', emit: logbook)


    script:
        """
        python ${baseDir}/CustomModules/CheckFingerprintVCF/check_fingerprint_vcf.py ${vcf_files} > logbook.txt
        """
}
