process MipsTrimDedup {
    // Custom process to run MIPS TrimDedup
    tag {"MIPS TrimDedup ${sample_id} - ${rg_id}"}
    label 'MIPS_1_0_1'
    label 'MIPS_1_0_1_TrimDedup'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        tuple(sample_id, rg_id, path(r1_fastqs), path(r2_fastqs))

    output:
        tuple(sample_id, rg_id, path('*_LMergedTrimmedDedup_R1_*.fastq.gz'), path('*_LMergedTrimmedDedup_R2_*.fastq.gz'), emit: fastq_files)

    script:
        def r1_args = r1_fastqs.collect{ "$it" }.join(" ")
        def r2_args = r2_fastqs.collect{ "$it" }.join(" ")

        rg_id = "${sample_id}_MergedTrimmedDedup"

        """
        python ${params.mips_trim_dedup_path}/mips_trim_dedup.py -d ${params.dxtracks_path}/${params.mips_design_file}  -l ${params.mips_uuid_length} -ur ${params.mips_uuid_read} -r1 ${r1_args} -r2 ${r2_args}
        """
}
