#!/usr/bin/env nextflow
nextflow.preview.dsl=2

process MosaicHunterGetGender {
    // Step 1: Process input files
    tag {"MosaicHunterGetGender ${sample_id}"}
    label 'MosaicHunterGetGender'
    container = 'quay.io/biocontainers/pysam:0.22.0--py38h15b938a_0'
    shell = ['/bin/bash', '-eo', 'pipefail']

    /*
    Define inputs.
    - Tuple consisting of a sample_id, a path to the .bam file, a path to the .bai file
    */
    input:
        tuple(val(sample_id), path(bam_files), path(bai_files))

    /*
    Define outputs.
    - A tuple containing respectively the number for the alpha and beta found in the sample.
    */
    output:
        path('gender_data_${sample_id}.tsv') into gender_data

    // The command to execute MosaicHunter Get Gender
    script:
    '''
    python ${projectDir}/CustomModules/MosaicHunter/1.0.0/get_gender_from_bam_chrx.py \
            ${sample_id} \
            ${bam_files} \
            ./ \
            ${mh_gender_mapping_qual} \
            $params.mh_gender_ratio_x_threshold_male \
            $params.mh_gender_ratio_x_threshold_female \
            $params.mh_gender_mapping_qual \
            $params.mh_gender_locus_x
    '''
}

process MosaicHunterQualityCorrection {
    // Step 1: Process input files
    tag {"MosaicHunterQualityCorrection ${sample_id}"}
    label 'MosaicHunterQualityCorrection'
    container = 'docker://umcugenbioinf/mosaic_hunter:1.0.0'
    shell = ['/bin/bash', '-euo', 'pipefail']

    /*
    Define inputs.
    - Tuple consisting of a sample_id, a path to the .bam file, a path to the .bai file
    - Path to the reference file
    - Path to the MosaicHunter common site filter bed file
    - Path to the MosaicHunter config file for the Quality Correction step
    */
    input:
        tuple(val(sample_id), path(bam_files), path(bai_files))
        path(mh_reference_file)
        path(mh_common_site_filter_bed_file)
        path(mh_config_file_one)
        path(gender_data) from gender_data

    /*
    Define outputs.
    - A tuple containing respectively the number for the alpha and beta found in the sample.
    */
    output:
        tuple(env(MHALPHA), env(MHBETA))

    // The command to execute MosaicHunter
    shell:
    '''
    SEX_STRING=$(awk 'NR>1 {print $2}' gender_data_!{sample_id})

    java -Xmx!{task.memory.toGiga()-4}G -jar /MosaicHunter/build/mosaichunter.jar \
-C !{mh_config_file_one} \
-P input_file=!{bam_files} \
-P mosaic_filter.sex=$SEX_STRING \
-P reference_file=!{mh_reference_file} \
-P common_site_filter.bed_file=!{mh_common_site_filter_bed_file} \
-P output_dir=./
    export MHALPHA="\$(grep -Po "(?<=alpha:\\s)\\w+" ./stdout*)"
    export MHBETA="\$(grep -Po "(?<=beta:\\s)\\w+" ./stdout*)"
    '''
}

process MosaicHunterMosaicVariantCalculation {
    // Caclulate the Mosaic Variants
    tag {"MosaicHunterMosaicVariantCalculation ${sample_id}"}
    label 'MosaicHunterMosaicVariantCalculation'
    container = 'docker://umcugenbioinf/mosaic_hunter:1.0.0'
    shell = ['/bin/bash', '-euo', 'pipefail']

    publishDir "QC/MosaicHunter", saveAs: { filename -> "${sample_id}_$filename" }, mode: 'copy'

    /*
    Define inputs.
    - Tuple consisting of a sample_id, a path to the .bam file, a path to the .bai file
    - Path to the reference file
    - Path to the MosaicHunter common site filter bed file
    - Path to the MosaicHunter config file for the Mosaic Variant Calculation step
    - The output of the MosaicHunterQualityCorrection step. This makes the environment variables available in this process
    - A tuple containing respectively the number for the alpha and beta found in the
      sample, which are stored in an environment variable.
    */
    input:
        tuple(val(sample_id), path(bam_files), path(bai_files))
        path(mh_reference_file)
        path(mh_common_site_filter_bed_file)
        path(mh_config_file_two)
        MosaicHunterQualityCorrection.out
        tuple(env(MHALPHA),env(MHBETA))
        path(gender_data) from gender_data

    // Final file, will be renamed to include sample_id and published to output directory
    output:
        path('final.passed.tsv')

    // The command to execute step two of MosaicHunter
    // First get the SEX_STRING from the sample
    shell:
    '''
    SEX_STRING=$(awk 'NR>1 {print $2}' gender_data_!{sample_id})

    java -Xmx!{task.memory.toGiga()-8}G -jar /MosaicHunter/build/mosaichunter.jar \
-C !{mh_config_file_two} \
-P mosaic_filter.alpha_param=$MHALPHA -P mosaic_filter.beta_param=$MHBETA \
-P input_file=!{bam_files} \
-P mosaic_filter.sex=$SEX_STRING \
-P reference_file=!{mh_reference_file} \
-P common_site_filter.bed_file=!{mh_common_site_filter_bed_file} \
-P output_dir=./
    '''
}