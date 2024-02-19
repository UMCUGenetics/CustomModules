#!/usr/bin/env nextflow
nextflow.preview.dsl=2

process MosaicHunterQualityCorrection {
    // Step 1: Process input files
    container = 'docker://umcugenbioinf/mosaic_hunter:1.0.0'
    label = 'MosaicHunterQualityCorrection'
    tag = "MosaicHunterQualityCorrection $sample_id"

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

    /*
    Define outputs.
    - A tuple containing respectively the number for the alpha and beta found in the sample.
    */
    output:
        tuple(env(MHALPHA), env(MHBETA))

    // The command to execute MosaicHunter
    shell:
    '''
    SEX_STRING=$(echo "!{sample_id}" | egrep -o '[MF]')

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
    container = 'docker://umcugenbioinf/mosaic_hunter:1.0.0'
    label = 'MosaicHunterMosaicVariantCalculation'
    tag = "MosaicHunterMosaicVariantCalculation $sample_id"

    publishDir "QC/MosaicHunter", saveAs: { filename -> "${sample_id}_$filename" }, mode: 'copy'

    /*
    Define inputs.
    - Tuple consisting of a sample_id, a path to the .bam file, a path to the .bai file
    - Path to the reference file
    - Path to the MosaicHunter common site filter bed file
    - Path to the MosaicHunter config file for the Mosaic Variant Calculation step
    - The output of the MosaicHunterQualityCorrection step. This makes the enviroment variables available in this process
    - A tuple containing respectively the number for the alpha and beta found in the
      sample, which are stored in an enivroment variable.
    */
    input:
        tuple(val(sample_id), path(bam_files), path(bai_files))
        path(mh_reference_file)
        path(mh_common_site_filter_bed_file)
        path(mh_config_file_two)
        MosaicHunterQualityCorrection.out
        tuple(env(MHALPHA),env(MHBETA))

    // Final file, will be renamed to include sample_id and published to output directory
    output:
        path('final.passed.tsv')

    // The command to execute step two of MosaicHunter
    shell:
    '''
    SEX_STRING=$(echo "!{sample_id}" | egrep -o '[MF]')

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