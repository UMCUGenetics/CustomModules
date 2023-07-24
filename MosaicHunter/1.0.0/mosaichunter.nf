#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// Step 1: Process input files
process MosaicHunterStepOne {
    container = 'docker://anabuurs/mosaichunter'
    label = 'MosaicHunterStepOne'
    tag = tag {"MosaicHunterStepOne $sample_id"}

    /*
    Define inputs.
    - input_bam_file -> a path to the .bam file
    - input_bai_file -> a path to the .bam.bai file
    - input_sex_string -> a string (must be "M" or "F") for the sex of the patient
    */
    input:
    tuple(sample_id, path(bam_files), path(bai_files), path(mh_reference_file), path(mh_common_site_filter_bed_file), path(mh_config_file))

    /*
    Define outputs.
    - a tuple containing respectively the number for the alpha and beta found in the sample.
    */
    output:
    tuple(env(MHALPHA), env(MHBETA))

    // The command to execute MosaicHunter
    """
    SEX_STRING=$(echo "$sample_id" | egrep -o '[MF]')

    java -Xmx${task.memory.toGiga()-4}G -jar /MosaicHunter/build/mosaichunter.jar \
-C $mh_config_file \
-P input_file=$input_bam_file \
-P mosaic_filter.sex=\$SEX_STRING \
-P reference_file=$mh_reference_file \
-P common_site_filter.bed_file=$mh_common_site_filter_bed_file \
-P output_dir=./
    export MHALPHA="\$(grep -Po "(?<=beta:\\s)\\w+" ./stdout*)"
    export MHBETA="\$(grep -Po "(?<=beta:\\s)\\w+" ./stdout*)"
    """
}

process MosaicHunterStepTwo {
    container = 'docker://anabuurs/mosaichunter'
    label = 'MosaicHunterStepTwo'
    tag = tag {"MosaicHunterStepTwo $sample_id"}

    /*
    Define inputs.
    - input_bam_file -> a path to the .bam file
    - input_bai_file -> a path to the .bam.bai file
    - input_sex_string -> a string (must be "M" or "F") for the sex of the patient
    - tuple val(env(MHALPHA),env(MHBETA)) -> a tuple containing respectively the number for the alpha and beta found in the sample, which are stored in an enivroment variable.
    */
    input:
    tuple(sample_id, path(bam_files), path(bai_files))
    tuple env(MHALPHA),env(MHBETA)
    path(mh_reference_file)
    path(mh_common_site_filter_bed_file)
    path(mh_config_file)
    MosaicHunterStepOne.out

    // Final file, will be published to output directory
    output:
    file 'final.passed.tsv'

    // The command to execute step two of MosaicHunter
    """
    SEX_STRING=$(echo "$sample_id" | egrep -o '[MF]')

    java -Xmx${task.memory.toGiga()-8}G -jar /MosaicHunter/build/mosaichunter.jar \
-C $mh_config_file \
-P mosaic_filter.alpha_param=\$MHALPHA -P mosaic_filter.beta_param=\$MHBETA \
-P input_file=$bam_files \
-P mosaic_filter.sex=\$SEX_STRING \
-P reference_file=$mh_reference_file \
-P common_site_filter.bed_file=$mh_common_site_filter_bed_file \
-P output_dir=./
    """
}