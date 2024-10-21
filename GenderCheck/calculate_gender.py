#! /usr/bin/env python3

import argparse
import pysam


def is_valid_read(read, min_mapping_qual):
    """
    Check if a read is properly mapped and thus considered valid.

    Args:
        read (AlignedSegment object): represents a single read from BAM file along with its fields and optional tags.
        min_mapping_qual (int): Minimum mapping quality of reads to be considered.

    Returns:
        Boolean: Whether the read is valid and properly mapped.
    """
    if (read.mapping_quality >= min_mapping_qual and read.reference_end and read.reference_start):
        return True
    return False


def translate_gender(gender):
    """
    Translate gender using a map.

    Args:
        gender (str): Gender to translate, as stated in our LIMS.

    Returns:
        str: If input gender exists in map, the translation value is returned. Otherwise the input gender is returned.
    """
    translation = {"Man": "male", "Vrouw": "female", "Onbekend": "unknown", "unknown": "not_detected"}
    if gender in translation:
        return translation[gender]
    return gender



def validate_gender(gender):
    """
    Raise ValueError in case input value is not valid.

    Args:
        gender (str): Gender to validate.

    Raises:
        ValueError: If input gender is not allowed, a value error is raised. Its message includes the allowed options.
    """
    allowed_gender_options = ["male", "female"]
    allowed_unknown_options = ["unknown", "not_detected"]
    if gender not in allowed_gender_options and gender not in allowed_unknown_options:
        raise ValueError(f"Provided gender {gender} is not allowed. Should be one of {allowed_gender_options + allowed_unknown_options}.")


def get_gender_from_bam(bam, min_mapping_qual, locus_y, ratio_y):
    """
    Calculates the gender based on a bam file and certain quality thresholds.

    Args:
        bam (str): the path to the bam file
        min_mapping_qual (int): the mapping quality
        locus_y (str): the locus of the chromosome
        ratio_y (float): the ratio threshold

    Returns:
        str: 'female' in case observerd percentage of Y reads is lower or equal to 'ratio_y'. Otherwise 'male'.
    """
    with pysam.AlignmentFile(bam, "rb") as bam_file:
        y_reads = float(
                      sum([is_valid_read(read, min_mapping_qual) for read in bam_file.fetch(region=locus_y)])
                  )
        total_reads = float(bam_file.mapped)
        y_ratio_perc = (y_reads / total_reads) * 100
    if y_ratio_perc <= ratio_y:
        return "female"
    else:
        return "male"


def compare_and_evaluate_gender(measured_gender, stated_gender):
    """
    Compare and evaluate measured and stated gender.

    Args:
        measured_gender (str): Gender as observed in the bam file.
        stated_gender (str): Expected gender as stated in our LIMS / metadata.

    Returns:
        tuple of two string:
            str 'qc': Whether the comparison is evaluated as passed ('PASS') or failed ('FAIL').
            str 'message': A message to explain comparison evaluation result. An empty string in case evaluation passed.
    """
    if measured_gender == stated_gender or stated_gender == "unknown":  # if gender is unknown/onbekend in database, pass
        return "PASS", ""
    elif stated_gender == "not_detected":  # not_detected in database considered failed
        return "FAIL", f"Gender has value '{stated_gender}' in LIMS. Observed gender '{measured_gender}' could not be verified."
    else:
        return "FAIL", f"Stated gender {stated_gender} does not equal observed gender {measured_gender}."


def write_qc_file(sample_id, analysis_id, measured_gender, stated_gender, status, message, outputfolder):
    """
    Write gender comparison and evaluation to a qc file.

    Args:
        sample_id (str): Identifier of the sample to evaluate.
        analysis_id (str): Identifier of the analysis that created the BAM file used in the comparison.
        measured_gender (str): Gender as observed in the bam file.
        stated_gender (str): Expected gender as stated in our LIMS / metadata.
        status (str): Whether the comparison is evaluated as passed ('PASS') or failed ('FAIL').
        message (str): A message to explain comparison evaluation result. An empty string in case evaluation passed.
        outputfolder (str): File path to a directory where output files are saved.
    """
    with open(f"{outputfolder}/{sample_id}_{analysis_id}_gendercheck.txt", 'w') as write_file:
        write_file.write("sample_id\tanalysis_id\tmeasured_gender\tstated_gender\tstatus\tmessage\n")
        write_file.write(f"{sample_id}\t{analysis_id}\t{measured_gender}\t{stated_gender}\t{status}\t{message}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('sample_id', help='sample_id')
    parser.add_argument('analysis_id', help='analysis_id')
    parser.add_argument('bam', help='path to bam file')
    parser.add_argument('outputfolder', help='path to output folder')
    parser.add_argument('stated_gender', help='gender regarded as the truth')
    parser.add_argument(
        "ratio_y",
        type=float,
        help="maximunum chromosome Y ratio for females [float]"
    )
    parser.add_argument('min_mapping_qual', type=int, help='minimum mapping quality of reads to be considered [int]')
    parser.add_argument('locus_y', help='Coordinates for includes region on chromosome Y (chr:start-stop)')
    args = parser.parse_args()

    stated_gender = translate_gender(args.stated_gender)
    validate_gender(stated_gender)

    measured_gender = get_gender_from_bam(args.bam, args.min_mapping_qual, args.locus_y, args.ratio_y)
    validate_gender(measured_gender)

    qc, msg = compare_and_evaluate_gender(measured_gender, stated_gender)
    write_qc_file(args.sample_id, args.analysis_id, measured_gender, stated_gender, qc, msg, args.outputfolder)
