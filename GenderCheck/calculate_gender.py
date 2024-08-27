#! /usr/bin/env python3
# Import statements, alphabetic order of main package.
import argparse

# Third party libraries alphabetic order of main package.
import pysam

# Custom libraries alphabetic order of main package.
from CustomModules.Utils.get_gender import get_ratio_from_bam, get_gender_on_y_ratio


def compare_gender(sample_id, analysis_id, test_gender, true_gender):
    # PASS qc if test and true gender are the same or if true gender is unknown/onbekend
    if test_gender == true_gender or true_gender == "unknown":
        qc = "PASS"
    # FAIL qc otherwise. not_detected as true gender is considered failed
    else:
        qc = "FAIL"
    return f"{sample_id}\t{analysis_id}\t{test_gender}\t{true_gender}\t{qc}\n"


def write_qc_file(sample_id, analysis_id, comparison, outputfolder):
    with open(f"{outputfolder}/{sample_id}_{analysis_id}_gendercheck.txt", 'w') as write_file:
        write_file.write("sample_id\tanalysis_id\ttest_gender\ttrue_gender\tstatus\n")
        write_file.write(comparison)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('sample_id', help='sample_id')
    parser.add_argument('analysis_id', help='analysis_id')
    parser.add_argument('bam', help='path to bam file')
    parser.add_argument('outputfolder', help='path to output folder')
    parser.add_argument('true_gender', help='gender regarded as the truth')
    parser.add_argument(
        "ratio_y",
        type=float,
        help="maximunum chromosome Y ratio for females [float]"
    )
    parser.add_argument('mapping_qual', type=int, help='minimum mapping quality of reads to be considered [int]')
    parser.add_argument('locus_y', help='Coordinates for includes region on chromosome Y (chr:start-stop)')
    args = parser.parse_args()

    translation = {"Man": "male", "Vrouw": "female", "Onbekend": "unknown", "unknown": "not_detected"}
    true_gender = args.true_gender
    if true_gender in translation:
        true_gender = translation[true_gender]

    ratio_y_reads = get_ratio_from_bam(args.bam, args.mapping_qual, args.locus_y)
    test_gender = get_gender_on_y_ratio(ratio_y_reads, args.ratio_y)
    comparison = compare_gender(args.sample_id, args.analysis_id, test_gender, true_gender)
    write_qc_file(args.sample_id, args.analysis_id, comparison, args.outputfolder)
