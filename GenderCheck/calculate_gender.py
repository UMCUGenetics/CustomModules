#! /usr/bin/env python3

import argparse
import pysam


def is_valid_read(read, mapping_qual):
    """Check if a read is properly mapped."""
    if (read.mapping_quality >= mapping_qual and read.reference_end and read.reference_start):
        return True
    return False


def get_gender_from_bam(bam, mapping_qual, locus_y, ratio_y):
    with pysam.AlignmentFile(bam, "rb") as bam_file:
        y_reads = float(
                      sum([is_valid_read(read, mapping_qual) for read in bam_file.fetch(region=locus_y)])
                  )
        total_reads = float(bam_file.mapped)
        y_ratio_perc = (y_reads / total_reads) * 100
    if y_ratio_perc <= ratio_y:
        return "female"
    else:
        return "male"


def compare_gender(test_gender, true_gender):
    if test_gender == true_gender or true_gender == "unknown":  # if gender is unknown/onbekend in database, pass
        qc = "PASS"
        msg = ""
    else:  # not_detected in database considered failed
        qc = "FAIL"
        msg = f"True gender {true_gender} does not equal estimated gender {test_gender}."
    return qc, msg


def write_qc_file(sample_id, analysis_id, test_gender, true_gender, status, message, outputfolder):

    with open(f"{outputfolder}/{sample_id}_{analysis_id}_gendercheck.txt", 'w') as write_file:
        write_file.write("sample_id\tanalysis_id\ttest_gender\ttrue_gender\tstatus\tmessage\n")
        write_file.write(f"{sample_id}\t{analysis_id}\t{test_gender}\t{true_gender}\t{status}\t{message}\n")


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

    test_gender = get_gender_from_bam(args.bam, args.mapping_qual, args.locus_y, args.ratio_y)
    qc, msg = compare_gender(test_gender, true_gender)
    write_qc_file(args.sample_id, args.analysis_id, test_gender, true_gender, qc, msg, args.outputfolder)
