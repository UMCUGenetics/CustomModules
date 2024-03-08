#! /usr/bin/env python3
import argparse
import pysam


def is_valid_read(read, mapping_qual):
    """Check if a read is properly mapped."""
    if (read.mapping_quality >= mapping_qual and read.reference_end and read.reference_start):
        return True
    return False


def get_gender_from_bam_chrx(bam_file_path, mapping_qual, locus_x, ratio_x_threshold_male, ratio_x_threshold_female):
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        reads = float(
            sum([is_valid_read(read, mapping_qual) for read in bam_file.fetch(region=locus_x)])
        )
        total_reads = float(bam_file.mapped)
        ratio_perc = (reads / total_reads) * 100

    # Check ratios
    if ratio_perc <= ratio_x_threshold_male:
        return "M"
    elif ratio_perc >= ratio_x_threshold_female:
        return "F"
    else:
        # Force Female if unknown
        return "F"


def write_to_file(sample_id, gender, output_folder):
    with open(f"{output_folder}/gender_data_{sample_id}.tsv", "w") as csv_file:
        csv_file.write("sample_id\tgender\n")
        csv_file.write(f"{sample_id}\t{gender}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('sample_id', help='sample_id')
    parser.add_argument('bam', help='path to bam file')
    parser.add_argument('output_folder', help='path to output folder')
    parser.add_argument(
        "ratio_x_threshold_male",
        type=float,
        help="maximunum chromosome X ratio for males [float]"
    )
    parser.add_argument(
        "ratio_x_threshold_female",
        type=float,
        help="maximunum chromosome X ratio for females [float]"
    )
    parser.add_argument('mapping_qual', type=int, help='minimum mapping quality of reads to be considered [int]')
    parser.add_argument('locus_x', help='Coordinates for includes region on chromosome Y (chr:start-stop)')
    args = parser.parse_args()

    gender = get_gender_from_bam_chrx(
        args.bam,
        args.mapping_qual,
        args.locus_x,
        args.ratio_x_threshold_male,
        args.ratio_x_threshold_female,
    )

    write_to_file(args.sample_id, gender, args.output_folder)
