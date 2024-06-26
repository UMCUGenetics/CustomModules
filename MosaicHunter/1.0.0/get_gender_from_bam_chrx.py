#! /usr/bin/env python3
import argparse
import pysam


def is_valid_read(read, mapping_qual):
    """Check if a read is properly mapped.

    Args:
        read (): a pysam.AlignmentFile read
        mapping_qual (int): the mapping quality

    Returns:
        bool: True if the read is valid
    """
    if (read.mapping_quality >= mapping_qual and read.reference_end and read.reference_start):
        return True
    return False


def get_gender_from_bam_chrx(bam_file_path, mapping_qual, locus_x, ratio_x_threshold_male, ratio_x_threshold_female):
    """
    Calculates the gender based on a bam file and thresholds

    Args:
        bam_file_path (str): the path to the bam file
        mapping_qual (int): the mapping quality
        locus_x (str): the locus of the chromosome
        ratio_x_threshold_male (float): the ratio threshold for males
        ratio_x_threshold_female (float): the ratio threshold for females

    Returns:
        A tuple that consist of two elements: the gender as a single character and if the gender was forced
    """
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        reads = float(
            sum([is_valid_read(read, mapping_qual) for read in bam_file.fetch(region=locus_x)])
        )
        total_reads = float(bam_file.mapped)
        ratio_perc = (reads / total_reads) * 100

    # Check ratios
    if ratio_perc <= ratio_x_threshold_male:
        return "M", False
    elif ratio_perc >= ratio_x_threshold_female:
        return "F", False
    else:
        # Force Female if unknown
        return "F", True


def write_genderdata_to_file(sample_id, gender_data, output_folder):
    """
    Write the gender data to a tsv file

    Args:
        sample_id (str): the id of the sample
        gender_data (tuple): the gender and the forced value as a tuple
        output_folder (str): the output folder
    """
    with open(f"{output_folder}/gender_data_{sample_id}.tsv", "w") as csv_file:
        gender, forced = gender_data
        csv_file.write("sample_id\tgender\tforced\n")
        csv_file.write(f"{sample_id}\t{gender}\t{forced}\n")


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

    gender_data = get_gender_from_bam_chrx(
        args.bam,
        args.mapping_qual,
        args.locus_x,
        args.ratio_x_threshold_male,
        args.ratio_x_threshold_female,
    )

    write_genderdata_to_file(args.sample_id, gender_data, args.output_folder)
