#! /usr/bin/env python3
# Import statements, alphabetic order of main package.
import argparse

# Third party libraries alphabetic order of main package.
import pysam

# Custom libraries alphabetic order of main package.
from CustomModules.Utils.get_gender import get_ratio_from_bam, get_gender_on_x_ratio


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

    ratio_x_reads = get_ratio_from_bam(args.bam, args.mapping_qual, args.locus_x)
    gender_data = get_gender_on_x_ratio(ratio_x_reads, args.ratio_x_threshold_male, args.ratio_x_threshold_female)

    write_genderdata_to_file(args.sample_id, gender_data, args.output_folder)
