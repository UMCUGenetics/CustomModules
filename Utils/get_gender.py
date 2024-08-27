#! /usr/bin/env python3
# Import statements, alphabetic order of main package.

# Third party libraries alphabetic order of main package.
import pysam

# Custom libraries alphabetic order of main package.

def is_valid_read(read, mapping_qual):
    """
    Check if a read is properly mapped.

    Args:
        read (): a pysam.AlignmentFile read
        mapping_qual (int): the mapping quality

    Returns:
        bool: True if the read is valid
    """
    if (read.mapping_quality >= mapping_qual and read.reference_end and read.reference_start):
        return True
    return False


def get_ratio_from_bam(bam_file_path, mapping_qual, locus):
    """
    Calculates ratio of reads from given locus to total number of reads in bam file

    Args:
        bam_file_path (str): the path to the bam file
        mapping_qual (int): the mapping quality
        locus (str): the locus of the chromosome of interest

    Returns:
        A tuple that consist of two elements: the gender as a single character and if the gender was forced
    """
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        reads = float(
            sum([is_valid_read(read, mapping_qual) for read in bam_file.fetch(region=locus)])
        )
        total_reads = float(bam_file.mapped)
    return (reads / total_reads) * 100



def get_gender_on_x_ratio(ratio_x_perc, ratio_x_threshold_male, ratio_x_threshold_female):
    """
    Get gender as a single character based on retrieved x ratio and given male/female thresholds.

    Args:
        ratio_x_perc (float): the ratio of reads on locus X to all reads
        ratio_x_threshold_male (float): the ratio threshold for males
        ratio_x_threshold_female (float): the ratio threshold for females


    Returns:
        A tuple that consist of two elements: the gender as a single character and whether the gender was forced (boolean)
    """
    # Check ratios
    if ratio_x_perc <= ratio_x_threshold_male:
        return "M", False
    elif ratio_x_perc >= ratio_x_threshold_female:
        return "F", False
    else:
        # Force Female if unknown
        return "F", True


def get_gender_on_y_ratio(ratio_y_perc, ratio_y_threshold):
    """
    Get gender as female/male based on retrieved xyratio and given thresholds.

    Args:
        ratio_x_perc (float): the ratio of reads on locus Y to all reads
        ratio_y_threshold (float): the ratio threshold to use. A ratio lower or equal to threshold is determined to be female.

    Returns:
        str: 'female' if ratio is lower or equal to threshold, otherwise 'male'
    """
    if ratio_y_perc <= ratio_y_threshold:
        return "female"
    else:
        return "male"
