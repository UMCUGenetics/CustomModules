#! /usr/bin/env python
import argparse


def parse_ped(ped_file):
    """
    Parse ped file to a samples dict where per sample (key) a metadata dict (value) is created
    with the family ID (string), parents (list) and children (list) as content.

    Args:
        ped_file (open file): Open file object to a table with sample metadata, including the following values:
            familyID (string): Unique identifier per family. In other words, s
                samples from the same family will have the same familyID.
            sampleID (string): Unique identifier of sample.
            father (string): reference to another sampleID.
            mother (string): reference to another sampleID.
            sex (int): 0 (unknown), 1 (male) or 2 (female)
            phenotype (int): 0 (unknown), 1 (unaffected) or 2 (affected)

    Returns:
        dict: Per sample (key) a metadata dict (value) with
            the family ID (string), parents (list) and children (list).
    """
    samples = {}  # 'sample_id': {'family': 'fam_id', 'parents': ['sample_id', 'sample_id']}
    for line in ped_file:
        ped_data = line.strip().split()
        family, sample, father, mother, sex, phenotype = ped_data

        # Create samples
        if sample not in samples:
            samples[sample] = {'family': family, 'parents': [], 'children': []}
        if father != '0' and father not in samples:
            samples[father] = {'family': family, 'parents': [], 'children': []}
        if mother != '0' and mother not in samples:
            samples[mother] = {'family': family, 'parents': [], 'children': []}

        # Save sample relations
        if father != '0':
            samples[sample]['parents'].append(father)
            samples[father]['children'].append(sample)
        if mother != '0':
            samples[sample]['parents'].append(mother)
            samples[mother]['children'].append(sample)
    return samples


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Check kinship output based on ped file.')
    parser.add_argument('ped_file', type=argparse.FileType('r'), help='PED file')
    parser.add_argument('samples_analysis', nargs='+', help='samples within the analysis (space seperated)')
    arguments = parser.parse_args()

    samples = parse_ped(arguments.ped_file)
    trio_sample = []
    for sample in samples:
        if len(samples[sample]['parents']) == 2:  # Sample = child with parents
            for sample_analysis in list(set(arguments.samples_analysis)):
                if sample in sample_analysis:
                    trio_sample.append(sample)

    print(",".join(trio_sample))
