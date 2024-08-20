#! /usr/bin/env python
# Import statements, alphabetic order of main package.
import argparse
from csv import writer
from errno import ENOENT as errno_ENOENT
from os import strerror as os_strerror
from pathlib import Path
from sys import argv, stdout
import tempfile

# Third party libraries alphabetic order of main package.
from pandas import DataFrame, read_table

# TODO: add docstrings
# TODO: add pytest

def non_empty_existing_path(file_or_dir):
    input_path = Path(file_or_dir)
    if not input_path.is_file() and not input_path.is_dir():
        raise FileNotFoundError(errno_ENOENT, os_strerror(errno_ENOENT), file_or_dir)
    elif not input_path.is_dir() and not input_path.stat().st_size:
        raise OSError(f'File {file_or_dir} is empty.')
    elif input_path.is_dir() and file_or_dir[::-1][0] != '/':
        return f'{file_or_dir}/'
    return file_or_dir


def parse_arguments_and_check(args_in):
    parser = argparse.ArgumentParser(description='Check kinship output based on ped file.')
    parser.add_argument('kinship_file', type=argparse.FileType('r'), help='Kinship file')
    parser.add_argument('ped_file', type=argparse.FileType('r'), help='PED file')
    parser.add_argument(
        '-p', '--output_path', type=non_empty_existing_path, default=None,
        help='Kinship output path where output file will be stored.',
    )
    parser.add_argument(
        '-o', '--output_prefix', type=str, default=None,
        help='Kinship output prefix for all output file names.'
    )
    parser.add_argument(
        '-s', '--kinship_settings', type=float, nargs=2, metavar=('minimum', 'maximum'), default=[0.177, 0.354],
        help='Kinship settings defining minimum and maximum threshold.'
    )
    arguments = parser.parse_args(args_in)
    return arguments


def parse_ped(ped_file):
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


def read_kinship(kinship_file, kinship_min, kinship_max):
    kinship_data = (
        # Read kinship data specific columns
        read_table(kinship_file, delimiter='\t', usecols=['FID1', 'FID2', 'Kinship'])
        # Rename columns
        .rename(columns={'FID1': 'sample_1', 'FID2': 'sample_2', 'Kinship': 'kinship'})
        # Add columns with default values
        .assign(
            related=None, type=None, status=None,
            thresholds=f'{kinship_min},{kinship_max}',
            message=''
        )
    )
    return kinship_data


def check_and_annotate_kinship(kinship_data, samples, kinship_min, kinship_max):
    for index, row in kinship_data.iterrows():
        status = 'OK'
        message = ''
        # Related
        if samples[row.sample_1]['family'] == samples[row.sample_2]['family']:
            related = True
            # Parent - child
            if row.sample_2 in samples[row.sample_1]['parents'] or row.sample_1 in samples[row.sample_2]['parents']:
                type = 'parent_child'
                if row.kinship <= kinship_min or row.kinship >= kinship_max:
                    status = 'FAIL'
                    expected_value_range = f'> {kinship_min} and < {kinship_max}'
            # Parent - Parent -> both samples have the same children
            elif samples[row.sample_1]['children'] and samples[row.sample_1]['children'] == samples[row.sample_2]['children']:
                type = 'parent_parent'
                if row.kinship > kinship_min:
                    status = 'FAIL'
                    expected_value_range = f'<= {kinship_min}'
            # Assume siblings
            else:
                type = 'sibling_sibling'
                if row.kinship <= kinship_min or row.kinship >= kinship_max:
                    status = 'FAIL'
                    expected_value_range = f'> {kinship_min} and < {kinship_max}'
        # Unrelated
        else:
            related = False
            type = 'unrelated'
            if row.kinship > kinship_min:
                status = 'FAIL'
                expected_value_range = f'<= {kinship_min}'

        # Create end user message if status is fail
        if status == 'FAIL':
            message = (
                f'Kinship value {row.kinship} between '
                f'{row.sample_1} ({samples[row.sample_1]['family']}) and {row.sample_2} ({samples[row.sample_2]['family']}) '
                f'is not between expected values for {type}: {expected_value_range}'
            )
        # Update row with retrieved related (boolean), relationship type, status (OK / FAIL) and message.
        kinship_data.loc[index, ['related', 'type', 'status', 'message']] = related, type, status, message
    return kinship_data


def write_kinship(df_kinship_out, output_path, output_prefix):
    # Collect all comments
    comments = []
    if any(df_kinship_out.status == 'FAIL'):
        comments.append('# WARNING: Kinship errors found.\n')
    else:
        comments.append('# No kinship errors found.\n')
    # Assume all row values of column thresholds are the same.
    comments.append(f'# Used kinship check settings: {df_kinship_out.thresholds[0]}\n')

    # Write to provided output settings or to a tempfile
    if output_path and output_prefix:
        file_out = open(f'{output_path}/{output_prefix}.kinship_check.out', 'a+')
    else:
        file_out = tempfile.TemporaryFile(mode='a+')

    # Write comments as header
    file_out.writelines(comments)
    # Append annotated kinship results
    df_kinship_out.to_csv(file_out, sep='\t', index=False, header=True)
    # Decide if file should be printed to stdout instead
    if not output_path or not output_prefix:
        file_out.seek(0)
        print(file_out.read())
    # Closing a tempfile will delete it as well
    file_out.close()


if __name__ == '__main__':
    arguments = parse_arguments_and_check(args_in=argv[1:])
    kinship_min, kinship_max = arguments.kinship_settings

    samples = parse_ped(arguments.ped_file)
    df_kinship_in = read_kinship(arguments.kinship_file, kinship_min, kinship_max)
    df_kinship_out = check_and_annotate_kinship(df_kinship_in, samples, kinship_min, kinship_max)
    write_kinship(df_kinship_out, arguments.output_path, arguments.output_prefix)
