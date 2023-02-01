import argparse
import os
from pathlib import Path
import shutil


def valid_outputfile(file):
    path_file = Path(file)
    if path_file.is_file():
        print(Warning("Output file already exists. Output will be added."))
    return open(file, "a")


def parse_arguments_and_check():
    parser = argparse.ArgumentParser(description='Check fingerprint vcf files.')
    parser.add_argument(
        '-o', '--output', type=valid_outputfile, help='Output filename to write results to.', default='logbook.txt'
    )
    parser.add_argument('fingerprint_vcf_files', type=argparse.FileType('r'), nargs='+', help='Fingerprint VCF')
    arguments = parser.parse_args()
    return arguments


def parse_vcf(vcf_file):
    refcov = 0  # Total reference reads for all homozygous (1/1) calls
    totalcov = 0  # Total coverage for all homozygous calls
    homaltcount = 0  # Number of homozygous calls
    ycount = 0  # Sum of coverage for two Y SNPs
    lowcovcount = 0  # Number of SNPs with <15X coverage
    disbalancecount = 0  # Number of heterozygous (0/1) calls with allelefrequency <0.2 or >0.8

    for line in vcf_file:
        if not line.startswith('#'):
            line = line.split()

            # Parse Genotype format
            gt_format = line[8].split(':')
            gt_index = gt_format.index('GT')

            # Parse sample genotype
            gt_values = line[9].split(':')
            gt_value = gt_values[gt_index]

            if line[0] == 'Y':
                if gt_value != './.':
                    ycount += int(gt_values[gt_format.index('DP')])
            elif gt_value == '1/1':
                homaltcount += 1
                if int(gt_values[gt_format.index('DP')]) < 15:
                    lowcovcount += 1
                refcov += int(gt_values[gt_format.index('AD')].split(',')[0])
                totalcov += int(gt_values[gt_format.index('DP')])
            elif gt_value != '1/1':
                if gt_value == './.':
                    lowcovcount += 1
                elif gt_value == '0/0':
                    if int(gt_values[gt_format.index('DP')]) < 15:
                        lowcovcount += 1
                else:
                    if int(gt_values[gt_format.index('DP')]) < 15:
                        lowcovcount += 1
            if gt_value == '0/1':
                af_value = int(gt_values[gt_format.index('AD')].split(',')[0]) / \
                    float(int(gt_values[gt_format.index('DP')]))
                if af_value > 0.8 or af_value < 0.2:
                    disbalancecount += 1
    contamination = round(refcov / float(totalcov), 6)
    return vcf_file.name, lowcovcount, homaltcount, contamination, vcf_file.name[8], ycount, disbalancecount


def get_state_and_warning(lowcovcount, disbalancecount, homaltcount, vcf_file_gender, ycount, vcf_file_name, gender):
    states = {
        'approved': 'approved',
        'disapproved': 'disapproved',
        'discuss': 'discuss with lab and disapprove if needed',
    }
    state = states['approved']
    warning = None
    if lowcovcount > 15:
        state = states['disapproved']
        warning = f'>15 SNPs with <15X coverage ({lowcovcount})'
    elif disbalancecount > 8:
        state = states['disapproved']
        warning = f'>8 heterozygous SNPs with <20% MAF ({disbalancecount})'
    elif homaltcount < 8:
        state = states['disapproved']
        warning = f'<8 homozygous ALT SNPs called ({homaltcount})'
    elif (vcf_file_gender == 'F' and ycount > 100) or (vcf_file_gender == 'M' and ycount < 100):
        state = states['discuss']
        warning = f'gender {vcf_file_gender} with {ycount} reads on chromosome Y'
    if vcf_file_gender not in gender:
        gender_warning = f'gender used in filename not allowed: {vcf_file_name}'
        if warning:
            state = ';'.join((state, states['discuss']))
            warning = ';'.join((warning, gender_warning))
        else:
            state = states['discuss']
            warning = gender_warning
    return state, warning


def write_result(output_file, lst_values):
    output_file.write('\t'.join(map(str, lst_values)) + '\n')


def check_fingerprint_vcf(output_file, fingerprint_vcf_files):
    gender = ['M', 'F', 'O']
    # Add header to output file
    disapproved_folder = 'disapprovedVCFs'
    approved_folder = 'approvedVCFs'
    os.mkdir(disapproved_folder)
    os.mkdir(approved_folder)
    write_result(
        output_file,
        [
            'filename', 'lowcovcount', 'homaltcount', 'contamination', 'vcf_file_gender', 'ycount', 'disbalancecount',
            'state', 'warning'
        ]
    )
    for vcf_file in fingerprint_vcf_files:
        vcf_file_name, lowcovcount, homaltcount, contamination, vcf_file_gender, ycount, disbalancecount = parse_vcf(vcf_file)
        state, warning = get_state_and_warning(
            lowcovcount, disbalancecount, homaltcount, vcf_file_gender, ycount, vcf_file_name, gender
        )
        # Move vcf file
        if state == 'disapproved':
            shutil.move(vcf_file_name, disapproved_folder)
        else:
            shutil.move(vcf_file_name, approved_folder)
        write_result(
            output_file,
            [vcf_file_name, lowcovcount, homaltcount, contamination, vcf_file_gender, ycount, disbalancecount, state, warning]
        )
    output_file.close()


if __name__ == '__main__':
    arguments = parse_arguments_and_check()
    check_fingerprint_vcf(output_file=arguments.output, fingerprint_vcf_files=arguments.fingerprint_vcf_files)
