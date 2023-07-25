#! /usr/bin/env python

import argparse

import vcfpy


def vcf2csv(args):
    reader = vcfpy.Reader.from_stream(args.vcf_file)

    # print header
    print('chrom', 'pos', 'id', 'genotype', 'sample', 'sequencing_run', sep=',')

    for record in reader:
        sample_call = record.calls[0]  # Assume single sample vcf
        if args.min_dp and sample_call.data['DP'] < args.min_dp:
            continue

        if args.min_gq and sample_call.is_variant and sample_call.data['GQ'] < args.min_gq:
            continue

        if args.min_rgq and not sample_call.is_variant and sample_call.data['RGQ'] < args.min_rgq:
            continue

        # Set genotype separator
        gt_sep = '/'
        if sample_call.is_phased:
            gt_sep = '|'

        print(
            record.CHROM,
            record.POS,
            record.ID[0] if record.ID else '',  # print record.ID unless empty
            gt_sep.join(sample_call.gt_bases),  # create genotype using genotype separator
            sample_call.sample,
            args.sequencing_run_id,
            sep=','
        )


if __name__ == "__main__":
    # Argument parser
    parser = argparse.ArgumentParser(description="")
    parser.set_defaults(func=vcf2csv)
    # Required arguments
    parser.add_argument("sequencing_run", help="")
    parser.add_argument("vcf_file", type=argparse.FileType('r', encoding='latin-1'), help="")
    # Optional arguments
    parser.add_argument("--min_dp", type=int, help="Minimum DP")
    parser.add_argument("--min_gq", type=int, help="Minimum GQ")
    parser.add_argument("--min_rgq", type=int, help="Minimum RGQ")

    args = parser.parse_args()
    args.func(args)
