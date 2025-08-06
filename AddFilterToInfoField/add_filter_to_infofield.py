#! /usr/bin/env python3

import os
import sys

from collections import OrderedDict
import vcfpy

def add_info_header(reader, ID, Number, Type, Description):
    reader.header.add_info_line(
        vcfpy.OrderedDict([
            ('ID', ID),
            ('Number', Number),
            ('Type', Type),
            ('Description', Description)
        ])
    )
    return reader


def convert_vcfs(input_path, output_path):
    with open(input_path, 'r') as input_vcf, open(output_path, 'w') as output_vcf:
            reader = vcfpy.Reader.from_stream(input_vcf)

            # Add filterFlag to INFO field
            reader = add_info_header(reader, 'filterFlag', ".", 'String', "FILTER flags")

            # Write new header in new output
            writer = vcfpy.Writer.from_stream(output_vcf, reader.header)

            # Get sampleID (assuming single sampl VCF!)
            sample = reader.header.samples.names[0]

            for record in reader:
                #ADD FILTER to INFO field
                record.INFO['filterFlag'] = record.FILTER

                writer.write_record(record)

if __name__ == "__main__":
    convert_vcfs(sys.argv[1], sys.argv[2])
