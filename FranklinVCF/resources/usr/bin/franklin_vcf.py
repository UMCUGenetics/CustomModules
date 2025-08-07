#!/usr/bin/env -S uv run --script --no-cache
# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "click",
#   "vcfpy",
# ]
# ///

import click
import vcfpy


def add_info_header(reader, ID, Number, Type, Description):
    reader.header.add_info_line(
        vcfpy.OrderedDict([("ID", ID), ("Number", Number), ("Type", Type), ("Description", Description)])
    )
    return reader


@click.command()
@click.argument("input_vcf", type=click.File("r"))
@click.argument("output_vcf", type=click.File("w"))
def main(input_vcf, output_vcf):
    reader = vcfpy.Reader.from_stream(input_vcf)

    # Add filterFlag to INFO field
    reader = add_info_header(reader, "filterFlag", ".", "String", "FILTER flags")

    # Write new header in new output
    writer = vcfpy.Writer.from_stream(output_vcf, reader.header)

    for record in reader:
        # ADD FILTER to INFO field
        record.INFO["filterFlag"] = record.FILTER

        writer.write_record(record)


if __name__ == "__main__":
    main()
