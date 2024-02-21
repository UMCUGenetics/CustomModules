#! /usr/bin/env python
# Import statements, alphabetic order of main package.
import argparse
from errno import ENOENT as errno_ENOENT
from os import strerror as os_strerror
from pathlib import Path
import re
import sys
import warnings

# Third party libraries alphabetic order of main package.
from pandas import concat, DataFrame, merge, read_csv
import yaml


def non_empty_existing_path(file_or_dir):
    """This function checks whether the provided file or dir exists and is not empty.

    Args:
        file_or_dir (string): Input file or directory

    Raises:
        FileNotFoundError: If input string file_or_dir is neither a file nor a dir.
        OSError: If input is not a dir and file is empty.

    Returns:
        string: Provided input file or directory. If dir, suffix '/' might be added.
    """    
    input_path = Path(file_or_dir)
    if not input_path.is_file() and not input_path.is_dir():
        raise FileNotFoundError(errno_ENOENT, os_strerror(errno_ENOENT), file_or_dir)
    elif not input_path.is_dir() and not input_path.stat().st_size:
        raise OSError(f"File {file_or_dir} is empty.")
    elif input_path.is_dir() and file_or_dir[::-1][0] != "/":
        return f"{file_or_dir}/"
    return file_or_dir


def parse_arguments_and_check(args_in):
    """Parses arguments and validates / checks format of input.

    Args:
        args_in (list of strings): Commandline input arguments.

    Returns:
        Namespace: Convert argument strings to objects and assign them as attributes of the namespace.
    """    
    parser = argparse.ArgumentParser(
        description="Check and summarize sample quality using qc metrics and their thresholds."
    )
    parser.add_argument(
        "settings",
        type=non_empty_existing_path,
        help="QC settings specifying at least the filename, which qc_col, the threshold and the operator."
    )
    parser.add_argument(
        "output_path",
        type=non_empty_existing_path,
        help="QC check output path for all output files."
    )
    parser.add_argument(
        "output_prefix",
        type=str,
        help="QC check output prefix for all output file names."
    )
    parser.add_argument(
        "input_files", nargs="+", type=non_empty_existing_path,
        help="Input files containing a QC metric."
    )
    args = parser.parse_args(args_in)
    return args


def read_yaml(yaml_file):
    """Read input yaml

    Args:
        yaml_file (string): String with path to yaml file

    Raises:
        ValueError: If reading the file returns None (not recognized as YAML format).

    Returns:
        Object: Content of the YAML file.
    """    
    yaml_loaded = yaml.safe_load(open(yaml_file))
    if not yaml_loaded:
        raise ValueError("Could not load YAML.")
    return yaml.safe_load(open(yaml_file))


def check_allowed_operators(qc_operator):
    """Check if provided qc_operator is allowed.

    Args:
        qc_operator (string): (Custom / math) operator

    Raises:
        ValueError: If provided qc_operator is invalid / unsupported.
    """    
    operators = ["<", "<=", ">", ">=", "==", "!=", "match"]
    if qc_operator not in operators:
        raise ValueError(f"Unsupported operator provided: {qc_operator}. Please select from: {operators}")


def check_required_keys_metrics(qc_metrics):
    """Check if all required settings are included in the qc_settings

    Args:
        qc_metrics (list with dicts): qc settings per qc metric

    Raises:
        KeyError: Required key is not provided for the qc metric settings
    """    
    for req_key in ["filename", "qc_col", "threshold", "operator", "report_cols"]:
        if any([req_key not in setting.keys() for setting in qc_metrics]):
            raise KeyError(f"Required key {req_key} not in all metrics settings.")


def select_metrics(filename, input_files):
    """Using regular expression to match the qc metric filename with the input files

    Args:
        filename (string): Filename of qc metric, could be regex.
        input_files (list): All qc metrics input files.

    Returns:
        list: Input files matching the given filename.
    """
    # If filename is string, change into regex to match absolute and relative paths in input_files.
    if filename.isalpha():
        filename=".*" + filename
    metrics = list(filter(re.compile(f"{filename}").match, input_files))
    if not metrics:
        warnings.warn(UserWarning(f"No input file provided with filename pattern {filename}"))
        return None
    return metrics


def get_columns_to_report(qc_report_cols, qc_metric_cols, qc_col):
    """Get valid columns to include in final output report

    Args:
        qc_report_cols (list, string): column name(s) to include in report.
        qc_metric_cols (list, string): column name(s) in qc metric.
        qc_col (string): column name that contains the qc value/score.

    Raises:
        TypeError: qc_report_cols is neither a string nor list.
        ValueError: Provided column names in qc_report_cols does not exist in the qc metric.

    Returns:
        list of strings: Valid column names to include in report.
    """    
    not_existing_cols = list(set(qc_report_cols) - set(qc_metric_cols))
    if qc_report_cols == "@all":
        qc_report_cols = qc_metric_cols
    elif not isinstance(qc_report_cols, str) and not isinstance(qc_report_cols, list):
        raise TypeError(f"{qc_report_cols} not string, list or '@all'")
    elif not_existing_cols:
        raise ValueError(f"Some column names provided as report_cols do not exists: {not_existing_cols}")
    # Rename qc_col with qc_value
    qc_report_cols = list(map(lambda x: x.replace(qc_col, "qc_value"), qc_report_cols))
    # Add column qc_title
    qc_report_cols.insert(0, "qc_title")
    return qc_report_cols


def add_and_rename_columns(qc_metric, qc_title, qc_col, qc_operator, qc_threshold):
    """Add and rename columns in qc_metric.

    Args:
        qc_metric (pandas DataFrame): DataFrame with columns required to judge qc values
        qc_title (string): Title of the qc check
        qc_col (string): qc column name with qc value/score
        qc_operator (string): (Custom / math) operator
        qc_threshold (string, int or float): qc threshold

    Returns:
        Pandas DataFrame: DataFrame with qc metric.
    """
    qc_metric_assigned = qc_metric.assign(
        qc_title=qc_title.lower(),
        qc_status="PASS",
        qc_check=f"{qc_threshold} {qc_operator} {qc_col}",
        qc_msg="",
    )
    qc_metric_out = qc_metric_assigned.rename(columns={qc_col: "qc_value"})
    return qc_metric_out


def get_failed_rows(qc_metric, qc_col, qc_operator, qc_threshold):
    """Get rows that fail provided qc threshold 

    Args:
        qc_metric (pandas DataFrame): DataFrame with columns required to judge qc values
        qc_col (string): column name that contains the qc value/score.
        qc_operator (string): (Custom / math) operator
        qc_threshold (string, int or float): qc threshold

    Raises:
        TypeError: If qc threshold is neither 'match', str, int or float.

    Returns:
        object: The indexes of failed rows in qc_metric
    """    
    # Select failed rows using qc_threshold regex pattern and qc_operator 'match'
    if qc_operator == "match" and isinstance(qc_threshold, str):
        return qc_metric[qc_col].str.match(qc_threshold)
    elif isinstance(qc_threshold, str):  # Add quotes areound qc_threshold if it is string.
        # Note: using `query` has the advantage to dynamically build the comparison condition.
        # Disadvantage: no boolean indexing available. Assumed it is ok to use 'index'
        return qc_metric.query(f"`{qc_col}` {qc_operator} '{qc_threshold}'").index
    elif isinstance(qc_threshold, int) or isinstance(qc_threshold, float):  # Query failed_rows using integers/floats
        # Note: using `query` has the advantage to dynamically build the comparison condition.
        # Disadvantage: no boolean indexing available. Assumed it is ok to use 'index'
        return qc_metric.query(f"`{qc_col}` {qc_operator} {qc_threshold}").index
    else:
        raise TypeError(f"QC threshold {qc_threshold} type not supported.")


def add_failed_samples_metric(qc_metric, failed_rows, report_cols, sample_cols):
    """Failed samples are added to the output metric, and removed from qc_metric.

    Args:
        qc_metric (pandas DataFrame): DataFrame with columns required to judge qc values
        failed_rows (object): Object with indexes of failed rows in qc_metric
        report_cols (list): Valid column names (strings) to include in report.
        sample_cols (list): Columnames (strings) of sample names.

    Returns:
        qc_metric (DataFrame): DataFrame of qc metric without failed rows
        qc_metric_out (DataFrame): DataFrame of qc metric to report with failed rows
    """    
    qc_metric_out = DataFrame(columns=["sample", "qc_check", "qc_status", "qc_msg", "qc_value"])
    failed_samples = []
    if failed_rows.to_list():
        failed_samples = list(qc_metric.loc[failed_rows, sample_cols].values.ravel())
        qc_metric.loc[failed_rows, "qc_status"] = "FAIL"
        # Concatenate columns to report into a single column 'qc_msg'
        qc_metric["qc_msg"] = qc_metric.loc[failed_rows, report_cols].astype(str).apply(" ".join, axis=1)
        # Add failed samples to output
        # A single qc metric could have multiple sample columns
        # If a qc check fails for a 'multiple sample check', each individual sample is flagged as "failed"
        for sample_col in sample_cols:
            qc_metric_out = concat([
                qc_metric_out,
                (
                    qc_metric
                    .rename(columns={sample_col: "sample"})
                    .loc[failed_rows, qc_metric_out.columns.to_list()]
                    .groupby(["sample", "qc_check", "qc_status"], dropna=False)
                    .agg(lambda val: ';'.join(val.astype(str)))
                    .reset_index()
                )
            ])
        # Drop failed samples current metric
        for sample_col in sample_cols:
            drop_index = qc_metric[qc_metric[sample_col].isin(set(failed_samples))].index
            if drop_index.to_list():
                qc_metric.drop(drop_index, inplace=True)
    return qc_metric, qc_metric_out


def add_passed_samples_metric(qc_metric, qc_metric_out, sample_cols):
    """Passed samples are added to the output metric

    Args:
        qc_metric (DataFrame): DataFrame of qc metric without failed rows
        qc_metric_out (DataFrame): DataFrame of qc metric to report with failed rows
        sample_cols (list): Columnames (strings) of sample names.

    Returns:
        pandas DataFrame: Sorted qc metric for both passed and failed samples, without duplicates.
    """    
    # Add passed samples to output
    for sample_col in sample_cols:
        qc_metric_out = concat([
            qc_metric_out,
            (
                qc_metric
                .rename(columns={sample_col: "sample"})
                .loc[:, qc_metric_out.columns]
            )
        ])
    # Try to convert column qc_value to float.
    # If ValueError is raised, probably because column is a string, continue.    
    try:
        qc_metric_out["qc_value"] = qc_metric_out["qc_value"].astype("float")
    except ValueError:
        pass
    # In case 'multiple sample qc check',
    # output could contain duplicate rows for individual samples used in multiple comparisons.
    return qc_metric_out.sort_values(by=["qc_check", "qc_status"]).drop_duplicates(keep="first")


def create_and_write_output(qc_output, output_path, output_prefix):
    """Joined qc metrics is created and written to output file. 

    Args:
        qc_output (pandas DataFrame): Sorted judged qc metric for both passed and failed samples, without duplicates.
        output_path (string): Relative or absolute path where output should be saved. 
        output_prefix (string): Output prefix for output file.
    """    
    # Add qc_summary
    qc_output.insert(1, "qc_summary", "PASS")
    qc_output.loc[qc_output.isin(["FAIL"]).any(axis=1), "qc_summary"] = "FAIL"
    # Write summary output
    qc_output.to_csv(output_path + output_prefix + "_summary.csv", index=False, header=True)


def read_and_judge_metrics(qc, metrics):
    """Read and judge each single qc metric and join results. 

    Args:
        qc (dict): qc settings of the metric
        metrics (list): List of input files specific for single qc metric

    Returns:
        pandas DataFrame: Joined and judged qc metrics.
    """    
    for qc_file in metrics:
        qc_metric_raw = read_csv(qc_file, comment=qc.get("comment", None), delimiter=qc.get("delim", "\t"), quotechar='"')
        report_cols = get_columns_to_report(qc["report_cols"], qc_metric_raw.columns.to_list(), qc["qc_col"])
        qc_metric_edit = add_and_rename_columns(qc_metric_raw, qc["title"], qc["qc_col"], qc["operator"], qc["threshold"])
        failed_rows = get_failed_rows(qc_metric_edit, "qc_value", qc["operator"], qc["threshold"])
        qc_metric_subset, qc_metric_judged = add_failed_samples_metric(
            qc_metric_edit, failed_rows, report_cols, qc["sample_cols"]
            )
        qc_metric_judged = add_passed_samples_metric(qc_metric_subset, qc_metric_judged, qc["sample_cols"])
        # Rename columns
        suffix = f"_{qc['title'].lower()}"
        qc_judged_renamed = qc_metric_judged.add_suffix(suffix).rename(columns={f"sample{suffix}": "sample"})
        # Concatenate/merge metric output
        if "output" not in locals():  # First time
            output = qc_judged_renamed
        else:
            is_duplicate_sample = False
            # Check for duplicate sampleIDs before merge.
            if any(qc_judged_renamed["sample"].isin(output["sample"])):
                is_duplicate_sample = True
            output = merge(output, qc_judged_renamed, on=output.columns.tolist(), how="outer")
            if is_duplicate_sample:
                dup_sampleIDs = output[output['sample'].duplicated()]['sample'].to_list()
                # Duplicate sampleIDs with different column values
                if output["sample"].nunique() != output.shape[0]:
                    # Warning to parse all qc values / samples.
                    msg = f"Different qc values for duplicated sample IDs in input: {dup_sampleIDs}"
                # Duplicate sampleIDs same column values
                else:
                    msg = f"Sample IDs occur multiple times in input: {dup_sampleIDs}"
                warnings.warn(UserWarning(msg))
    return output


def check_qc(input_files, settings, output_path, output_prefix):
    """
    Main function to judge input files on configured qc settings. 
    It creates a single results table, each row representing 
        sample (string): sample name
        qc_summary: Summarized status of all qcs for sample (pass or fail)
        qc columns (5 per each qc metric);
            qc_check (string): QC check consiting of qc title, operator and threshold
            qc_status (string): Status of performed qc check (pass or fail)
            qc_msg (string): String with human readable message if sample failed qc check, empty if passed.  
            qc_value (string, int, float): qc value/score to check.

    Args:
        input_files (list): All qc metrics input files.
        settings (string): Path to yaml file
        output_path (string): Relative or absolute path where output should be saved. 
        output_prefix (string): Output prefix for output file.

    Raises:
        ValueError: No input files found to match any qc metric patterns defined in settings.
        ValueError: Duplicated samples with different values found in some of the input files. 
    """    
    # A single qc metric file can be used multiple times, by defining a metric section for each check in the qc settings.
    qc_settings = read_yaml(settings)
    check_required_keys_metrics(qc_settings["metrics"])
    duplicated_sample_file = []
    for qc_metric_settings in qc_settings["metrics"]:
        check_allowed_operators(qc_metric_settings["operator"])
        metric_files = select_metrics(qc_metric_settings["filename"], input_files)
        if not metric_files:
            continue
        # Join multiple metrices files into single table
        metric_out = read_and_judge_metrics(qc_metric_settings, metric_files)
        if any(metric_out.duplicated(subset="sample")):
            duplicated_sample_file.append(qc_metric_settings["filename"])
            continue
        if "merged_out" not in locals():
            merged_out = metric_out
        else:
            # Join all metrics output to single table.
            merged_out = merge(merged_out, metric_out, on="sample", how="outer")

    if "metric_out" not in locals():
        raise ValueError("No input files found to match any qc metric pattern.")
    if duplicated_sample_file:
        raise ValueError(f"Duplicated samples with different values found in files matching {duplicated_sample_file}.")
    create_and_write_output(merged_out, output_path, output_prefix)


if __name__ == "__main__":
    args = parse_arguments_and_check(args_in=sys.argv[1:])
    check_qc(args.input_files, args.settings, args.output_path, args.output_prefix)
