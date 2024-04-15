#!/usr/bin/env python
# Import statements, alphabetic order of main package.
from itertools import combinations
from pathlib import Path
import pytest
from pytest_unordered import unordered

# Third party libraries alphabetic order of main package.
from pandas import DataFrame, read_csv

# Custom libraries alphabetic order of main package.
import check_qc


@pytest.fixture(scope="module", autouse=True)
def setup_test_path(tmp_path_factory):
    test_tmp_path = str(tmp_path_factory.mktemp("test")) + "/"
    # Create empty files
    open(str(test_tmp_path) + "/empty.txt", "a").close()
    open(str(test_tmp_path) + "/empty.yaml", "a").close()
    return test_tmp_path


@pytest.fixture(scope="class")
def mock_settings(class_mocker, dataset_class):
    return class_mocker.patch("check_qc.read_yaml", return_value=dataset_class["settings_single_metric"])


class TestNonEmptyExistingPath():
    def test_existing_dir(self, setup_test_path):
        file_or_dir = check_qc.non_empty_existing_path(setup_test_path)
        assert file_or_dir

    def test_not_file_not_dir(self):
        fake_string = "fake_string"
        with pytest.raises(FileNotFoundError) as file_dir_error:
            check_qc.non_empty_existing_path(fake_string)
        assert fake_string in str(file_dir_error.value)

    def test_empty_file(self, setup_test_path):
        with pytest.raises(OSError) as empty_error:
            check_qc.non_empty_existing_path(setup_test_path + "empty.txt")
        assert f"File {setup_test_path}empty.txt is empty." in str(empty_error.value)

    def test_append_suffix(self, setup_test_path):
        dir_without_suffix = setup_test_path.rstrip("/")
        dir_with_suffix = check_qc.non_empty_existing_path(dir_without_suffix)
        assert dir_without_suffix[-1] != "/"
        assert dir_with_suffix[-1] == "/"


class TestReadYaml():
    def test_empty_yaml(self, setup_test_path):
        with pytest.raises(ValueError) as empty_error:
            check_qc.read_yaml(setup_test_path + "empty.yaml")
        assert "Could not load YAML." in str(empty_error.value)


class TestCheckAllowedOperators():
    def test_not_existing_operator(self):
        with pytest.raises(ValueError) as ops_error:
            check_qc.check_allowed_operators("=!")
        assert "Unsupported operator provided" in str(ops_error.value)
        assert "=!" in str(ops_error.value)


class TestCheckRequiredKeysMetrics():
    def test_required_keys_present(self):
        qc_metrics = [
            {"filename": "fake", "qc_col": "fake", "threshold": "fake", "operator": "fake",  "report_cols": "fake"},
        ]
        check_qc.check_required_keys_metrics(qc_metrics)
        assert True

    @pytest.mark.parametrize(
        "incomplete_qc_metrics",
        [
            [{"filename": "fakename"}],
            [
                {"filename": "fake", "qc_col": "fake", "threshold": "fake", "operator": "fake",  "report_cols": "fake"},
                # Missing report_cols
                {"filename": "fake", "qc_col": "fake", "threshold": "fake", "operator": "fake"},
            ],
            [
                # Missing report_cols
                {"filename": "fake", "qc_col": "fake", "threshold": "fake", "operator": "fake"},
                # Missing report_cols
                {"filename": "fake", "qc_col": "fake", "threshold": "fake", "operator": "fake"},
            ]
        ]
    )
    def test_missing_keys(self, incomplete_qc_metrics):
        with pytest.raises(KeyError) as required_error:
            check_qc.check_required_keys_metrics(incomplete_qc_metrics)
        error_val = str(required_error.value)
        assert "not in all metrics settings." in error_val
        assert error_val.split(" ")[2] in ["filename", "qc_col", "threshold", "operator", "report_cols"]


class TestSelectMetrics():
    @pytest.mark.parametrize("filename_or_regex,input_files,expected", [
        # Multi match
        ("test", ["test1.txt", "test2.txt"], ["test1.txt", "test2.txt"]),
        # Single match
        ("test", ["test1.txt", "fake2.txt"],  ["test1.txt"]),
        # Match with relative path
        ("test", ["./random/path/to/test1.txt"],  ["./random/path/to/test1.txt"]),
        # Match with absolute path
        ("test", ["/random/path/to/test1.txt"],  ["/random/path/to/test1.txt"]),
        # Match regex: kinship file suffix
        (
            ".*.kinship_check.out$",
            ["240101_A00295_0001_AHWCFKDSX7_CREv4_1.kinship_check.out"],
            ["240101_A00295_0001_AHWCFKDSX7_CREv4_1.kinship_check.out"]
        ),
        # Match on word truth and SNP
        (
            ".*truth.*SNP",
            [
                '2/U000000CFGIAB12878a_GIAB12878_nist2.19_truth_fix_header.vcf.gz_SNP_ALL.csv',
                '12/U000000CFGIAB12878b_U000000CFGIAB12878c_SNP_ALL.csv'
            ],
            ['2/U000000CFGIAB12878a_GIAB12878_nist2.19_truth_fix_header.vcf.gz_SNP_ALL.csv']
        ),
        # Match on word truth and SNP
        (
            ".*truth.*SNP",
            [
                'U000000CFGIAB12878a_GIAB12878_nist2.19_truth_fix_header.vcf.gz_SNP_ALL.csv',
                'U000000CFGIAB12878b_U000000CFGIAB12878c_SNP_ALL.csv'
            ],
            ['U000000CFGIAB12878a_GIAB12878_nist2.19_truth_fix_header.vcf.gz_SNP_ALL.csv']
        ),
        # Match if 'truth' is absent and contains 'SNP'
        # ?: Match expression but do not capture it
        # ?! Match if 'truth' is absent.
        (
            "(?:(?!truth).)*SNP.*$",
            [
                '2/U000000CFGIAB12878a_GIAB12878_nist2.19_truth_fix_header.vcf.gz_SNP_ALL.csv',
                '12/U000000CFGIAB12878b_U000000CFGIAB12878c_SNP_ALL.csv'
            ],
            ['12/U000000CFGIAB12878b_U000000CFGIAB12878c_SNP_ALL.csv']
        ),

    ])
    def test_select_metric(self, filename_or_regex, input_files, expected):
        metrics = check_qc.select_metrics(filename_or_regex, input_files)
        assert metrics == expected

    @pytest.mark.parametrize("filename_or_regex,input_files", [
        # No match
        ("test", ["fake1.txt", "fake2.txt"]),
        # No match, filename with specialchar assumed to be a regex
        ("specialchar_@", ["12/specialchar_@.txt"]),
    ])
    def test_no_match(self, filename_or_regex, input_files):
        with pytest.warns(UserWarning) as match_warning:
            return_val = check_qc.select_metrics("test", ["fake1.txt", "fake2.txt"])
        warn_msg = match_warning[0].message.args[0]
        assert "No input file provided with filename pattern" in warn_msg
        assert "test" in warn_msg
        assert not return_val


class TestGetColumnsToReport():
    @pytest.mark.parametrize("report_cols,metric_cols,qc_col,expected", [
        (["col1"], ["col1"], "col1", ["qc_title", "qc_value"]),
        # Additional report col
        (["col1", "col2"], ["col1", "col2"], "col1", ["qc_title", "qc_value", "col2"]),
        # Different order output
        (["col1", "col2"], ["col1", "col2"], "col2", ["qc_title", "qc_value", "col1"]),
        # Additional metric col
        (["col1"], ["col1", "col3"], "col1", ["qc_title", "qc_value"]),
        # Special @all option
        ("@all", ["col1", "col2"], "col1", ["qc_title", "qc_value", "col2"]),
    ])
    def test_get_columns_to_report(self, report_cols, metric_cols, qc_col, expected):
        qc_report_cols = check_qc.get_columns_to_report(report_cols, metric_cols, qc_col)
        assert qc_report_cols == unordered(expected)

    def test_wrong_type(self):
        with pytest.raises(TypeError) as type_error:
            check_qc.get_columns_to_report({"col1"}, ["col1"], "col1")
        assert "{'col1'}" in str(type_error.value)

    def test_non_existing_cols(self):
        with pytest.raises(ValueError) as not_exists_error:
            check_qc.get_columns_to_report(["col1", "col3"], ["col1"], "col1")
        assert "not exists" in str(not_exists_error.value)
        assert "col3" in str(not_exists_error.value)


class TestAddAndRenameColumns():
    def test_add_and_rename_columns(self):
        fake_qc_metric = DataFrame({"sample": ["sample1"], "fake_qc_col": ["0.01"]})
        qc_metric_out = check_qc.add_and_rename_columns(fake_qc_metric, "FAKE_title", "fake_qc_col", "fake_op", "fake_thres")
        # Assert expected column values
        assert qc_metric_out["qc_title"].values == "fake_title"
        assert qc_metric_out["qc_status"].values == "PASS"
        assert qc_metric_out["qc_check"].values == "fake_thres fake_op fake_qc_col"
        # Assert all expected columns exist
        assert not list(
            set(['sample', 'qc_value', 'qc_title', 'qc_status', 'qc_check', 'qc_msg']) - set(qc_metric_out.columns)
        )
        assert "fake_qc_col" not in qc_metric_out.columns


class TestGetFailedRows():
    @pytest.mark.parametrize("qc_op,qc_thres", [
        # Test match
        ("match", "fake_thres"),
        # Test string
        ("==", "FAIL"),
        # Test int
        ("==", 1),
        # Test float
        ("==", 0.1),
    ])
    def test_correct(self, qc_op, qc_thres):
        fake_qc_metric = DataFrame({"sample": ["sample1"], "fake_qc_col": [qc_thres]})
        rows = check_qc.get_failed_rows(fake_qc_metric, "fake_qc_col", qc_op, qc_thres)
        assert len(rows) == 1

    def test_not_support_type_threshold(self):
        with pytest.raises(TypeError) as type_error:
            check_qc.get_failed_rows(None, None, None, [])
        assert "[] type not supported" in str(type_error.value)


class TestAddFailedSamplesMetric():
    def test_add_failed_samples_single_sample_col(self):
        fake_qc_metric = DataFrame({
            "sample_col": ["sample1", "sample2"],
            "qc_check": [None, None],
            "qc_status": ["PASS", "PASS"],
            "qc_msg": [None, None],
            "qc_value": [0.1, 0.2],
        })
        failed_rows = fake_qc_metric.loc[fake_qc_metric["sample_col"] == "sample2"].index
        qc_metric, qc_metric_out = check_qc.add_failed_samples_metric(
            fake_qc_metric, failed_rows, fake_qc_metric.columns.to_list(), ["sample_col"])
        # Test rename column
        assert "sample" in qc_metric_out.columns.to_list()
        # Test removal failed sample
        assert "sample2" not in qc_metric["sample_col"].to_list()
        # Test added failed sample
        assert "sample2" in qc_metric_out["sample"].to_list()
        assert len(qc_metric) == 1 and len(qc_metric_out) == 1
        assert qc_metric_out["qc_status"].values == "FAIL"

    def test_add_failed_samples_multi_sample_col(self):
        fake_kinship_metric = DataFrame(
            data=list(combinations(["sample1", "sample2", "sample3", "sample4", "sample5"], 2)),
            columns=["sample_col1", "sample_col2"]
        )
        fake_kinship_metric = fake_kinship_metric.assign(qc_check="checks", qc_value="wrong")
        # Define sample1 vs sample2 and sample1 vs sample3 as failed
        failed_rows = fake_kinship_metric.iloc[0:2].index
        qc_metric, qc_metric_out = check_qc.add_failed_samples_metric(
            fake_kinship_metric, failed_rows, fake_kinship_metric.columns.to_list(), ["sample_col1", "sample_col2"])
        for failed_sample in ["sample1", "sample2", "sample3"]:
            # Test removal failed sample
            assert failed_sample not in list(qc_metric[["sample_col1", "sample_col2"]].values.ravel())
            # Test added failed sample
            assert failed_sample in qc_metric_out["sample"].to_list()
        assert qc_metric_out["qc_status"].values.all() == "FAIL"
        assert len(qc_metric) == 1
        twice_failed = qc_metric_out.loc[qc_metric_out["sample"] == "sample1"]
        # Assert join with ';' on column qc_value
        assert "wrong;wrong" == twice_failed["qc_value"].item()
        # Assert join with ';' on column qc_msg
        assert "sample1 sample2 checks wrong;sample1 sample3 checks wrong" == twice_failed["qc_msg"].item()
        for passed_sample in ["sample4", "sample5"]:
            assert passed_sample in list(qc_metric[["sample_col1", "sample_col2"]].values.ravel())

    def test_only_passed_rows(self):
        fake_qc_metric = DataFrame({"sample": ["sample1"], "fake_qc_col": ["0.1"]})
        failed_rows = DataFrame().index
        qc_metric, qc_metric_out = check_qc.add_failed_samples_metric(fake_qc_metric, failed_rows, None, None)
        assert qc_metric_out.empty
        assert fake_qc_metric.equals(qc_metric)


class TestAddPassedSamplesMetric():
    def test_add_passed_samples_multi_sample_col(self):
        fake_sample_qc = DataFrame([
                ("s1", "kinship", "FAIL", "s1 s2 kinship wrong;s1 s3 kinship wrong", "wrong;wrong"),
                ("s2", "kinship", "FAIL", "s1 s2 kinship wrong", "wrong"),
                ("s3", "kinship", "FAIL", "s1 s3 kinship wrong", "wrong"),
        ], columns=["sample", "qc_check", "qc_status", "qc_msg", "qc_value"])
        fake_qc_metric = DataFrame(
            [
                ("s4", "s5", "kinship", "ok", "PASS", "", "empty"),
                ("s4", "s6", "kinship", "ok", "PASS", "", "empty"),
                ("s5", "s6", "kinship", "ok", "PASS", "", "empty"),
            ],
            columns=["sample_col1", "sample_col2", "qc_check", "qc_value", "qc_status", "qc_msg", "new_col"]
        )
        qc_metric_out = check_qc.add_passed_samples_metric(
            fake_qc_metric, fake_sample_qc, ["sample_col1", "sample_col2"])
        # Test rename column
        assert "sample" in qc_metric_out.columns.to_list()
        # Test removal duplicates
        assert qc_metric_out["sample"].to_list().count("s4") == 1
        # Test additional columns ignored
        assert "new_col" not in qc_metric_out.columns.to_list()


class TestCreateAndWriteOutput():
    @pytest.mark.parametrize("exp_summary,qc_output", [
        # All qc checks passed
        ("PASS", DataFrame({"sample": ["s1"], "qc_status_cov": ["PASS"], "qc_status_kinship": ["PASS"]})),
        # Single qc check failed
        ("FAIL", DataFrame({"sample": ["s1"], "qc_status_cov": ["PASS"], "qc_status_kinship": ["FAIL"]})),
        # All qc check failed
        ("FAIL", DataFrame({"sample": ["s1"], "qc_status_cov": ["FAIL"], "qc_status_kinship": ["FAIL"]})),
        # Not restricted to qc_status_<check> column name
        ("PASS", DataFrame({"sample": ["s1"], "random_col1": ["PASS"], "random_col2": ["PASS"]})),
    ])
    def test_create_and_write_output(self, setup_test_path, exp_summary, qc_output):
        prefix = "test_output"
        check_qc.create_and_write_output(qc_output, setup_test_path, prefix)
        expected_output = Path(f"{setup_test_path}/{prefix}_summary.csv")
        assert expected_output.exists()
        out = read_csv(expected_output)
        assert "qc_summary" in out.columns.to_list()
        assert out["qc_summary"].values == exp_summary


class TestGetOutputMetrics():
    @pytest.mark.parametrize("data_in,nr_rows", [
        # Single sample
        (["sample1_fake_check.txt"], 1),
        # Multiple single samples
        (["sample1_fake_check.txt", "sample2_fake_check.txt"], 2),
        # Single multi samples
        (["240101_fake_check.txt"], 2),
        # Multiple multi samples
        (["240101_fake_check.txt", "240102_fake_check.txt"], 4),
        # Multi and single sample
        (["sample1_fake_check.txt", "240101_fake_check.txt"], 3),
    ])
    def test_input_ok(self, data_in, nr_rows, dataset, datadir):
        datadir_files = [f"{datadir}/{filename}" for filename in data_in]
        df_output = check_qc.read_and_judge_metrics(dataset["settings_single_metric"]["metrics"][0], datadir_files)
        assert not df_output.empty
        observed_cols = df_output.columns.to_list()
        # Shape results in tuple with no. rows and no. cols
        assert df_output.shape[0] == nr_rows
        assert len(observed_cols) == 5
        assert observed_cols == ['sample', 'qc_check_fc', 'qc_status_fc', 'qc_msg_fc', 'qc_value_fc']

    @pytest.mark.parametrize("data_in,nr_rows,exp_warn_msg", [
        # Single sample duplicate
        (["sample1_fake_check.txt"]*2, 1, "Sample IDs occur multiple times in input:"),
        # Single multi samples duplicate
        (["240101_fake_check.txt"]*2, 2, "Sample IDs occur multiple times in input:"),
        # Multiple multi samples, duplicate samples
        (["240101_fake_check.txt", "240101_v2_fake_check.txt"], 4, "Different qc values for duplicated sample IDs in input:"),
    ])
    def test_input_warn(self, data_in, nr_rows, exp_warn_msg, dataset, datadir):
        datadir_files = [f"{datadir}/{filename}" for filename in data_in]
        with pytest.warns(UserWarning) as match_warning:
            df_output = check_qc.read_and_judge_metrics(dataset["settings_single_metric"]["metrics"][0], datadir_files)
        warn_msg = match_warning[0].message.args[0]
        assert exp_warn_msg in warn_msg
        assert not df_output.empty
        observed_cols = df_output.columns.to_list()
        # Shape: tuple with no. rows and no. cols
        assert df_output.shape[0] == nr_rows
        assert len(observed_cols) == 5
        assert observed_cols == ['sample', 'qc_check_fc', 'qc_status_fc', 'qc_msg_fc', 'qc_value_fc']


class TestCheckQc():
    @pytest.mark.parametrize("settings,data_in,exp_shape", [
        # Single metric, single sample input
        ("settings_single_metric", ["sample1_fake_check.txt"], (1, 5)),
        # Two metrics, single sample input
        ("settings_two_metrics", ["sample1_fake_check.txt"], (1, 9)),
        # Single metric, multiple samples input
        ("settings_single_metric", ["240101_fake_check.txt"], (2, 5)),
        ("settings_single_metric", ["240101_fake_check.txt", "240102_fake_check.txt"], (4, 5)),
        # Two metrics, multiple sample input
        ("settings_two_metrics", ["240101_fake_check.txt", "240102_fake_check.txt"], (4, 9)),
        # Two metricS, multi and single sample input
        ("settings_two_metrics", ["sample1_fake_check.txt", "240101_fake_check.txt"], (3, 9)),
    ])
    def test_ok(self, settings, data_in, exp_shape, datadir, dataset, mocker, ):
        datadir_files = [f"{datadir}/{filename}" for filename in data_in]
        mocker.patch("check_qc.read_yaml", return_value=dataset[settings])
        mock_write_output = mocker.patch("check_qc.create_and_write_output")
        check_qc.check_qc(input_files=datadir_files, settings="", output_path="", output_prefix="")
        mock_write_output.assert_called_once()
        # Shape: tuple with no. rows and no. cols
        assert mock_write_output.call_args[0][0].shape == exp_shape
        mock_write_output.reset_mock()

    def test_no_match_input_error(self, mocker, mock_settings):
        mock_select_metrics = mocker.patch("check_qc.select_metrics", return_value=None)
        mock_get_output = mocker.patch("check_qc.read_and_judge_metrics")
        with pytest.raises(ValueError) as no_match_error:
            check_qc.check_qc(input_files=[], settings="", output_path="", output_prefix="")
        mock_select_metrics.assert_called_once()
        assert not mock_get_output.called
        assert "No input files found to match any qc metric pattern." == str(no_match_error.value)
        mock_settings.reset_mock()

    def test_duplicate_samples_error(self, datadir, mocker, mock_settings):
        mock_pandas_merge = mocker.patch("pandas.merge")
        with pytest.raises(ValueError) as duplicate_error:
            check_qc.check_qc(input_files=[f"{datadir}/240101_fake_check.txt", f"{datadir}/240101_v2_fake_check.txt"],
                              settings="", output_path="", output_prefix="")
        assert "Duplicated samples with different values found in files matching" in str(duplicate_error.value)
        assert "fake_check.txt" in str(duplicate_error.value)
        assert not mock_pandas_merge.called
        mock_settings.reset_mock()
