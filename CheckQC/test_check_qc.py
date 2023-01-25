#!/usr/bin/env python
# Import statements, alphabetic order of main package.
from itertools import combinations
import pytest
from pytest_unordered import unordered

# Third party libraries alphabetic order of main package.
from pandas import DataFrame

# Custom libraries alphabetic order of main package.
import check_qc


@pytest.fixture(scope="module", autouse=True)
def setup_test_path(tmp_path_factory):
    test_tmp_path = str(tmp_path_factory.mktemp("test")) + "/"
    # create empty files
    open(str(test_tmp_path) + "/empty.txt", "a").close()
    open(str(test_tmp_path) + "/empty.yaml", "a").close()
    return test_tmp_path


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
        assert "File is empty." in str(empty_error.value)

    def test_append_suffix(self, setup_test_path):
        dir_without_suffix = setup_test_path.rstrip("/")
        dir_with_suffix = check_qc.non_empty_existing_path(dir_without_suffix)
        assert dir_without_suffix[-1] != "/"
        assert dir_with_suffix[-1] == "/"


class TestReadYaml():
    def test_empty_yaml(self, setup_test_path):
        with pytest.raises(ValueError) as empty_error:
            check_qc.read_yaml(setup_test_path + "empty.yaml")
        assert "File is empty." in str(empty_error.value)


class TestCheckAllowedOperators():
    def test_not_existing_operator(self):
        with pytest.raises(ValueError) as ops_error:
            check_qc.check_allowed_operators("=!")
        assert "Unsupported operator provided" in str(ops_error.value)
        assert "=!" in str(ops_error.value)


class TestCheckRequiredKeysMetrics():
    def test_required_keys_present(self):
        qc_settings = {"metrics": [
            {"filename": "fake", "qc_col": "fake", "threshold": "fake", "operator": "fake",  "report_cols": "fake"},
        ]}
        check_qc.check_required_keys_metrics(qc_settings)
        assert True

    @pytest.mark.parametrize(
        "incomplete_qc_settings",
        [
            {"metrics": [{"filename": "fakename"}]},
            {"metrics": [
                {"filename": "fake", "qc_col": "fake", "threshold": "fake", "operator": "fake",  "report_cols": "fake"},
                {"filename": "fake", "qc_col": "fake", "threshold": "fake", "operator": "fake"},  # missing report_cols
            ]},
            {"metrics": [
                {"filename": "fake", "qc_col": "fake", "threshold": "fake", "operator": "fake"},  # missing report_cols
                {"filename": "fake", "qc_col": "fake", "threshold": "fake", "operator": "fake"},  # missing report_cols
            ]}
        ]
    )
    def test_missing_keys(self, incomplete_qc_settings):
        with pytest.raises(KeyError) as required_error:
            check_qc.check_required_keys_metrics(incomplete_qc_settings)
        error_val = str(required_error.value)
        assert "not in all metrics settings." in error_val
        assert error_val.split(" ")[2] in ["filename", "qc_col", "threshold", "operator", "report_cols"]
    

class TestSelectMetrics():
    @pytest.mark.parametrize("input_files,expected", [
        (["test1.txt", "test2.txt"], ["test1.txt", "test2.txt"]),  # multi match
        (["test1.txt", "fake2.txt"],  ["test1.txt"]),  # single match
    ])
    def test_select_metric(self, input_files, expected):
        metrics = check_qc.select_metrics("test", input_files)
        assert metrics == expected

    def test_no_match(self):
        with pytest.raises(ValueError) as match_error:
            check_qc.select_metrics("test", ["fake1.txt", "fake2.txt"])
        error_val = str(match_error.value)
        assert "No input file provided with filename pattern" in error_val
        assert "test" in error_val


class TestGetColumnsToReport():
    @pytest.mark.parametrize("report_cols,metric_cols,qc_col,expected", [
        (["col1"], ["col1"], "col1", ["qc_title", "qc_value"]),
        (["col1", "col2"], ["col1", "col2"], "col1", ["qc_title", "qc_value", "col2"]),  # additional report col
        (["col1", "col2"], ["col1", "col2"], "col2", ["qc_title", "qc_value", "col1"]),  # additional report col, different order output
        (["col1"], ["col1", "col3"], "col1", ["qc_title", "qc_value"]),  # additional metric col

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
        # assert expected column values
        assert qc_metric_out["qc_title"].values == "fake_title"
        assert qc_metric_out["qc_status"].values == "PASS"
        assert qc_metric_out["qc_check"].values == "fake_thres fake_op fake_qc_col"
        # assert all expected columns exist
        assert not list(set(['sample', 'qc_value', 'qc_title', 'qc_status', 'qc_check', 'qc_msg']) - set(qc_metric_out.columns))
        assert "fake_qc_col" not in qc_metric_out.columns


class TestGetFailedRows():
    @pytest.mark.parametrize("qc_op,qc_thres", [
        ("match", "fake_thres"),  # test match
        ("==", "FAIL"),  # test string
        ("==", "1"),  # test int
        ("==", ""),  # test float
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
        assert "sample" in qc_metric_out.columns.to_list()  # test rename column
        assert "sample2" not in qc_metric["sample_col"].to_list()  # test removal failed sample
        assert "sample2" in qc_metric_out["sample"].to_list()  # test added failed sample
        assert len(qc_metric) == 1 and len(qc_metric_out) == 1
        assert qc_metric_out["qc_status"].values == "FAIL"

    def test_add_failed_samples_multi_sample_col(self):
        fake_kinship_metric = DataFrame(
            data=list(combinations(["sample1", "sample2", "sample3", "sample4", "sample5"], 2)), 
            columns=["sample_col1", "sample_col2"]
        )
        fake_kinship_metric = fake_kinship_metric.assign(qc_check="val op thres", qc_value="wrong")
        failed_rows = fake_kinship_metric.iloc[0:2].index  # define sample1 vs sample2 and sample1 vs sample3 as failed
        qc_metric, qc_metric_out = check_qc.add_failed_samples_metric(
            fake_kinship_metric, failed_rows, fake_kinship_metric.columns.to_list(), ["sample_col1", "sample_col2"])
        for failed_sample in ["sample1", "sample2", "sample3"]:
            assert failed_sample not in list(qc_metric[["sample_col1", "sample_col2"]].values.ravel())  # test removal failed sample
            assert failed_sample in qc_metric_out["sample"].to_list()  # test added failed sample
        assert qc_metric_out["qc_status"].values.all() == "FAIL"
        assert len(qc_metric) == 1
        for passed_sample in ["sample4", "sample5"]:
            assert passed_sample in list(qc_metric[["sample_col1", "sample_col2"]].values.ravel())

    def test_only_passed_rows(self):
        fake_qc_metric = DataFrame({"sample": ["sample1"], "fake_qc_col": ["0.1"]})
        failed_rows = DataFrame().index
        qc_metric, qc_metric_out = check_qc.add_failed_samples_metric(fake_qc_metric, failed_rows, None, None)
        assert qc_metric_out.empty
        assert fake_qc_metric.equals(qc_metric)

