#!/usr/bin/env python
# Import statements, alphabetic order of main package.
import pytest

# Custom libraries alphabetic order of main package.
from CustomModules.GenderCheck import calculate_gender


class TestCompareGender():
    @pytest.mark.parametrize("sample_id,analysis_id,test_gender,true_gender,expected", [
        # test_gender and true_gender identical, should be PASS
        ("test_sample", "test_analyse", "male", "male", "test_sample\ttest_analyse\tmale\tmale\tPASS\n"),
        # test_gender and true_gender not identical , should be FAIL
        ("test_sample", "test_analyse", "male", "female", "test_sample\ttest_analyse\tmale\tfemale\tFAIL\n"),
        # true_gender unknown, should be PASS
        ("test_sample", "test_analyse", "male", "unknown", "test_sample\ttest_analyse\tmale\tunknown\tPASS\n"),
        # true_gender not_detected, should be FAIL
        ("test_sample", "test_analyse", "male", "not_detected", "test_sample\ttest_analyse\tmale\tnot_detected\tFAIL\n"),
    ])
    def test_compare_gender(self, sample_id, analysis_id, test_gender, true_gender, expected):
        assert expected == calculate_gender.compare_gender(sample_id, analysis_id, test_gender, true_gender)


def test_write_qc_file(tmp_path):
    path = tmp_path / "qc_folder"
    path.mkdir()
    qc_file = path / "test_sample_test_analyse_gendercheck.txt"
    calculate_gender.write_qc_file("test_sample", "test_analyse", "test_sample\ttest_analyse\tmale\tmale\tPASS\n", path)
    message = "sample_id\tanalysis_id\ttest_gender\ttrue_gender\tstatus\ntest_sample\ttest_analyse\tmale\tmale\tPASS\n"
    assert message in qc_file.read_text()
