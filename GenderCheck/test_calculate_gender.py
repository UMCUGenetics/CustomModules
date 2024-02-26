import calculate_gender

import pytest


class TestIsValidRead():


    class MyObject:
        def __init__(self, qual, start, end):
            self.mapping_quality = qual
            self.reference_start = start
            self.reference_end = end


    @pytest.mark.parametrize("read,mapping_qual,expected", [
        (MyObject(19, True, True), 20, False),  # mapping quality is below the threshold
        (MyObject(20, True, True), 20, True),  # mapping quality is equal to the threshold
        (MyObject(20, True, True), 19, True),  # mapping quality is higher than the threshold
        (MyObject(20, False, True), 20, False),  # reference_end is false
        (MyObject(20, True, False), 20, False),  # reference_start is false
    ])
    def test_is_valid_read(self, read, mapping_qual, expected):
        assert expected == calculate_gender.is_valid_read(read, mapping_qual)


class TestGetGenderFromBam():
    @pytest.mark.parametrize("bam,mapping_qual,locus_y,ratio_y_female,ratio_y_male,expected", [
        ("./test_bam.bam", 20, "Y:2649520-59034050", 0.02, 0.12, "male"),  # output male
        ("./test_bam.bam", 20, "Y:2649520-59034050", 0.22, 0.32, "female"),  # output female
        ("./test_bam.bam", 20, "Y:2649520-59034050", 0.02, 0.32, "unknown"),  # output unknown
    ])
    def test_get_gender_from_bam(self, bam, mapping_qual, locus_y, ratio_y_female, ratio_y_male, expected):
        assert expected == calculate_gender.get_gender_from_bam(bam, mapping_qual, locus_y, ratio_y_female, ratio_y_male)


class TestCompareGender():
    @pytest.mark.parametrize("sample_id,analysis_id,test_gender,true_gender,expected", [
        ("test_sample", "test_analyse", "male", "male", "test_sample\ttest_analyse\tmale\tmale\tPASS\n"),  # test_gender and true_gender identical, PASS
        ("test_sample", "test_analyse", "male", "female", "test_sample\ttest_analyse\tmale\tfemale\tFAIL\n"),  # test_gender and true_gender not identical , FAIL
        ("test_sample", "test_analyse", "male", "unknown", "test_sample\ttest_analyse\tmale\tunknown\tPASS\n"),  # true_gender unknown, PASS
        ("test_sample", "test_analyse", "male", "not_detected", "test_sample\ttest_analyse\tmale\tnot_detected\tFAIL\n"),  # true_gender not_detected, FAIL
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
