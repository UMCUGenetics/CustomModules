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


@pytest.mark.parametrize("input_gender,exp_output", [
    # Given gender should return translated gender
    ("Man", "male"),
    ("Vrouw", "female"),
    ("Onbekend", "unknown"),
    ("unknown", "not_detected"),
    # Given gender should be returned
    ("MAN", "MAN")
])
def test_translate_gender(input_gender, exp_output):
    out = calculate_gender.translate_gender(input_gender)
    assert out == exp_output


class TestValidateGender():
    @pytest.mark.parametrize("input_gender", [("male"), ("female"), ("unknown"), ("not_detected")])
    def test_allowed_genders(self, input_gender):
        try:
            calculate_gender.validate_gender(input_gender)
        except ValueError:
            assert False
        assert True

    @pytest.mark.parametrize("input_gender", [("Man"), ("Vrouw"), ("Onbekend"), ("fakeGender")])
    def test_not_allowed_genders(self, input_gender):
        with pytest.raises(ValueError) as value_error:
            calculate_gender.validate_gender(input_gender)
        assert f"Provided gender {input_gender} is not allowed. Should be one of ['male', 'female', 'unknown', 'not_detected']." == str(value_error.value)


class TestGetGenderFromBam():
    @pytest.mark.parametrize("bam,mapping_qual,locus_y,ratio_y,expected", [
        ("./test_bam.bam", 20, "Y:2649520-59034050", 0.02, "male"),  # output male below
        ("./test_bam.bam", 20, "Y:2649520-59034050", 0.22, "female"),  # output female
    ])
    def test_get_gender_from_bam(self, bam, mapping_qual, locus_y, ratio_y, expected):
        assert expected == calculate_gender.get_gender_from_bam(bam, mapping_qual, locus_y, ratio_y)


class TestCompareGender():
    @pytest.mark.parametrize("test_gender,true_gender,expected_qc,expected_msg", [
        # test_gender and true_gender identical, should be PASS
        ("male", "male", "PASS", ""),
        # test_gender and true_gender not identical , should be FAIL
        ("male", "female", "FAIL", "True gender female does not equal estimated gender male."),
        # true_gender unknown, should be PASS
        ("male", "unknown", "PASS", ""),
        # true_gender not_detected, should be FAIL
        (
            "male", "not_detected", "FAIL",
            "Gender has value 'not_detected' in LIMS. Observed gender 'male' could not be verified."
        ),
    ])
    def test_compare_gender(self, test_gender, true_gender, expected_qc, expected_msg):
        status, message = calculate_gender.compare_gender(test_gender, true_gender)
        assert status == expected_qc
        assert message == expected_msg


def test_write_qc_file(tmp_path):
    path = tmp_path / "qc_folder"
    path.mkdir()
    qc_file = path / "test_sample_test_analyse_gendercheck.txt"
    calculate_gender.write_qc_file("test_sample", "test_analyse", "male", "male", "PASS", "", path)
    lines = qc_file.read_text().splitlines()
    assert len(lines) == 2
    assert lines[0] == "sample_id\tanalysis_id\ttest_gender\ttrue_gender\tstatus\tmessage"
    assert lines[1] == "test_sample\ttest_analyse\tmale\tmale\tPASS\t"
