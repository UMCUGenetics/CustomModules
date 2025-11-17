#!/usr/bin/env python
# Import statements, alphabetic order of main package.
# Third party libraries alphabetic order of main package.
import pytest

# Custom libraries alphabetic order of main package.
import calculate_gender


class TestIsValidRead():

    class MyObject:
        def __init__(self, qual, start, end):
            self.mapping_quality = qual
            self.reference_start = start
            self.reference_end = end

    @pytest.mark.parametrize("read,min_mapping_qual,expected", [
        (MyObject(19, True, True), 20, False),  # Mapping quality is below the threshold
        (MyObject(20, True, True), 20, True),  # Mapping quality is equal to the threshold
        (MyObject(20, True, True), 19, True),  # Mapping quality is higher than the threshold
        (MyObject(20, False, True), 20, False),  # Reference_end is false
        (MyObject(20, True, False), 20, False),  # Reference_start is false
    ])
    def test_is_valid_read(self, read, min_mapping_qual, expected):
        assert expected == calculate_gender.is_valid_read(read, min_mapping_qual)


@pytest.mark.parametrize("input_gender,exp_output", [
    # Given gender should return translated gender
    ("Man", "male"),
    ("Vrouw", "female"),
    ("Onbekend", "unknown"),
    ("unknown", "gender_data_not_found"),
    ("no_sample_found", "sample_not_found"),
    ("multiple_values_for_udf", "multiple_values_for_udf"),
    # Given gender should be returned
    ("MAN", "MAN")
])
def test_translate_gender(input_gender, exp_output):
    out = calculate_gender.translate_gender(input_gender)
    assert out == exp_output


class TestValidateGender():
    @pytest.mark.parametrize(
        "input_gender",
        [("male"), ("female"), ("unknown"), ("gender_data_not_found"), ("sample_not_found"), ("multiple_values_for_udf")]
    )
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
        assert (
            f"Provided gender {input_gender} is not allowed. "
            "Should be one of ['male', 'female', 'unknown', 'gender_data_not_found', 'sample_not_found', "
            "'multiple_values_for_udf']."
        ) == str(value_error.value)


class TestGetGenderFromBam():
    @pytest.mark.parametrize("bam,min_mapping_qual,locus_y,ratio_y,expected", [
        ("test_bam.bam", 20, "Y:2649520-59034050", 0.02, "male"),  # Output male below
        ("test_bam.bam", 20, "Y:2649520-59034050", 0.22, "female"),  # Output female
    ])
    def test_get_gender_from_bam(self, bam, min_mapping_qual, locus_y, ratio_y, expected, datadir):
        assert expected == calculate_gender.get_gender_from_bam(f"{datadir}/{bam}", min_mapping_qual, locus_y, ratio_y)


class TestCompareGender():
    @pytest.mark.parametrize("measured_gender,stated_gender,expected_qc,expected_msg", [
        # measured_gender and stated_gender identical, should be PASS
        ("male", "male", "PASS", ""),
        # measured_gender and stated_gender not identical (upper vs lowercase), should fail
        ("MALE", "male", "FAIL", "Stated gender male does not equal observed gender MALE."),
        # measured_gender and stated_gender not identical (upper vs lowercase), should fail
        ("FEMALE", "female", "FAIL", "Stated gender female does not equal observed gender FEMALE."),
        # measured_gender and stated_gender not identical , should be FAIL
        ("male", "female", "FAIL", "Stated gender female does not equal observed gender male."),
        # measured_gender and stated_gender not identical , should be FAIL
        ("fakegender", "female", "FAIL", "Stated gender female does not equal observed gender fakegender."),
        # stated_gender unknown, should be PASS
        ("male", "unknown", "PASS", ""),
        # stated_gender gender_data_not_found, should be FAIL
        (
            "male", "gender_data_not_found", "FAIL",
            "Gender has value 'gender_data_not_found' in LIMS. Observed gender 'male' could not be verified."
        ),
        # stated_gender sample_not_found, should be FAIL
        (
            "male", "sample_not_found", "FAIL",
            "Gender has value 'sample_not_found' in LIMS. Observed gender 'male' could not be verified."
        ),
        # stated_gender multiple_values_for_udf, should be FAIL
        (
            "male", "multiple_values_for_udf", "FAIL",
            "Gender has value 'multiple_values_for_udf' in LIMS. Observed gender 'male' could not be verified."
        )
    ])
    def test_compare_and_evaluate_gender(self, measured_gender, stated_gender, expected_qc, expected_msg):
        status, message = calculate_gender.compare_and_evaluate_gender(measured_gender, stated_gender)
        assert status == expected_qc
        assert message == expected_msg


def test_write_qc_file(tmp_path):
    path = tmp_path / "qc_folder"
    path.mkdir()
    qc_file = path / "test_sample_test_analyse_gendercheck.txt"
    calculate_gender.write_qc_file("test_sample", "test_analyse", "male", "male", "PASS", "", path)
    lines = qc_file.read_text().splitlines()
    assert len(lines) == 2
    assert lines[0] == "sample_id\tanalysis_id\tmeasured_gender\tstated_gender\tstatus\tmessage"
    assert lines[1] == "test_sample\ttest_analyse\tmale\tmale\tPASS\t"
