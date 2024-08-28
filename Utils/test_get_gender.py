#!/usr/bin/env python
# Import statements, alphabetic order of main package.
import pytest

# Custom libraries alphabetic order of main package.
import get_gender

class TestIsValidRead():
    class ValidReadObject:
        def __init__(self, qual, start, end):
            self.mapping_quality = qual
            self.reference_start = start
            self.reference_end = end

    @pytest.mark.parametrize("read,mapping_qual,expected", [
        (ValidReadObject(19, True, True), 20, False),  # mapping quality is below the threshold
        (ValidReadObject(20, True, True), 20, True),  # mapping quality is equal to the threshold
        (ValidReadObject(20, True, True), 19, True),  # mapping quality is higher than the threshold
        (ValidReadObject(20, False, True), 20, False),  # reference_end is false
        (ValidReadObject(20, True, False), 20, False),  # reference_start is false
    ])
    def test_is_valid_read(self, read, mapping_qual, expected):
        assert expected == get_gender.is_valid_read(read, mapping_qual)


class TestGetRatioFromBam:
    @pytest.mark.parametrize("bam,mapping_qual,locus_x,expected_outcome", [
        ("test_bam.bam", 20, "X:2699520-154931044", 4.75),
        ("test_bam.bam", 20, "Y:2649520-59034050", 0.22)
    ])
    def test_get_gender_from_bam(self, bam, mapping_qual, locus_x, expected_outcome, shared_datadir):
        assert expected_outcome == round(
            get_gender.get_ratio_from_bam(f"{shared_datadir}/{bam}", mapping_qual, locus_x),
            ndigits=2
        )


class TestGetGenderOnXRatio:
    @pytest.mark.parametrize("ratio,ratio_x_threshold_male,ratio_x_threshold_female,expected_outcome", [
        (4.75, 3.5, 4.5, ("F", False)),
        (4.75, 5.5, 7.5, ("M", False)),
        (4.75, 4.5, 6.5, ("F", True)),
    ])
    def test_get_gender_from_bam(self, ratio, ratio_x_threshold_male, ratio_x_threshold_female, expected_outcome):
        assert expected_outcome == get_gender.get_gender_on_x_ratio(ratio, ratio_x_threshold_male, ratio_x_threshold_female)


class TestGetGenderOnYRatio:
    @pytest.mark.parametrize("ratio,ratio_y_threshold,expected_outcome", [
        (0.22, 0.02, "male"),
        (0.22, 0.22, "female"),
    ])
    def test_get_gender_from_bam(self, ratio, ratio_y_threshold, expected_outcome):
        assert expected_outcome == get_gender.get_gender_on_y_ratio(ratio, ratio_y_threshold)
