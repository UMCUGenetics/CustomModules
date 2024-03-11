import get_gender_from_bam_chrx
import pytest


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
        assert expected == get_gender_from_bam_chrx.is_valid_read(read, mapping_qual)


class TestGetGenderFromBam:
    @pytest.mark.parametrize("bam,mapping_qual,locus_x,ratio_x_threshold_male,ratio_x_threshold_female", [
        ("./test_bam.bam", 20, "X:2699520-154931044", 5.5, 6.5),
        ("./test_bam.bam", 20, "X:2699520-154931044", 8, 6.5),
    ])
    def test_get_gender_from_bam(self, bam, mapping_qual, locus_x, ratio_x_threshold_male, ratio_x_threshold_female):
        assert ("F", False) == get_gender_from_bam_chrx.get_gender_from_bam_chrx(
            bam, mapping_qual, locus_x, ratio_x_threshold_male, ratio_x_threshold_female)


class TestWriteGenderDataToFile:
    @pytest.mark.parametrize("sample_id,gender_data", [
        ("ThisIsASampleID1", ("M", False)),
        ("ThisIsASampleID2", ("M", True)),
        ("ThisIsASampleID3", ("F", False)),
        ("ThisIsASampleID4", ("F", True))
    ])
    def test_write_gender_data_to_file(self, tmpdir, sample_id, gender_data):
        print(tmpdir)

        get_gender_from_bam_chrx.write_genderdata_to_file(sample_id, gender_data, tmpdir)
        file = tmpdir.join(f"/gender_data_{sample_id}.tsv")
        assert f"{sample_id}\t{gender_data[0]}\t{gender_data[1]}\n" in file.read()
