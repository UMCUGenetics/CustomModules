#!/usr/bin/env python
# Import statements, alphabetic order of main package.
import pytest

# Custom libraries alphabetic order of main package.
import get_gender_from_bam_chrx

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
