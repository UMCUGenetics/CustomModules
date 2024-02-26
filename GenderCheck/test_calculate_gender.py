import calculate_gender


class MyObject:
    def __init__(self, qual, start, end):
        self.mapping_quality = qual
        self.reference_start = start
        self.reference_end = end


def test_is_valid_read():
    assert calculate_gender.is_valid_read(MyObject(19, True, True), 20) is False
    assert calculate_gender.is_valid_read(MyObject(20, True, True), 20) is True
    assert calculate_gender.is_valid_read(MyObject(20, True, True), 19) is True
    assert calculate_gender.is_valid_read(MyObject(20, False, True), 20) is False
    assert calculate_gender.is_valid_read(MyObject(20, True, False), 20) is False


def test_get_gender_from_bam():
    assert calculate_gender.get_gender_from_bam("./test_bam.bam", 20, "Y:2649520-59034050", 0.02, 0.12) == "male"
    assert calculate_gender.get_gender_from_bam("./test_bam.bam", 20, "Y:2649520-59034050", 0.22, 0.32) == "female"
    assert calculate_gender.get_gender_from_bam("./test_bam.bam", 20, "Y:2649520-59034050", 0.02, 0.32) == "unknown"


def test_compare_gender():
    assert (
        calculate_gender.compare_gender("test_sample", "test_analyse", "male", "male") ==
        "test_sample\ttest_analyse\tmale\tmale\tPASS\n"
    )
    assert (
        calculate_gender.compare_gender("test_sample", "test_analyse", "male", "female") ==
        "test_sample\ttest_analyse\tmale\tfemale\tFAIL\n"
    )
    assert (
        calculate_gender.compare_gender("test_sample", "test_analyse", "male", "unknown") ==
        "test_sample\ttest_analyse\tmale\tunknown\tPASS\n"
    )
    assert (
        calculate_gender.compare_gender("test_sample", "test_analyse", "male", "not_detected") ==
        "test_sample\ttest_analyse\tmale\tnot_detected\tFAIL\n"
    )


def test_write_qc_file(tmp_path):
    path = tmp_path / "qc_folder"
    path.mkdir()
    qc_file = path / "test_sample_test_analyse_gendercheck.txt"
    calculate_gender.write_qc_file("test_sample", "test_analyse", "test_sample\ttest_analyse\tmale\tmale\tPASS\n", path)
    message = "sample_id\tanalysis_id\ttest_gender\ttrue_gender\tstatus\ntest_sample\ttest_analyse\tmale\tmale\tPASS\n"
    assert message in qc_file.read_text()
