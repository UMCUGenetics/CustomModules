import calculate_gender


# Test is_valid_read

class MyObject:
    def __init__(self, qual, start, end):
        self.mapping_quality = qual
        self.reference_start = start
        self.reference_end = end


def test_is_valid_read_1a():
    assert calculate_gender.is_valid_read(MyObject(19, True, True), 20) is False


def test_is_valid_read_1b():
    assert calculate_gender.is_valid_read(MyObject(20, True, True), 20) is True


def test_is_valid_read_1c():
    assert calculate_gender.is_valid_read(MyObject(20, True, True), 19) is True


def test_is_valid_read_1d():
    assert calculate_gender.is_valid_read(MyObject(20, True, True), 21) is False


def test_is_valid_read_1e():
    assert calculate_gender.is_valid_read(MyObject(21, True, True), 20) is True


def test_is_valid_read_2():
    assert calculate_gender.is_valid_read(MyObject(20, False, True), 20) is False


def test_is_valid_read_3():
    assert calculate_gender.is_valid_read(MyObject(20, True, False), 20) is False


# Test get_gender_from_bam


def test_get_gender_from_bam_1():
    assert calculate_gender.get_gender_from_bam("./test_bam.bam", 20, "Y:2649520-59034050", 0.02, 0.12) == "male"


def test_get_gender_from_bam_2():
    assert calculate_gender.get_gender_from_bam("./test_bam.bam", 20, "Y:2649520-59034050", 0.22, 0.32) == "female"


def test_get_gender_from_bam_3():
    assert calculate_gender.get_gender_from_bam("./test_bam.bam", 20, "Y:2649520-59034050", 0.02, 0.32) == "unknown"


# Test compare_gender


def test_compare_gender_1():
    assert (
        calculate_gender.compare_gender("test_sample", "test_analyse", "male", "male") ==
        "test_sample\ttest_analyse\tmale\tmale\tPASS\n"
    )


def test_compare_gender_2():
    assert (
        calculate_gender.compare_gender("test_sample", "test_analyse", "male", "female") ==
        "test_sample\ttest_analyse\tmale\tfemale\tFAIL\n"
    )


def test_compare_gender_3():
    assert (
        calculate_gender.compare_gender("test_sample", "test_analyse", "male", "unknown") ==
        "test_sample\ttest_analyse\tmale\tunknown\tPASS\n"
    )


# Test write_qc_file


def test_write_qc_file(tmp_path):
    path = tmp_path / "qc_folder"
    path.mkdir()
    qc_file = path / "test_sample_test_analyse_gendercheck.txt"
    calculate_gender.write_qc_file("test_sample", "test_analyse", "test_sample\ttest_analyse\tmale\tmale\tPASS\n", path)
    message = "sample_id\tanalysis_id\ttest_gender\ttrue_gender\tstatus\ntest_sample\ttest_analyse\tmale\tmale\tPASS\n"
    assert message in qc_file.read_text()
