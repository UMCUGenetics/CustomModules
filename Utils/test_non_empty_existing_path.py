#!/usr/bin/env python
# Import statements, alphabetic order of main package.
import pytest

# Third party libraries alphabetic order of main package.

# Custom libraries alphabetic order of main package.
from CustomModules.Utils.non_empty_existing_path import non_empty_existing_path


@pytest.fixture(scope="module", autouse=True)
def setup_test_path(tmp_path_factory):
    test_tmp_path = str(tmp_path_factory.mktemp("test")) + "/"
    # create empty files
    open(str(test_tmp_path) + "/empty.txt", "a").close()
    return test_tmp_path


class TestNonEmptyExistingPath():
    def test_existing_dir(self, setup_test_path):
        file_or_dir = non_empty_existing_path(setup_test_path)
        assert file_or_dir

    def test_not_file_not_dir(self):
        fake_string = "fake_string"
        with pytest.raises(FileNotFoundError) as file_dir_error:
            non_empty_existing_path(fake_string)
        assert fake_string in str(file_dir_error.value)

    def test_empty_file(self, setup_test_path):
        with pytest.raises(OSError) as empty_error:
            non_empty_existing_path(setup_test_path + "empty.txt")
        assert f"File {setup_test_path}empty.txt is empty." in str(empty_error.value)

    def test_append_suffix(self, setup_test_path):
        dir_without_suffix = setup_test_path.rstrip("/")
        dir_with_suffix = non_empty_existing_path(dir_without_suffix)
        assert dir_without_suffix[-1] != "/"
        assert dir_with_suffix[-1] == "/"
