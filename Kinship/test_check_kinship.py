#!/usr/bin/env python
# Import statements, alphabetic order of main package.
from pathlib import Path


# Third party libraries alphabetic order of main package.
from pandas import DataFrame
import pytest
from pytest_unordered import unordered

# Custom libraries alphabetic order of main package.
import check_kinship


@pytest.fixture(scope="module", autouse=True)
def setup_test_path(tmp_path_factory):
    test_tmp_path = str(tmp_path_factory.mktemp("test")) + "/"
    open(str(test_tmp_path) + "/empty.txt", "a").close()
    return test_tmp_path


@pytest.fixture(scope="module", autouse=True)
def kinship_settings():
    return 0.177, 0.354

class TestNonEmptyExistingPath():
    def test_existing_dir(self, setup_test_path):
        file_or_dir = check_kinship.non_empty_existing_path(setup_test_path)
        assert file_or_dir

    def test_not_file_not_dir(self):
        fake_string = "fake_string"
        with pytest.raises(FileNotFoundError) as file_dir_error:
            check_kinship.non_empty_existing_path(fake_string)
        assert fake_string in str(file_dir_error.value)

    def test_empty_file(self, setup_test_path):
        with pytest.raises(OSError) as empty_error:
            check_kinship.non_empty_existing_path(setup_test_path + "empty.txt")
        assert f"File {setup_test_path}empty.txt is empty." in str(empty_error.value)

    def test_append_suffix(self, setup_test_path):
        dir_without_suffix = setup_test_path.rstrip("/")
        dir_with_suffix = check_kinship.non_empty_existing_path(dir_without_suffix)
        assert dir_without_suffix[-1] != "/"
        assert dir_with_suffix[-1] == "/"


@pytest.mark.parametrize("input_file,exp_dict_samples", [
    (
        "multi_family.ped",
        {
            "2024D00001": {'family': 'U000001', 'parents': [], 'children': ["2024D00003"]},
            "2024D00002": {'family': 'U000001', 'parents': [], 'children': ["2024D00003"]},
            "2024D00003": {'family': 'U000001', 'parents': ["2024D00001", "2024D00002"], 'children': []},
            "2024D00004": {'family': 'U000002', 'parents': [], 'children': ["2024D00006"]},
            "2024D00005": {'family': 'U000002', 'parents': [], 'children': ["2024D00006"]},
            "2024D00006": {'family': 'U000002', 'parents': ["2024D00004", "2024D00005"], 'children': []}
        },
    ),
    (
        "multi_siblings.ped",
        {
            "2024D00001": {'family': 'U000001', 'parents': [], 'children': ["2024D00003", "2024D00004"]},
            "2024D00002": {'family': 'U000001', 'parents': [], 'children': ["2024D00003", "2024D00004"]},
            "2024D00003": {'family': 'U000001', 'parents': ["2024D00001", "2024D00002"], 'children': []},
            "2024D00004": {'family': 'U000001', 'parents': ["2024D00001", "2024D00002"], 'children': []}
        },
    ),
    (
        "multi_unrelated_samples.ped",
        {
            "2024D00001": {'family': 'U000001', 'parents': [], 'children': []},
            "2024D00002": {'family': 'U000002', 'parents': [], 'children': []}
        },
    ),
    (
        "single_family.ped",
        {
            "2024D00001": {'family': 'U000001', 'parents': [], 'children': ["2024D00003"]},
            "2024D00002": {'family': 'U000001', 'parents': [], 'children': ["2024D00003"]},
            "2024D00003": {'family': 'U000001', 'parents': ["2024D00001", "2024D00002"], 'children': []}
        },
    ),
])
def test_parse_ped(input_file, exp_dict_samples, datadir):
    dict_samples = check_kinship.parse_ped(f"{datadir}/{input_file}")
    assert dict_samples.keys() == unordered(exp_dict_samples.keys())
    for sample, meta in dict_samples.items():
        assert meta.get("family") == exp_dict_samples.get(sample).get("family")
        assert meta.get("parents") == unordered(exp_dict_samples.get(sample).get("parents"))
        assert meta.get("children") == unordered(exp_dict_samples.get(sample).get("children"))


def test_read_kinship_trio(datadir, kinship_settings):
    kinship_min, kinship_max = kinship_settings
    df_out = check_kinship.read_kinship(f"{datadir}/trio.kinship", kinship_min, kinship_max)
    # Assert succeeded read_table
    assert not df_out.empty
    # Assert renaming columns
    assert set(["sample_1", "sample_2", "kinship"]) < set(df_out.columns)
    assert not set(["FID1", "FID2", "Kinship"]) <= set(df_out.columns)
    # Assert added columns
    assert set(["related", "type", "status", "thresholds", "message"]) < set(df_out.columns)
    # Assert threshold added as string
    assert df_out.thresholds.iloc[0] == f"{kinship_min},{kinship_max}" and isinstance(df_out.thresholds.iloc[0], str)


class TestCheckAndAnnotateKinship():
    @pytest.mark.parametrize("samples,file,related,type", [
        # Type parent_child with kinship values between min and max threshold
        (
            {
                "2024D00003-2024D00004": {'family': 'U000001', 'parents': ["2024D00001-2024D00002"], 'children': []},
                "2024D00001-2024D00002": {'family': 'U000001', 'parents': [], 'children': ["2024D00003-2024D00004"]},
            },
            "single_row_kin.kinship", True, 'parent_child'
        ),
        # Type parent_parent with kinship values lower than min threshold
        (
            {
                "2024D00003-2024D00004": {'family': 'U000001', 'parents': [], 'children': ["2024D00005-2024D00006"]},
                "2024D00001-2024D00002": {'family': 'U000001', 'parents': [], 'children': ["2024D00005-2024D00006"]},
            },
            "single_row.kinship", True, 'parent_parent'
        ),
        # Type sibling_sibling with kinship values between min and max threshold
        (
            {
                "2024D00003-2024D00004": {'family': 'U000001', 'parents': [], 'children': []},
                "2024D00001-2024D00002": {'family': 'U000001', 'parents': [], 'children': []},
            },
            "single_row_kin.kinship", True, 'sibling_sibling'
        ),
        # Type unrelated with kinship values lower than min threshold
        (
            {
                "2024D00003-2024D00004": {'family': 'U000001', 'parents': [], 'children': []},
                "2024D00001-2024D00002": {'family': 'U000002', 'parents': [], 'children': []},
            },
            "single_row.kinship", False, 'unrelated'
        ),
    ])
    def test_single_rows_ok(self, samples, file, related, type, datadir, kinship_settings):
        kinship_min, kinship_max = kinship_settings
        df_in = check_kinship.read_kinship(f"{datadir}/{file}", kinship_min, kinship_max)
        df_out = check_kinship.check_and_annotate_kinship(df_in, samples, kinship_min, kinship_max)
        assert df_out.shape == (1, 8)
        assert df_out.loc[0, "related"] == related
        assert df_out.loc[0, "type"] == type
        assert df_out.loc[0, "status"] == "OK"
        assert df_out.loc[0, "message"] == ''

    @pytest.mark.parametrize("samples,file,related,type,msg", [
        # Type parent_child with kinship values lower than min threshold
        (
            {
                "2024D00003-2024D00004": {'family': 'U000001', 'parents': ["2024D00001-2024D00002"], 'children': []},
                "2024D00001-2024D00002": {'family': 'U000001', 'parents': [], 'children': ["2024D00003-2024D00004"]},
            },
            "single_row.kinship", True, 'parent_child', '> 0.177 and < 0.354'
        ),
        # Type parent_child with kinship values higher than min threshold
        (
            {
                "2024D00003-2024D00004": {'family': 'U000001', 'parents': ["2024D00001-2024D00002"], 'children': []},
                "2024D00001-2024D00002": {'family': 'U000001', 'parents': [], 'children': ["2024D00003-2024D00004"]},
            },
            "single_row_duplo.kinship", True, 'parent_child', '> 0.177 and < 0.354'
        ),
        # Type parent_parent with kinship values between min and max threshold
        (
            {
                "2024D00003-2024D00004": {'family': 'U000001', 'parents': [], 'children': ["2024D00005-2024D00006"]},
                "2024D00001-2024D00002": {'family': 'U000001', 'parents': [], 'children': ["2024D00005-2024D00006"]},
            },
            "single_row_kin.kinship", True, 'parent_parent', '<= 0.177'
        ),
        # Type sibling_sibling with kinship values lower than min threshold
        (
            {
                "2024D00003-2024D00004": {'family': 'U000001', 'parents': [], 'children': []},
                "2024D00001-2024D00002": {'family': 'U000001', 'parents': [], 'children': []},
            },
            "single_row.kinship", True, 'sibling_sibling', '> 0.177 and < 0.354'
        ),
        # Type sibling_sibling with kinship values higher than max threshold
        (
            {
                "2024D00003-2024D00004": {'family': 'U000001', 'parents': [], 'children': []},
                "2024D00001-2024D00002": {'family': 'U000001', 'parents': [], 'children': []},
            },
            "single_row_duplo.kinship", True, 'sibling_sibling', '> 0.177 and < 0.354'
        ),
        # Type unrelated with kinship values between min and max threshold
        (
            {
                "2024D00003-2024D00004": {'family': 'U000001', 'parents': [], 'children': []},
                "2024D00001-2024D00002": {'family': 'U000002', 'parents': [], 'children': []},
            },
            "single_row_kin.kinship", False, 'unrelated', '<= 0.177'
        ),
    ])
    def test_single_rows_fail(self, samples, file, related, type, msg, datadir, kinship_settings):
        kinship_min, kinship_max = kinship_settings
        df_in = check_kinship.read_kinship(f"{datadir}/{file}", kinship_min, kinship_max)
        df_out = check_kinship.check_and_annotate_kinship(df_in, samples, kinship_min, kinship_max)
        assert df_out.shape == (1, 8)
        assert df_out.loc[0, "related"] == related
        assert df_out.loc[0, "type"] == type
        assert df_out.loc[0, "status"] == "FAIL"

        msg_out = df_out.loc[0, "message"]
        assert str(df_in.loc[0, "kinship"]) in msg_out
        for sample, meta in samples.items():
            assert sample in msg_out
            assert meta.get("family") in msg_out
        assert type in msg_out
        assert msg in msg_out

    @pytest.mark.parametrize("samples,file", [
        # OK: trio
        (
            {
                "2024D00005-2024D00006": {'family': 'U000001', 'parents': [], 'children': ["2024D00001-2024D00002"]},
                "2024D00003-2024D00004": {'family': 'U000001', 'parents': [], 'children': ["2024D00001-2024D00002"]},
                "2024D00001-2024D00002": {
                    'family': 'U000001', 'parents': ["2024D00003-2024D00004", "2024D00005-2024D00006"], 'children': []
                },
            },
            "trio.kinship"
        ),
    ])
    def test_multi_rows_trio_ok(self, samples, file, datadir, kinship_settings):
        kinship_min, kinship_max = kinship_settings
        df_in = check_kinship.read_kinship(f"{datadir}/{file}", kinship_min, kinship_max)
        df_out = check_kinship.check_and_annotate_kinship(df_in, samples, kinship_min, kinship_max)
        assert df_out.shape == (3, 8)
        assert "parent_child" in df_out.type.values
        assert "parent_parent" in df_out.type.values
        assert "OK" == df_out.status.unique()
        assert not df_out.message.unique()


class TestWriteKinship():
    @pytest.mark.parametrize("input_dict,error_comment", [
        (
            {
                "sample_1": ["2024D00001-2024D00002"],
                "sample_2": ["2024D00003-2024D00004"],
                "status": ["OK"],
                "thresholds": ["0.177, 0.354"]
            },
            "# No kinship errors found."
        ),
        (
            {
                "sample_1": ["2024D00001-2024D00002"],
                "sample_2": ["2024D00003-2024D00004"],
                "status": ["FAIL"],
                "thresholds": ["0.177, 0.354"]
            },
            "# WARNING: Kinship errors found."
        ),
    ])
    def test_file(self, input_dict, error_comment, setup_test_path):
        df_in = DataFrame(input_dict)
        prefix = f"fake_prefix_{input_dict['status']}"
        check_kinship.write_kinship(df_in, setup_test_path, prefix)

        output_file = Path(f"{setup_test_path}/{prefix}.kinship_check.out")
        assert output_file.exists() and output_file.is_file()
        with open(output_file) as file:
            file_content = file.read().rstrip().splitlines()
        assert len(file_content) == 4
        assert file_content[0] == error_comment
        assert file_content[1] == "# Used kinship check settings: 0.177, 0.354"
        assert file_content[2].split("\t") == df_in.columns.to_list()
        assert file_content[3].split("\t") == df_in.iloc[0].tolist()

    def test_stdout(self, kinship_settings, capsys):
        kinship_min, kinship_max = kinship_settings
        df_in = DataFrame({
            "sample_1": ["2024D00001-2024D00002"],
            "sample_2": ["2024D00003-2024D00004"],
            "status": ["OK"],
            "thresholds": [f"{kinship_min}, {kinship_max}"],
        })

        check_kinship.write_kinship(df_kinship_out=df_in, output_path=None, output_prefix=None)
        out, err = capsys.readouterr()
        lst_lines = out.rstrip().splitlines()
        # Check number of expected lines
        assert len(lst_lines) == 4
        # Check comments / headers
        assert lst_lines[0] == "# No kinship errors found."
        assert lst_lines[1] == "# Used kinship check settings: 0.177, 0.354"
        # Check table
        assert lst_lines[2].split("\t") == df_in.columns.to_list()
        assert lst_lines[3].split("\t") == df_in.iloc[0].tolist()
