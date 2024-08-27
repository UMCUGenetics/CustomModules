#!/usr/bin/env python
# Import statements, alphabetic order of main package.

# Third party libraries alphabetic order of main package.
import pytest
from pytest_unordered import unordered

# Custom libraries alphabetic order of main package.
from CustomModules.Utils.parse_child_from_fulltrio import parse_ped


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
def test_parse_ped(input_file, exp_dict_samples, shared_datadir):
    dict_samples = parse_ped(open(f"{shared_datadir}/{input_file}", "r"))
    assert dict_samples.keys() == unordered(exp_dict_samples.keys())
    for sample, meta in dict_samples.items():
        assert meta.get("family") == exp_dict_samples.get(sample).get("family")
        assert meta.get("parents") == unordered(exp_dict_samples.get(sample).get("parents"))
        assert meta.get("children") == unordered(exp_dict_samples.get(sample).get("children"))
