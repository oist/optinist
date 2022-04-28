import os

import pytest

from optinist.api.dir_path import DIRPATH
from optinist.api.snakemake.smk import Rule


@pytest.mark.parametrize(
    ("original_rule", "cluster_converted_rule"),
    [
        (
            Rule(
                rule_file="test",
                input=[os.path.join(DIRPATH.INPUT_DIR, "Sue", "Sue.tiff")],
                return_arg=(),
                params={},
                output=os.path.join(
                    DIRPATH.OUTPUT_DIR, "unique_id", "algorithm", "Sue.pkl"
                ),
                type="image",
            ),
            Rule(
                rule_file="test",
                input=[os.path.join(DIRPATH.CLUSTER_INPUT_DIR, "Sue", "Sue.tiff")],
                return_arg=(),
                params={},
                output=os.path.join(
                    DIRPATH.CLUSTER_OUTPUT_DIR, "unique_id", "algorithm", "Sue.pkl"
                ),
                type="image",
            ),
        ),
        (
            Rule(
                rule_file="test",
                input=[os.path.join(DIRPATH.OUTPUT_DIR, "Sue", "Sue.tiff")],
                return_arg=(),
                params={},
                output=os.path.join(
                    DIRPATH.OUTPUT_DIR, "unique_id", "algorithm", "Sue.pkl"
                ),
                type="image",
            ),
            Rule(
                rule_file="test",
                input=[os.path.join(DIRPATH.CLUSTER_OUTPUT_DIR, "Sue", "Sue.tiff")],
                return_arg=(),
                params={},
                output=os.path.join(
                    DIRPATH.CLUSTER_OUTPUT_DIR, "unique_id", "algorithm", "Sue.pkl"
                ),
                type="image",
            ),
        ),
    ],
)
def test_snakemake_rule_object_conversion(original_rule, cluster_converted_rule):
    assert original_rule.to_cluster_rule() == cluster_converted_rule
