import pytest

from studio.app.common.core.workflow.workflow_params import get_typecheck_params


def test_unknown_param_raises_named_yaml_error():
    with pytest.raises(KeyError) as e:
        get_typecheck_params({"transpose": {"type": "child", "value": False}}, "dpca")
    detail = e.value.args[0]
    assert detail.startswith("Workflow yaml error, see FAQ")
    assert "'transpose'" in detail
    assert "dpca" in detail
