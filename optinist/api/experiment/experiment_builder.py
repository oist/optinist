from optinist.api.experiment.experiment import ExptConfig


class ExptConfigBuilder:
    def __init__(self):
        self._started_at = None
        self._finished_at = None
        self._success = None
        self._name = None
        self._unique_id = None
        self._hasNWB = False
        self._function = {}
        self._nodeDict = None
        self._edgeDict = None

    def set_config(self, config: ExptConfig) -> "ExptConfigBuilder":
        self._started_at = config.started_at
        self._finished_at = config.finished_at
        self._success = config.success
        self._name = config.name
        self._unique_id = config.unique_id
        self._hasNWB = config.hasNWB
        self._function = config.function
        self._nodeDict = config.nodeDict
        self._edgeDict = config.edgeDict
        return self

    def set_started_at(self, timestamp) -> "ExptConfigBuilder":
        self._started_at = timestamp
        return self

    def set_success(self, success: str) -> "ExptConfigBuilder":
        self._success = success
        return self

    def set_name(self, name) -> "ExptConfigBuilder":
        self._name = name
        return self

    def set_unique_id(self, unique_id) -> "ExptConfigBuilder":
        self._unique_id = unique_id
        return self

    def set_hasNWB(self, hasNWB) -> "ExptConfigBuilder":
        self._hasNWB = hasNWB
        return self

    def set_function(self, function) -> "ExptConfigBuilder":
        self._function = function
        return self

    def set_nodeDict(self, nodeDict) -> "ExptConfigBuilder":
        self._nodeDict = nodeDict
        return self

    def set_edgeDict(self, edgeDict) -> "ExptConfigBuilder":
        self._edgeDict = edgeDict
        return self

    def build(self) -> ExptConfig:
        return ExptConfig(
            started_at=self._started_at,
            finished_at=self._finished_at,
            success=self._success,
            name=self._name,
            unique_id=self._unique_id,
            hasNWB=self._hasNWB,
            function=self._function,
            nodeDict=self._nodeDict,
            edgeDict=self._edgeDict,
        )
