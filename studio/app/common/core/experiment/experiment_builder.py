from studio.app.common.core.experiment.experiment import ExptConfig


class ExptConfigBuilder:
    def __init__(self):
        self._workspace_id = None
        self._unique_id = None
        self._name = None
        self._started_at = None
        self._finished_at = None
        self._success = None
        self._hasNWB = False
        self._function = {}
        self._nwbfile = None
        self._snakemake = None

    def set_config(self, config: ExptConfig) -> "ExptConfigBuilder":
        self._workspace_id = config.workspace_id
        self._unique_id = config.unique_id
        self._name = config.name
        self._started_at = config.started_at
        self._finished_at = config.finished_at
        self._success = config.success
        self._hasNWB = config.hasNWB
        self._function = config.function
        self._nwbfile = config.nwb
        self._snakemake = config.snakemake
        return self

    def set_workspace_id(self, workspace_id) -> "ExptConfigBuilder":
        self._workspace_id = workspace_id
        return self

    def set_unique_id(self, unique_id) -> "ExptConfigBuilder":
        self._unique_id = unique_id
        return self

    def set_name(self, name) -> "ExptConfigBuilder":
        self._name = name
        return self

    def set_started_at(self, timestamp) -> "ExptConfigBuilder":
        self._started_at = timestamp
        return self

    def set_success(self, success: str) -> "ExptConfigBuilder":
        self._success = success
        return self

    def set_hasNWB(self, hasNWB) -> "ExptConfigBuilder":
        self._hasNWB = hasNWB
        return self

    def set_function(self, function) -> "ExptConfigBuilder":
        self._function = function
        return self

    def set_nwbfile(self, nwbfile) -> "ExptConfigBuilder":
        self._nwbfile = nwbfile
        return self

    def set_snakemake(self, snakemake) -> "ExptConfigBuilder":
        self._snakemake = snakemake
        return self

    def build(self) -> ExptConfig:
        return ExptConfig(
            workspace_id=self._workspace_id,
            unique_id=self._unique_id,
            name=self._name,
            started_at=self._started_at,
            finished_at=self._finished_at,
            success=self._success,
            hasNWB=self._hasNWB,
            function=self._function,
            nwb=self._nwbfile,
            snakemake=self._snakemake,
        )
