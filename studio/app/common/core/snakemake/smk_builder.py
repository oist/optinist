from studio.app.common.core.snakemake.smk import Rule


class RuleBuilder:
    def __init__(self) -> None:
        self._input = None
        self._return_arg = None
        self._params = None
        self._output = None
        self._type = None
        self._nwbfile = None
        self._hdf5Path = None
        self._matPath = None
        self._path = None

    def set_input(self, input, workspace_id=None) -> "RuleBuilder":
        if workspace_id:
            if isinstance(input, list):
                self._input = [f"{workspace_id}/{i}" for i in input]
            else:
                self._input = f"{workspace_id}/{input}"
        else:
            self._input = input
        return self

    def set_return_arg(self, return_arg) -> "RuleBuilder":
        self._return_arg = return_arg
        return self

    def set_params(self, params) -> "RuleBuilder":
        self._params = params
        return self

    def set_output(self, output) -> "RuleBuilder":
        self._output = output
        return self

    def set_type(self, type) -> "RuleBuilder":
        self._type = type
        return self

    def set_nwbfile(self, nwbfile) -> "RuleBuilder":
        self._nwbfile = nwbfile
        return self

    def set_hdf5Path(self, hdf5Path) -> "RuleBuilder":
        self._hdf5Path = hdf5Path
        return self

    def set_matPath(self, matPath) -> "RuleBuilder":
        self._matPath = matPath
        return self

    def set_path(self, path) -> "RuleBuilder":
        self._path = path
        return self

    def build(self) -> Rule:
        return Rule(
            input=self._input,
            return_arg=self._return_arg,
            params=self._params,
            output=self._output,
            type=self._type,
            nwbfile=self._nwbfile,
            hdf5Path=self._hdf5Path,
            matPath=self._matPath,
            path=self._path,
        )
