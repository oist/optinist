import {
  createStructuredFileNode,
  FileNodeConfig,
} from "components/Workspace/FlowChart/FlowChartNode/BaseStructuredFileNode"
import { FILE_TYPE_SET } from "config/fileTypes.config"
import { getHDF5Tree } from "store/slice/HDF5/HDF5Action"
import {
  selectHDF5IsLoading,
  selectHDF5Nodes,
} from "store/slice/HDF5/HDF5Selectors"
import {
  selectHDF5InputNodeSelectedFilePath,
  selectInputNodeHDF5Path,
} from "store/slice/InputNode/InputNodeSelectors"
import { setInputNodeHDF5Path } from "store/slice/InputNode/InputNodeSlice"

const hdf5Config: FileNodeConfig = {
  fileType: FILE_TYPE_SET.HDF5,
  handleId: "hdf5",
  handleType: "HDF5Data",
  treeKeyPrefix: "hdf5tree",
  selectFilePath: selectHDF5InputNodeSelectedFilePath,
  selectStructurePath: selectInputNodeHDF5Path,
  setStructurePath: setInputNodeHDF5Path,
  getTree: getHDF5Tree,
  selectTree: selectHDF5Nodes,
  selectIsLoading: selectHDF5IsLoading,
}

export const HDF5FileNode = createStructuredFileNode(hdf5Config)
