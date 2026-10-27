import {
  createStructuredFileNode,
  FileNodeConfig,
} from "components/Workspace/FlowChart/FlowChartNode/BaseStructuredFileNode"
import { FILE_TYPE_SET } from "config/fileTypes.config"
import {
  selectInputNodeMatlabPath,
  selectMatlabInputNodeSelectedFilePath,
} from "store/slice/InputNode/InputNodeSelectors"
import { setInputNodeMatlabPath } from "store/slice/InputNode/InputNodeSlice"
import { getMatlabTree } from "store/slice/Matlab/MatlabAction"
import {
  selectMatlabIsLoading,
  selectMatlabNodes,
} from "store/slice/Matlab/MatlabSelectors"

const matlabConfig: FileNodeConfig = {
  fileType: FILE_TYPE_SET.MATLAB,
  handleId: "matlab",
  handleType: "MatlabData",
  treeKeyPrefix: "matlabtree",
  selectFilePath: selectMatlabInputNodeSelectedFilePath,
  selectStructurePath: selectInputNodeMatlabPath,
  setStructurePath: setInputNodeMatlabPath,
  getTree: getMatlabTree,
  selectTree: selectMatlabNodes,
  selectIsLoading: selectMatlabIsLoading,
}

export const MatlabFileNode = createStructuredFileNode(matlabConfig)
