import { Node } from 'react-flow-renderer'

import {
  AlgorithmNodePostData,
  InputNodePostData,
  NodeDict,
  RunPostData,
} from 'api/run/Run'

import { RootState } from 'store/store'

import {
  selectAlgorithmFunctionPath,
  selectAlgorithmIsUpdated,
  selectAlgorithmName,
  selectAlgorithmParams,
} from 'store/slice/AlgorithmNode/AlgorithmNodeSelectors'
import {
  selectEdgeListForRun,
  selectFlowElements,
} from 'store/slice/FlowElement/FlowElementSelectors'
import {
  isAlgorithmNodeData,
  isNodeData,
} from 'store/slice/FlowElement/FlowElementUtils'
import { selectNwbParams } from 'store/slice/NWB/NWBSelectors'
import { selectPipelineNodeResultStatus } from 'store/slice/Pipeline/PipelineSelectors'
import { NODE_RESULT_STATUS } from 'store/slice/Pipeline/PipelineType'
import { selectSnakemakeParams } from 'store/slice/Snakemake/SnakemakeSelectors'
import { NODE_TYPE_SET } from 'store/slice/FlowElement/FlowElementType'
import {
  selectInputNodeFileType,
  selectInputNodeHDF5Path,
  selectInputNodeParam,
  selectInputNodeSelectedFilePath,
} from 'store/slice/InputNode/InputNodeSelectors'

export const selectRunPostData = (state: RootState) => {
  const nwbParam = selectNwbParams(state)
  const snakemakeParam = selectSnakemakeParams(state)
  const edgeListForRun = selectEdgeListForRun(state)
  const nodeDictForRun = selectNodeDictForRun(state)
  const forceRunList = selectForceRunList(state)
  const runPostData: Omit<RunPostData, 'name'> = {
    nwbParam,
    snakemakeParam,
    edgeList: edgeListForRun,
    nodeDict: nodeDictForRun,
    forceRunList,
  }
  return runPostData
}

/**
 * 前回の結果で、エラーまたはParamに変更があるnodeのリストを返す
 */
const selectForceRunList = (state: RootState) => {
  const elements = selectFlowElements(state)
  return elements
    .filter(isAlgorithmNodeData)
    .filter((node) => {
      const isUpdated = selectAlgorithmIsUpdated(node.id)(state)
      const status = selectPipelineNodeResultStatus(node.id)(state)
      return isUpdated || status === NODE_RESULT_STATUS.ERROR
    })
    .map((node) => ({
      nodeId: node.id,
      name: selectAlgorithmName(node.id)(state),
    }))
}

const selectNodeDictForRun = (state: RootState): NodeDict => {
  const elements = selectFlowElements(state)
  const nodeDict: NodeDict = {}
  elements.filter(isNodeData).forEach((node) => {
    if (isAlgorithmNodeData(node)) {
      const param = selectAlgorithmParams(node.id)(state) ?? {}
      const functionPath = selectAlgorithmFunctionPath(node.id)(state)
      const algorithmNodePostData: Node<AlgorithmNodePostData> = {
        ...node,
        data: {
          ...node.data,
          label: node.data?.label ?? '',
          type: NODE_TYPE_SET.ALGORITHM,
          path: functionPath,
          param,
        },
      }
      nodeDict[node.id] = algorithmNodePostData
    } else {
      const filePath = selectInputNodeSelectedFilePath(node.id)(state)
      const fileType = selectInputNodeFileType(node.id)(state)
      const param = selectInputNodeParam(node.id)(state)
      const hdf5Path = selectInputNodeHDF5Path(node.id)(state)
      const inputNodePosyData: Node<InputNodePostData> = {
        ...node,
        data: {
          ...node.data,
          label: node.data?.label ?? '',
          type: NODE_TYPE_SET.INPUT,
          path: filePath ?? '',
          param,
          hdf5Path: hdf5Path,
          fileType,
        },
      }
      nodeDict[node.id] = inputNodePosyData
    }
  })
  return nodeDict
}
