import { Node } from 'react-flow-renderer'

import {
  AlgorithmNodePostData,
  InputNodePostData,
  NodePostDataType,
  RunPostData,
} from 'api/run/Run'

import { RootState } from 'store/store'
import {
  selectEdgeListForRun,
  selectFlowElements,
} from '../FlowElement/FlowElementSelectors'
import { NODE_TYPE_SET } from '../FlowElement/FlowElementType'
import {
  isAlgorithmNodeData,
  isNodeData,
} from '../FlowElement/FlowElementUtils'
import { selectNwbParams } from '../NWB/NWBSelectors'
import { selectSnakemakeParams } from '../Snakemake/SnakemakeSelectors'
import {
  NodeResult,
  NodeResultPending,
  NodeResultSuccess,
  RUN_STATUS,
} from './PipelineType'
import {
  isNodeResultPending,
  isStartedPipeline,
  isNodeResultSuccess,
} from './PipelineUtils'

import {
  selectAlgorithmFunctionPath,
  selectAlgorithmParams,
} from '../AlgorithmNode/AlgorithmNodeSelectors'
import {
  selectInputNodeFileType,
  selectInputNodeHDF5Path,
  selectInputNodeParam,
  selectInputNodeSelectedFilePath,
} from '../InputNode/InputNodeSelectors'

export const selectPipelineLatestUid = (state: RootState) => {
  return state.pipeline.currentPipeline?.uid
}

export const selectRunPostData = (state: RootState) => {
  const nwbParam = selectNwbParams(state)
  const snakemakeParam = selectSnakemakeParams(state)
  const edgeListForRun = selectEdgeListForRun(state)
  const nodePostDataList = selectNodePostDataListForRun(state)
  const runPostData: RunPostData = {
    nwbParam,
    snakemakeParam,
    edgeList: edgeListForRun,
    nodeList: nodePostDataList,
  }
  return runPostData
}

export const selectNodePostDataListForRun = (
  state: RootState,
): Node<NodePostDataType>[] => {
  const elements = selectFlowElements(state)
  const nodeList = elements.filter(isNodeData).map((node) => {
    console.log('node:', node)
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
      return algorithmNodePostData
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
      return inputNodePosyData
    }
  })
  return nodeList
}

export const selectStartedPipeline = (state: RootState) => {
  return state.pipeline.run
}

export const selectRunResultPendingList = (state: RootState) => {
  const pipeline = selectStartedPipeline(state)
  if (isStartedPipeline(pipeline)) {
    return Object.values(pipeline.runResult).filter(isNodeResultPending)
  } else {
    return []
  }
}

export const selectRunResultPendingNodeIdList = (state: RootState) => {
  const pipeline = selectStartedPipeline(state)
  if (isStartedPipeline(pipeline)) {
    return Object.entries(pipeline.runResult)
      .map(([nodeId, nodeResult]) => ({ nodeId, nodeResult }))
      .filter(isNodeResultPendingAndNodeId)
      .map(({ nodeId }) => nodeId)
  } else {
    return []
  }
}

function isNodeResultPendingAndNodeId(arg: {
  nodeId: string
  nodeResult: NodeResult
}): arg is {
  nodeId: string
  nodeResult: NodeResultPending
} {
  return isNodeResultPending(arg.nodeResult)
}

export const selectPipelineStatus = (state: RootState) => {
  const pipeline = selectStartedPipeline(state)
  return pipeline.status
}

export const selectPipelineIsCanceled = (state: RootState) => {
  const pipeline = selectStartedPipeline(state)
  return pipeline.status === RUN_STATUS.CANCELED
}

export const selectPipelineIsStartedSuccess = (state: RootState) => {
  const pipeline = selectStartedPipeline(state)
  return pipeline.status === RUN_STATUS.START_SUCCESS
}

// export const selectPipelineRunResult = (state: RootState) => {
//   const pipeline = selectStartedPipeline(state)
//   if (isStartedPipeline(pipeline)) {
//     return pipeline.runResult
//   } else {
//     throw new Error("todo")
//   }
// }

export const selectPipelineNodeResultSuccessList = (state: RootState) => {
  const pipeline = selectStartedPipeline(state)
  if (isStartedPipeline(pipeline)) {
    return Object.entries(pipeline.runResult)
      .map(([nodeId, nodeResult]) => {
        return {
          nodeId,
          nodeResult,
        }
      })
      .filter(isNodeResultSuccessAndNodeId)
  } else {
    return []
  }
}

// selectPipelineNodeResultSuccessListの返り値の型を正しく認識させるためだけに作った
function isNodeResultSuccessAndNodeId(arg: {
  nodeId: string
  nodeResult: NodeResult
}): arg is {
  nodeId: string
  nodeResult: NodeResultSuccess
} {
  return isNodeResultSuccess(arg.nodeResult)
}

export const selectPipelineNodeResultStatus =
  (nodeId: string) => (state: RootState) => {
    const pipeline = selectStartedPipeline(state)
    if (isStartedPipeline(pipeline)) {
      if (Object.keys(pipeline.runResult).includes(nodeId)) {
        return pipeline.runResult[nodeId].status
      }
    }
    return null
  }

export const selectPipelineNodeResultMessage =
  (nodeId: string) => (state: RootState) => {
    const pipeline = selectStartedPipeline(state)
    if (isStartedPipeline(pipeline)) {
      if (Object.keys(pipeline.runResult).includes(nodeId)) {
        return pipeline.runResult[nodeId].message
      }
    }
    return null
  }

export const selectPipelineNodeResultOutputKeyList =
  (nodeId: string) => (state: RootState) => {
    const pipeline = selectStartedPipeline(state)
    if (isStartedPipeline(pipeline)) {
      const nodeResult = pipeline.runResult[nodeId]
      if (
        Object.keys(pipeline.runResult).includes(nodeId) &&
        isNodeResultSuccess(nodeResult)
      ) {
        return Object.keys(nodeResult.outputPaths)
      }
    }
    return []
  }
