import { RootState } from 'store/store'

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

export const selectPipelineLatestUid = (state: RootState) => {
  return state.pipeline.currentPipeline?.uid
}

export const selectStartedPipeline = (state: RootState) => {
  return state.pipeline.run
}

export const selectPipelineRunBtn = (state: RootState) => {
  return state.pipeline.runBtn
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

const selectPipelineNodeResultOutputPaths =
  (nodeId: string) => (state: RootState) => {
    const pipeline = selectStartedPipeline(state)
    if (isStartedPipeline(pipeline)) {
      const nodeResult = pipeline.runResult[nodeId]
      if (
        Object.keys(pipeline.runResult).includes(nodeId) &&
        isNodeResultSuccess(nodeResult)
      ) {
        return nodeResult.outputPaths
      }
    }
    throw new Error(`key error. nodeId:${nodeId}`)
  }

export const selectPipelineNodeResultOutputFilePath =
  (nodeId: string, outputKey: string) => (state: RootState) => {
    const outputPaths = selectPipelineNodeResultOutputPaths(nodeId)(state)
    if (Object.keys(outputPaths).includes(outputKey)) {
      return outputPaths[outputKey].path
    } else {
      throw new Error(`key error. outputKey:${outputKey}`)
    }
  }

export const selectPipelineNodeResultOutputFileDataType =
  (nodeId: string, outputKey: string) => (state: RootState) => {
    const outputPaths = selectPipelineNodeResultOutputPaths(nodeId)(state)
    if (Object.keys(outputPaths).includes(outputKey)) {
      return outputPaths[outputKey].type
    } else {
      throw new Error(`key error. outputKey:${outputKey}`)
    }
  }
