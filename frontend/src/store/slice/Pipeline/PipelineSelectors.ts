import { RunPostData } from 'api/Run/Run'
import { RootState } from 'store/store'
import { selectElementListForRun } from '../FlowElement/FlowElementSelectors'
import { selectNwbParams } from '../NWB/NWBSelectors'
import { selectSnakemakeParams } from '../Snakemake/SnakemakeSelectors'
import {
  NodeResult,
  NodeResultPending,
  NodeResultSuccess,
  RUN_STATUS,
  StartedPipeline,
} from './PipelineType'
import {
  isNodeResultPending,
  isStartedPipeline,
  isNodeResultSuccess,
} from './PipelineUtils'

export const selectPipelineLatestUid = (state: RootState) => {
  const history = state.pipeline.uidHistory
  if (history.length > 0) {
    return history.slice(-1)[0]
  } else {
    return undefined
  }
}

export const selectRunPostData = (state: RootState) => {
  const nwbParam = selectNwbParams(state)
  const snakemakeParam = selectSnakemakeParams(state)
  const elementListForRun = selectElementListForRun(state)
  const runPostData: RunPostData = {
    nwbParam,
    snakemakeParam,
    ...elementListForRun,
  }
  return runPostData
}

export const selectStartedPipeline = (uid: string) => (state: RootState) => {
  const pipeline = Object.values(state.pipeline.pipelines).find(
    (value) => isStartedPipeline(value) && value.uid === uid,
  )
  if (pipeline != null) {
    return pipeline as StartedPipeline
  } else {
    throw new Error(`invalid uid: ${uid}`)
  }
}

export const selectRunResultPendingList =
  (uid: string) => (state: RootState) => {
    const pipeline = selectStartedPipeline(uid)(state)
    return Object.values(pipeline.runResult).filter(isNodeResultPending)
  }

export const selectRunResultPendingNodeIdList =
  (uid: string) => (state: RootState) => {
    const pipeline = selectStartedPipeline(uid)(state)
    return Object.entries(pipeline.runResult)
      .map(([nodeId, nodeResult]) => ({ nodeId, nodeResult }))
      .filter(isNodeResultPendingAndNodeId)
      .map(({ nodeId }) => nodeId)
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

export const selectPipelineStatus = (uid: string) => (state: RootState) => {
  const pipeline = selectStartedPipeline(uid)(state)
  return pipeline.status
}

export const selectPipelineIsCanceled = (uid: string) => (state: RootState) => {
  const pipeline = selectStartedPipeline(uid)(state)
  return pipeline.status === RUN_STATUS.CANCELED
}

export const selectPipelineIsStartedSuccess =
  (uid: string) => (state: RootState) => {
    const pipeline = selectStartedPipeline(uid)(state)
    return pipeline.status === RUN_STATUS.START_SUCCESS
  }

export const selectPipelineRunResult = (uid: string) => (state: RootState) => {
  const pipeline = selectStartedPipeline(uid)(state)
  return pipeline.runResult
}

export const selectPipelineNodeResultSuccessList =
  (uid: string) => (state: RootState) => {
    const pipeline = selectStartedPipeline(uid)(state)
    return Object.entries(pipeline.runResult)
      .map(([nodeId, nodeResult]) => {
        return {
          nodeId,
          nodeResult,
        }
      })
      .filter(isNodeResultSuccessAndNodeId)
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
  (uid: string, nodeId: string) => (state: RootState) => {
    const pipeline = selectStartedPipeline(uid)(state)
    if (Object.keys(pipeline.runResult).includes(nodeId)) {
      return pipeline.runResult[nodeId].status
    } else {
      return null
    }
  }

export const selectPipelineNodeResultMessage =
  (uid: string, nodeId: string) => (state: RootState) => {
    const pipeline = selectStartedPipeline(uid)(state)
    if (Object.keys(pipeline.runResult).includes(nodeId)) {
      return pipeline.runResult[nodeId].message
    } else {
      return null
    }
  }

export const selectPipelineNodeResultOutputKeyList =
  (uid: string, nodeId: string) => (state: RootState) => {
    const pipeline = selectStartedPipeline(uid)(state)
    if (Object.keys(pipeline.runResult).includes(nodeId)) {
      const nodeResult = pipeline.runResult[nodeId]
      if (isNodeResultSuccess(nodeResult)) {
        return Object.keys(nodeResult.outputPaths)
      }
    }
    return []
  }
