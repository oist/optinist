import { createSlice, isAnyOf, PayloadAction } from '@reduxjs/toolkit'

import { convertToParamMap, getChildParam } from 'utils/param/ParamUtils'
import {
  deleteFlowElements,
  deleteFlowElementsById,
} from '../FlowElement/FlowElementSlice'
import { isNodeData } from '../FlowElement/FlowElementUtils'
import { NODE_TYPE_SET } from '../FlowElement/FlowElementType'
import { importExperimentByUid } from '../Experiments/ExperimentsActions'
import { getAlgoParams } from './AlgorithmNodeActions'
import { ALGORITHM_NODE_SLICE_NAME, AlgorithmNode } from './AlgorithmNodeType'
import { isAlgorithmNodePostData } from 'api/run/RunUtils'
import { run, runByCurrentUid } from '../Pipeline/PipelineActions'
import { addAlgorithmNode } from '../FlowElement/FlowElementActions'

const initialState: AlgorithmNode = {}

export const algorithmNodeSlice = createSlice({
  name: ALGORITHM_NODE_SLICE_NAME,
  initialState,
  reducers: {
    updateParam: (
      state,
      action: PayloadAction<{
        nodeId: string
        path: string
        newValue: unknown
      }>,
    ) => {
      const { nodeId, path, newValue } = action.payload
      const param = state[nodeId].params
      if (param != null) {
        const target = getChildParam(path, param)
        if (target != null) {
          target.value = newValue
          state[nodeId].isUpdated = true
        }
      }
    },
  },
  extraReducers: (builder) => {
    builder
      .addCase(getAlgoParams.fulfilled, (state, action) => {
        const { nodeId } = action.meta.arg
        state[nodeId].params = convertToParamMap(action.payload)
      })
      .addCase(addAlgorithmNode.fulfilled, (state, action) => {
        const { node, functionPath, name } = action.meta.arg
        const params = action.payload
        if (node.data?.type === NODE_TYPE_SET.ALGORITHM) {
          state[node.id] = {
            functionPath,
            name,
            params: convertToParamMap(params),
            isUpdated: false,
          }
        }
      })
      .addCase(deleteFlowElements, (state, action) => {
        action.payload
          .filter((node) => isNodeData(node))
          .forEach((node) => {
            if (node.data?.type === NODE_TYPE_SET.ALGORITHM) {
              delete state[node.id]
            }
          })
      })
      .addCase(deleteFlowElementsById, (state, action) => {
        if (Object.keys(state).includes(action.payload)) {
          delete state[action.payload]
        }
      })
      .addCase(importExperimentByUid.fulfilled, (_, action) => {
        const newState: AlgorithmNode = {}
        action.payload.nodeList
          .filter(isAlgorithmNodePostData)
          .forEach((node) => {
            if (node.data != null) {
              newState[node.id] = {
                name: node.data.label,
                functionPath: node.data.path,
                params: node.data.param,
                isUpdated: false,
              }
            }
          })
        return newState
      })
      .addMatcher(
        isAnyOf(run.fulfilled, runByCurrentUid.fulfilled),
        (state, action) => {
          const runPostData = action.meta.arg.runPostData
          runPostData.nodeList
            .filter(isAlgorithmNodePostData)
            .forEach((node) => {
              state[node.id].isUpdated = false
            })
        },
      )
  },
})

export const { updateParam } = algorithmNodeSlice.actions
export default algorithmNodeSlice.reducer
