import { createSlice, PayloadAction } from '@reduxjs/toolkit'

import { convertToParamMap, getChildParam } from 'store/utils/param/ParamUtils'
import {
  addFlowElementNode,
  deleteFlowElements,
  deleteFlowElementsById,
} from '../FlowElement/FlowElementSlice'
import { isNodeData } from '../FlowElement/FlowElementUtils'
import { NODE_TYPE_SET } from '../FlowElement/FlowElementType'
import { importExperimentByUid } from '../Experiments/ExperimentsActions'
import { getAlgoParams } from './AlgorithmNodeActions'
import { ALGORITHM_NODE_SLICE_NAME, AlgorithmNode } from './AlgorithmNodeType'
import { isAlgorithmNodePostData } from 'api/run/RunUtils'

const initialState: AlgorithmNode = {}

export const algorithmNodeSlice = createSlice({
  name: ALGORITHM_NODE_SLICE_NAME,
  initialState,
  reducers: {
    updateParam: (
      state,
      action: PayloadAction<{
        nodeId: string
        // paramKey: string
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
      .addCase(addFlowElementNode, (state, action) => {
        const { node, algoNodeInfo } = action.payload
        if (
          node.data?.type === NODE_TYPE_SET.ALGORITHM &&
          algoNodeInfo != null
        ) {
          state[node.id] = {
            ...algoNodeInfo,
            params: null,
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
              }
            }
          })
        return newState
      })
  },
})

export const { updateParam } = algorithmNodeSlice.actions
export default algorithmNodeSlice.reducer
