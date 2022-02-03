import { createSlice, PayloadAction } from '@reduxjs/toolkit'

import { convertToParamMap, getChildParam } from 'store/utils/param/ParamUtils'
import {
  addFlowElementNode,
  deleteFlowElements,
  deleteFlowElementsById,
} from '../FlowElement/FlowElementSlice'
import { isNodeData } from '../FlowElement/FlowElementUtils'
import { NODE_TYPE_SET } from '../FlowElement/FlowElementType'
import { getAlgoParams } from './AlgorithmNodeActions'
import { ALGORITHM_NODE_SLICE_NAME, AlgorithmNode } from './AlgorithmNodeType'

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
    setSelectedOutputKey: (
      state,
      action: PayloadAction<{
        nodeId: string
        selectedOutputKey: string
      }>,
    ) => {
      const { nodeId, selectedOutputKey } = action.payload
      state[nodeId].selectedOutputKey = selectedOutputKey
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
            selectedOutputKey: null,
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
  },
})

export const { updateParam, setSelectedOutputKey } = algorithmNodeSlice.actions
export default algorithmNodeSlice.reducer
