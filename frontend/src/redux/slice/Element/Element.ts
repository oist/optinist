import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { Elements, Node } from 'react-flow-renderer'
import { initialElements, INITIAL_ALGO_ELEMENT_ID } from 'const/flowchart'
import { NodeData, NODE_DATA_TYPE_SET } from 'const/NodeData'
import { isNodeData } from 'utils/ElementUtils'
import {
  clickNode,
  runPipeline,
  stopPipeline,
  getAlgoParams,
} from './ElementAction'
import { Element, ELEMENT_SLICE_NAME, RUN_STATUS } from './ElementType'
import { uploadImageFile } from '../ImageIndex/ImageIndexAction'

const initialState: Element = {
  flowElements: initialElements,
  clickedNodeId: null,
  currentAlgoId: INITIAL_ALGO_ELEMENT_ID,
  algoParams: {},
  runStatus: RUN_STATUS.STOPPED,
}

export const elementSlice = createSlice({
  name: ELEMENT_SLICE_NAME,
  initialState,
  reducers: {
    setFlowElements: (state, action: PayloadAction<Elements>) => {
      state.flowElements = action.payload
    },
    addFlowElement: (state, action: PayloadAction<Node<NodeData>>) => {
      state.flowElements.push(action.payload)
    },
    updateParam: (
      state,
      action: PayloadAction<{ paramKey: string; newValue: unknown }>,
    ) => {
      const { paramKey, newValue } = action.payload
      state.algoParams[state.currentAlgoId].param[paramKey] = newValue
    },
  },
  extraReducers: (builder) => {
    builder
      .addCase(clickNode, (state, action) => {
        state.clickedNodeId = action.payload.id
        if (action.payload.type === NODE_DATA_TYPE_SET.ALGO) {
          state.currentAlgoId = action.payload.id
        }
      })
      .addCase(getAlgoParams.fulfilled, (state, action) => {
        const { id, algoName } = action.meta.arg
        state.algoParams[id] = {
          name: algoName,
          param: action.payload,
        }
      })
      .addCase(uploadImageFile.fulfilled, (state, action) => {
        const { elementId } = action.meta.arg
        const { tiffPath } = action.payload
        var idx = state.flowElements.findIndex((e) => e.id === elementId)
        const node = state.flowElements[idx]
        if (isNodeData(node) && node.data) {
          node.data = {
            ...node.data,
            type: NODE_DATA_TYPE_SET.DATA,
            path: tiffPath,
          }
        }
      })
      .addCase(runPipeline.pending, (state) => {
        state.runStatus = RUN_STATUS.RUNNING
      })
      .addCase(runPipeline.rejected, (state, action) => {
        state.runStatus = RUN_STATUS.FAILED
        state.runMessage = action.error.message
      })
      .addCase(runPipeline.fulfilled, (state, action) => {
        if (action.payload.message === 'success') {
          state.runStatus = RUN_STATUS.SUCCESS
        } else {
          state.runStatus = RUN_STATUS.FAILED
        }
        state.runMessage = action.payload.message
      })
      .addCase(stopPipeline, (state) => {
        state.runStatus = RUN_STATUS.STOPPED
      })
  },
})

export const { setFlowElements, updateParam, addFlowElement } =
  elementSlice.actions

export default elementSlice.reducer
