import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { Elements, Node, Position } from 'react-flow-renderer'
import {
  initialElements,
  INITIAL_ALGO_STYLE,
  INITIAL_DATA_STYLE,
} from 'const/flowchart'
import { NodeData, NODE_DATA_TYPE_SET } from 'const/NodeData'
import { isAlgoNodeData, isInputNodeData, isNodeData } from 'utils/ElementUtils'
import { runPipeline, stopPipeline } from './ElementAction'
import { Element, ELEMENT_SLICE_NAME, RUN_STATUS } from './ElementType'
import { uploadImageFile } from '../UploadImage/UploadImageAction'

const initialState: Element = {
  flowElements: initialElements,
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
      let node = action.payload
      if (isAlgoNodeData(node)) {
        node = {
          ...node,
          style: {
            ...node.style,
            ...INITIAL_ALGO_STYLE,
          },
          targetPosition: Position.Left,
          sourcePosition: Position.Right,
        }
      } else if (isInputNodeData(node)) {
        node = {
          ...node,
          style: {
            ...node.style,
            ...INITIAL_DATA_STYLE,
          },
          targetPosition: Position.Left,
          sourcePosition: Position.Right,
        }
      }
      state.flowElements.push(node)
    },
  },
  extraReducers: (builder) => {
    builder
      .addCase(uploadImageFile.fulfilled, (state, action) => {
        const { nodeId } = action.meta.arg
        const { tiffFilePath } = action.payload
        const idx = state.flowElements.findIndex((e) => e.id === nodeId)
        const node = state.flowElements[idx]
        if (isNodeData(node) && node.data) {
          node.data = {
            ...node.data,
            type: NODE_DATA_TYPE_SET.DATA,
            path: tiffFilePath,
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

export const { setFlowElements, addFlowElement } = elementSlice.actions

export default elementSlice.reducer
