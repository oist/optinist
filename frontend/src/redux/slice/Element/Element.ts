import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { Elements, Node, Position } from 'react-flow-renderer'
import { initialElements, INITIAL_ALGO_STYLE } from 'const/flowchart'
import { NodeData, NODE_DATA_TYPE_SET } from 'const/NodeData'
import { isAlgoNodeData, isNodeData } from 'utils/ElementUtils'
import { clickNode, runPipeline, stopPipeline } from './ElementAction'
import { Element, ELEMENT_SLICE_NAME, RUN_STATUS } from './ElementType'
import { uploadImageFile } from '../ImageIndex/ImageIndexAction'

const initialState: Element = {
  flowElements: initialElements,
  clickedNodeId: null,
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
      }
      state.flowElements.push(node)
    },
  },
  extraReducers: (builder) => {
    builder
      .addCase(clickNode, (state, action) => {
        state.clickedNodeId = action.payload.id
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

export const { setFlowElements, addFlowElement } = elementSlice.actions

export default elementSlice.reducer
