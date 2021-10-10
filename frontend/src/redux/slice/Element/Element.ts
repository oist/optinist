import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { Elements, Node } from 'react-flow-renderer'
import { initialElements, INITIAL_ALGO_ELEMENT_ID } from 'const/flowchart'
import { clickNode, getAlgoParams } from './ElementAction'
import { Element, ELEMENT_SLICE_NAME } from './ElementType'
import { NodeData } from 'const/NodeData'
import { uploadImageFile } from '../ImageIndex/ImageIndexAction'
import { isNodeData } from '../../../utils/ElementUtils'

const initialState: Element = {
  flowElements: initialElements,
  clickedNodeId: null,
  currentAlgoId: INITIAL_ALGO_ELEMENT_ID,
  algoParams: {},
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
        if (action.payload.type === 'algo') {
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
        const { elementId, fileName: path } = action.meta.arg
        var idx = state.flowElements.findIndex((e) => e.id === elementId)
        const node = state.flowElements[idx]
        if (isNodeData(node) && node.data) {
          node.data = {
            ...node.data,
            type: 'data',
            path: path,
          }
        }
      })
  },
})

export const { setFlowElements, updateParam, addFlowElement } =
  elementSlice.actions

export default elementSlice.reducer
