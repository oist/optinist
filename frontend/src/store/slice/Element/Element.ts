import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { Elements, Node, Position } from 'react-flow-renderer'
import {
  initialElements,
  INITIAL_ALGO_STYLE,
  INITIAL_DATA_STYLE,
} from 'const/flowchart'
import { NodeData, NODE_DATA_TYPE_SET } from 'const/NodeData'
import { isAlgoNodeData, isInputNodeData, isNodeData } from 'utils/ElementUtils'
import { Element, ELEMENT_SLICE_NAME } from './ElementType'
import { uploadImageFile } from '../ImageFile/ImageFileAction'
import { selectImageFile } from '../ImageFile/ImageFile'

const initialState: Element = {
  flowElements: initialElements,
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
        setFilePathFn(state, nodeId, tiffFilePath)
      })
      .addCase(uploadImageFile.pending, (state, action) => {
        const { nodeId } = action.meta.arg
        setFilePathFn(state, nodeId, undefined)
      })
      .addCase(selectImageFile, (state, action) => {
        const { nodeId, path } = action.payload
        setFilePathFn(state, nodeId, path)
      })
  },
})

function setFilePathFn(
  state: Element,
  nodeId: string,
  path: string | undefined,
) {
  const idx = state.flowElements.findIndex((e) => e.id === nodeId)
  const node = state.flowElements[idx]
  if (isNodeData(node) && node.data) {
    node.data = {
      ...node.data,
      type: NODE_DATA_TYPE_SET.DATA,
      path,
    }
  }
}

export const { setFlowElements, addFlowElement } = elementSlice.actions

export default elementSlice.reducer
